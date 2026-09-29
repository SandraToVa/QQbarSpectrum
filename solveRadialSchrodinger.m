function [E, W, x, info] = solveRadialSchrodinger(Vfun, m, nstates, opts)
% [E, W, x, info] = solveRadialSchrodinger(Vfun, m, nstates, opts)
% Lowest eigenvalues and wave functions of the coupled-channel radial problem
%
%     -u''(r) + V(r) u(r) = m*E u(r),   u(0) = 0,  u(inf) = 0
%
% using a Lagrange-Laguerre mesh (see radialMesh). The whole spectrum comes
% from a single symmetric eigenvalue problem, so there is no shooting,
% no bracketing and no overflow; degenerate levels are returned correctly.
%
% INPUT:  Vfun    : handle, Vfun(r) -> n x n x numel(r) potential matrix
%                   (includes centrifugal terms and the factor m)
%         m       : quark mass; the eigenvalues are returned as E = lambda/m
%         nstates : number of states to return (lowest ones)
%         opts    : (optional) struct with fields
%                   .N     number of mesh points per channel (default 100)
%                   .rmax  last mesh point in GeV^-1           (default 25)
%                   .x     points where W is evaluated (default linspace(0,22,1000))
%                   .H     precomputed Hamiltonian (skips building it)
%                   .mesh  mesh used to build opts.H
% OUTPUT: E    : 1 x nstates energies (GeV)
%         W    : n x numel(x) x nstates radial wave functions u_c(r),
%                normalized to sum_c int |u_c|^2 dr = 1, sign fixed so that
%                the largest component is positive
%         x    : 1 x numel(x) points where W is evaluated
%         info : struct with the mesh, the eigenvectors on the mesh (C), the
%                channel weights of each state (weights, n x nstates) and
%                the number of channels n

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'N'),    opts.N = 100;  end
if ~isfield(opts, 'rmax'), opts.rmax = 25; end
if ~isfield(opts, 'x'),    opts.x = linspace(0, 22, 1000); end

if isfield(opts, 'H') && ~isempty(opts.H)
    H = opts.H;
    mesh = opts.mesh;
    n = size(H,1) / mesh.N;
else
    mesh = radialMesh(opts.N, opts.rmax);
    [H, n] = radialHamiltonian(Vfun, mesh);
end
N = mesh.N;

nstates = min(nstates, n*N);
[C, lam] = eig(H, 'vector');
[lam, idx] = sort(lam);
lam = lam(1:nstates);
C = C(:, idx(1:nstates));
E = lam(:)' / m;

% Channel weights (probability in each channel), C is normalized to 1
weights = zeros(n, nstates);
for s = 1:nstates
    Cs = reshape(C(:,s), n, N);
    weights(:,s) = sum(Cs.^2, 2);
end

info.mesh = mesh;
info.n = n;
info.C = C;
info.weights = weights;
info.H = H;

if nargout < 2
    return;
end

% Wave functions evaluated at the requested points
x = opts.x(:)';
G = lagrangeInterp(mesh, x);             % numel(x) x N
scale = 1 ./ sqrt(mesh.h * mesh.lambda); % u(r_i) = c_i / sqrt(h*lambda_i)
W = zeros(n, numel(x), nstates);
for s = 1:nstates
    Cs = reshape(C(:,s), n, N) .* scale';  % values u_c(r_i)
    Ws = Cs * G';                          % n x numel(x)
    [~, imax] = max(abs(Ws(:)));
    W(:,:,s) = sign(Ws(imax)) * Ws;
end
end


function G = lagrangeInterp(mesh, r)
% G(k,j) = value at r(k) of the Lagrange function that is 1 at r_j and 0 at
% the other mesh points (regularized Laguerre: it behaves as r near 0)
xs = r(:) / mesh.h;
xj = mesh.x';
LN = scaledLaguerreN(mesh.N, xs);          % L_N(x) e^{-x/2}
D = xs - xj;
G = (xs ./ xj) .* LN ./ (D .* mesh.dLN');
% points that coincide with a mesh point
[kk, jj] = find(abs(D) < 1e-12);
for q = 1:numel(kk)
    G(kk(q),:) = 0;
    G(kk(q),jj(q)) = 1;
end
end


function L = scaledLaguerreN(N, x)
Lm1 = zeros(size(x));
L   = exp(-x/2);
for n = 0:N-1
    Lp1 = ((2*n+1-x).*L - n*Lm1) / (n+1);
    Lm1 = L;
    L   = Lp1;
end
end
