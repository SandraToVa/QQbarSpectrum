function [H, n] = radialHamiltonian(Vfun, mesh)
% [H, n] = radialHamiltonian(Vfun, mesh)
% Builds the (n*N x n*N) symmetric Hamiltonian matrix of the coupled-channel
% radial problem  -u'' + V(r) u = lambda u  on a Lagrange-Laguerre mesh.
%
% INPUT:  Vfun : function handle, Vfun(r) with r a row vector of length K
%                returns the n x n x K potential matrix (it must already
%                contain the centrifugal terms and the factor m)
%         mesh : output of radialMesh
% OUTPUT: H    : Hamiltonian. Ordering of the unknowns: the channel index
%                runs fastest, i.e. unknown (i-1)*n + c = channel c at r_i
%         n    : number of channels

Vr = Vfun(mesh.r');
n = size(Vr, 1);
N = mesh.N;
if size(Vr, 3) ~= N
    error('radialHamiltonian:size', ...
        'The potential matrix returned %d points, expected %d.', size(Vr,3), N);
end
if ~all(isfinite(Vr(:)))
    error('radialHamiltonian:nonfinite', 'The potential matrix contains NaN/Inf values.');
end
if ~any(Vr(:))
    error('radialHamiltonian:zero', ...
        'The potential matrix is identically zero: check the mix/hf flags (that branch may not be implemented).');
end

% Block-diagonal potential: block i is V(r_i)
Vblk = cell(1, N);
for i = 1:N
    Vblk{i} = Vr(:,:,i);
end
H = kron(mesh.T, eye(n)) + blkdiag(Vblk{:});
H = (H + H') / 2;
end
