function mesh = radialMesh(N, rmax)
% mesh = radialMesh(N, rmax)
% Regularized Lagrange-Laguerre mesh for the radial Schrodinger equation
%
%     -u''(r) + V(r) u(r) = lambda u(r),   u(0) = 0,  u(inf) = 0
%
% The basis functions behave as r near the origin (so u(0)=0 is built in)
% and decay exponentially at large r. With the Gauss quadrature
% approximation the kinetic matrix is known analytically and the potential
% is diagonal (just evaluated at the mesh points), so the Hamiltonian is a
% plain symmetric matrix. Accuracy converges exponentially with N.
% See D. Baye, "The Lagrange-mesh method", Phys. Rep. 565 (2015) 1, Sec. 3.4.5.
%
% INPUT:  N    : number of mesh points (per channel)
%         rmax : position of the last mesh point (sets the scale h)
% OUTPUT: mesh.r      : mesh points (N x 1), r = h*x
%         mesh.x      : zeros of the Laguerre polynomial L_N (N x 1)
%         mesh.h      : scale factor
%         mesh.lambda : Gauss-Laguerre weights including exp(x)
%         mesh.T      : kinetic matrix for -d^2/dr^2 (N x N)
%         mesh.dLN    : L_N'(x_i)*exp(-x_i/2), used for interpolation

% Zeros of L_N: eigenvalues of the symmetric Jacobi (Golub-Welsch) matrix
k = (1:N)';
J = diag(2*k-1) + diag(k(1:end-1),1) + diag(k(1:end-1),-1);
x = sort(eig(J));

% L_{N-1}(x_i)*exp(-x_i/2) from the scaled three-term recurrence
% (the exp(-x/2) factor keeps everything finite for large x)
[~, LNm1] = scaledLaguerre(N, x);
dLN = -N .* LNm1 ./ x;          % L_N'(x_i) e^{-x_i/2}, since L_N(x_i)=0
lambda = 1 ./ (x .* dLN.^2);    % Gauss weights times exp(x_i)

% Kinetic matrix -d^2/dx^2 (regularized Laguerre mesh, Gauss approximation)
[XI, XJ] = ndgrid(x, x);
sgn = (-1).^(k - k');
T = sgn .* (XI + XJ) ./ (sqrt(XI.*XJ) .* (XI - XJ).^2);
T(1:N+1:end) = -(x.^2 - 2*(2*N+1)*x - 4) ./ (12*x.^2);

h = rmax / x(end);

mesh.N = N;
mesh.x = x;
mesh.h = h;
mesh.r = h*x;
mesh.lambda = lambda;
mesh.dLN = dLN;
mesh.T = T / h^2;
end


function [LN, LNm1] = scaledLaguerre(N, x)
% Returns L_N(x)*exp(-x/2) and L_{N-1}(x)*exp(-x/2)
Lm1 = zeros(size(x));
L   = exp(-x/2);
for n = 0:N-1
    Lp1 = ((2*n+1-x).*L - n*Lm1) / (n+1);
    Lm1 = L;
    L   = Lp1;
end
LN = L;
LNm1 = Lm1;
end
