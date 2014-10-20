function [Pc2f, Pf2c] = make_projection_g2g_h_gen(N, xi_c)
% [Pc2f, Pf2c] = make_projection_g2g_h_gen(N, xi_c)
%
% Generate projection operators to go from the f grid to the c grid both with
% the same order of accuracy N:
%
%    c      f
%    o      o
%    |      |
%    |      |
%    |      |
%    |      o
%    |      |
%    |      |
%    |      |
%    |      o
%    |      |
%    |      |
%    |      |
%    o      o
%
%  xi_c(i,j) is the xi component of the i'th vertex of the j'th cell for the
%            left mesh (assumed to be ordered from smallest to largest).
%
%  It is assumed that xi_c(1,1) and xi_c(2,end) are the components of the
%  f cell.

[r, ~] = Gaussian_quad(N);

K = size(xi_c, 2);

VX = [xi_c(1,:), xi_c(2,end)];
fa = xi_c(1,1);
fb = xi_c(2,end);


VX = 2/(fb-fa) * VX + 1 - 2*fb/(fb-fa);

va = (1:K)'; vb = (2:K+1)';
x = ones(N+1,1)*VX(va) + 0.5*(r+1)*(VX(vb)-VX(va));

Jc = (VX(vb) - VX(va))/2;
Jf = 1;

Pc2f = zeros(length(x(:)), length(r(:)));

m = (2./(2*(0:N)+1));
Vr  = Legendre_Vandermonde(r, N);
sq_invm_invVr = diag(sqrt(1./m)) / Vr;
for k = 1:K
  Pc2f((k-1)*(N+1)+1:k*(N+1), :) = ...
      sq_invm_invVr * Legendre_Vandermonde(x(:,k),N) * diag(sqrt(m));
end

Pf2c = (1/Jf) * diag(1./m) * Pc2f' * kron(diag(Jc), diag(m));
