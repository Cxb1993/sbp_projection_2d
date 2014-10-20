function [Vxi_g,Pl2g,Pg2l,Pr2g,Pg2r] = make_projection_g2g_hr(N, Vxi_l, Vxi_r)
% [Vxi_g, Pl2g, Pg2l, Pr2g, Pg2r] = make_projection_g2g_hr(N, xi_l, xi_r)
%
% Generate projection operators to go from the {l,r} grid to g and from g to
% {l,r} same order of accuracy N:
%
%    l   g   r
%    o   o   o
%    |   |   |
%    o   o   |
%    |   |   |
%    |   o   o
%    |   |   |
%    |   |   |
%    o   o   |
%    |   |   |
%    |   o   o
%    |   |   |
%    |   |   |
%    |   |   |
%    o   o   o
%
% Here
%
%  Vxi_l(i,j) is the right mesh of xi cell edges
%  Vxi_r(i,j) is the left mesh of xi cell edges
%
%  Vxi_g(i,j) if the glue mesh xi cell edges

% absolute tolerance equality
% isequalAbs = @(x,y,tol) ( abs(x-y) <= tol );

% relative tolerance equality
isequalRel = @(x,y,tol) ( abs(x-y) <= ( tol*max(abs(x),abs(y)) + eps) );

Np = N+1;

Vxi_l = Vxi_l(:)';
Vxi_r = Vxi_r(:)';

Nv_l = length(Vxi_l);
Nv_r = length(Vxi_r);

K_l = Nv_l - 1;
K_r = Nv_r - 1;

tol = 100*eps;

% Check to make sure the first and last grid points lineup
assert(isequalRel(Vxi_l(1),   Vxi_r(1),   tol));
assert(isequalRel(Vxi_l(end), Vxi_r(end), tol));

Vxi_g = sort([Vxi_l Vxi_r]);
dVxi_g = arrayfun(isequalRel, Vxi_g(1:end-1), Vxi_g(2:end), ...
                  tol*ones(1,length(Vxi_g)-1));
Vxi_g(dVxi_g) = [];

xi_g = [Vxi_g(1:end-1); Vxi_g(2:end)];
K_g = size(xi_g, 2);

Np_l = K_l * Np;
Np_r = K_r * Np;
Np_g = K_g * Np;

%% Loop through the elements
Pg2l = sparse(Np_l, Np_g);
Pl2g = sparse(Np_g, Np_l);

for k_l = 1:K_l
  xia_l = Vxi_l(k_l);
  xib_l = Vxi_l(k_l+1);

  k_g = find(all([Vxi_g < xib_l - tol;  Vxi_g >= xia_l - tol]));

  [Pc2f, Pf2c] = make_projection_g2g_h_gen(N, xi_g(:, k_g));

  idx_l = (k_l   -1) * Np + 1;
  idx_g = (k_g(1)-1) * Np + 1;

  Pg2l(idx_l:idx_l+Np-1, idx_g:idx_g+Np*length(k_g)-1) = Pf2c;
  Pl2g(idx_g:idx_g+Np*length(k_g)-1, idx_l:idx_l+Np-1) = Pc2f;
end

Pg2r = sparse(Np_r, Np_g);
Pr2g = sparse(Np_g, Np_r);

for k_r = 1:K_r
  xia_r = Vxi_r(k_r);
  xib_r = Vxi_r(k_r+1);

  k_g = find(all([Vxi_g < xib_r - tol;  Vxi_g >= xia_r - tol]));

  [Pc2f, Pf2c] = make_projection_g2g_h_gen(N, xi_g(:, k_g));

  idx_r = (k_r   -1) * Np + 1;
  idx_g = (k_g(1)-1) * Np + 1;

  Pg2r(idx_r:idx_r+Np-1, idx_g:idx_g+Np*length(k_g)-1) = Pf2c;
  Pr2g(idx_g:idx_g+Np*length(k_g)-1, idx_r:idx_r+Np-1) = Pc2f;
end

return
