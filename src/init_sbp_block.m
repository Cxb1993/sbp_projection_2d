% INIT_SBP_BLOCK initialize an SBP block
% [OP] = init_sbp_block(rho,lam,N1,N2,p,grid,q,alpha)
%
% inputs:
%   rho:  density
%   lam:  Lame's parameter
%   N1:   N in the first dimension
%   N2:   N in the second dimension
%   p:    FD interior order
%   grid: grid data structure (see transfinite.m)
%   q:    solution data structure (see modal_solution.m)
%   alpha: upwinding parameter
%
% output:
%   OP: block data structure for this block

function [OP] = init_sbp_block(rho,lam,N1,N2,p,grid,q,alpha)
  OP.sbp = true;
  OP.lam = lam;
  OP.rho = rho;
  OP.cp  = sqrt(OP.lam/OP.rho);
  OP.Z   = OP.rho*OP.cp;
  OP.N1  = N1;
  OP.N2  = N2;
  OP.p   = p;
  OP.pg  = p;
  OP.alpha = alpha;

  % Set up reference grid
  r1 = linspace(-1,1,OP.N1+1)';OP.h1 = r1(2)-r1(1);
  r2 = linspace(-1,1,OP.N2+1)';OP.h2 = r2(2)-r1(1);
  [r1,r2] = meshgrid(r1,r2);
  OP.r1 = r1(:);
  OP.r2 = r2(:);
  OP.x1 = grid.x1(OP.r1,OP.r2);
  OP.x2 = grid.x2(OP.r1,OP.r2);

  % plot(OP.x1 ,OP.x2 ,'.k');
  % hold on

  % evalute the metric relations
  OP.x1_1 = grid.x1_1(OP.r1,OP.r2);
  OP.x1_2 = grid.x1_2(OP.r1,OP.r2);
  OP.x2_1 = grid.x2_1(OP.r1,OP.r2);
  OP.x2_2 = grid.x2_2(OP.r1,OP.r2);

  % set up the Jacobian
  OP.J = OP.x1_1.*OP.x2_2-OP.x1_2.*OP.x2_1;
  OP.JI = 1./OP.J;

  % set up the inverse metrics
  OP.r1_1 =  OP.JI.*OP.x2_2;
  OP.r1_2 = -OP.JI.*OP.x1_2;
  OP.r2_1 = -OP.JI.*OP.x2_1;
  OP.r2_2 =  OP.JI.*OP.x1_1;

  % make these all sparse diagonal matrices
  OP.J    = diag(sparse(   OP.J));
  OP.JI   = diag(sparse(  OP.JI));
  OP.x1_1 = diag(sparse(OP.x1_1));
  OP.x1_2 = diag(sparse(OP.x1_2));
  OP.x2_1 = diag(sparse(OP.x2_1));
  OP.x2_2 = diag(sparse(OP.x2_2));
  OP.r1_1 = diag(sparse(OP.r1_1));
  OP.r1_2 = diag(sparse(OP.r1_2));
  OP.r2_1 = diag(sparse(OP.r2_1));
  OP.r2_2 = diag(sparse(OP.r2_2));

  % set up the difference operators
  [OP.D1,OP.Hi1] = diagonal_sbp(OP.p, OP.N1);
  OP.D1  = kron( OP.D1,speye(OP.N2+1))/OP.h1;
  OP.H1  = kron(diag(sparse(1./diag(OP.Hi1))),speye(OP.N2+1))*OP.h1;
  OP.Hi1 = kron(OP.Hi1,speye(OP.N2+1))/OP.h1;

  [OP.D2,OP.Hi2] = diagonal_sbp(OP.p, OP.N2);
  OP.D2  = kron(speye(OP.N1+1), OP.D2)/OP.h2;
  OP.H2  = kron(speye(OP.N1+1),diag(sparse(1./diag(OP.Hi2))))*OP.h2;
  OP.Hi2 = kron(speye(OP.N1+1),OP.Hi2)/OP.h2;

  % Lift operators
  OP.Lw = kron(sparse(      1,1,1,OP.N1+1,1),speye(OP.N2+1));
  OP.Le = kron(sparse(OP.N1+1,1,1,OP.N1+1,1),speye(OP.N2+1));
  OP.Ls = kron(speye(OP.N1+1),sparse(      1,1,1,OP.N2+1,1));
  OP.Ln = kron(speye(OP.N1+1),sparse(OP.N2+1,1,1,OP.N2+1,1));

  % Setup the fluxes
  OP.f1_v1_pr = -OP.J*OP.r1_1;
  OP.f1_v2_pr = -OP.J*OP.r1_2;
  OP.f1_pr_v1 = -OP.lam*OP.J*OP.r1_1;
  OP.f1_pr_v2 = -OP.lam*OP.J*OP.r1_2;

  OP.f2_v1_pr = -OP.J*OP.r2_1;
  OP.f2_v2_pr = -OP.J*OP.r2_2;
  OP.f2_pr_v1 = -OP.lam*OP.J*OP.r2_1;
  OP.f2_pr_v2 = -OP.lam*OP.J*OP.r2_2;

  % Setup the fluxes on each boundary
  OP.fw_v1_pr = OP.Lw'*OP.f1_v1_pr;
  OP.fw_v2_pr = OP.Lw'*OP.f1_v2_pr;
  OP.fw_pr_v1 = OP.Lw'*OP.f1_pr_v1;
  OP.fw_pr_v2 = OP.Lw'*OP.f1_pr_v2;

  OP.fe_v1_pr = OP.Le'*OP.f1_v1_pr;
  OP.fe_v2_pr = OP.Le'*OP.f1_v2_pr;
  OP.fe_pr_v1 = OP.Le'*OP.f1_pr_v1;
  OP.fe_pr_v2 = OP.Le'*OP.f1_pr_v2;

  OP.fs_v1_pr = OP.Ls'*OP.f2_v1_pr;
  OP.fs_v2_pr = OP.Ls'*OP.f2_v2_pr;
  OP.fs_pr_v1 = OP.Ls'*OP.f2_pr_v1;
  OP.fs_pr_v2 = OP.Ls'*OP.f2_pr_v2;

  OP.fn_v1_pr = OP.Ln'*OP.f2_v1_pr;
  OP.fn_v2_pr = OP.Ln'*OP.f2_v2_pr;
  OP.fn_pr_v1 = OP.Ln'*OP.f2_pr_v1;
  OP.fn_pr_v2 = OP.Ln'*OP.f2_pr_v2;

  % compute the normals
  tmp = sqrt(diag(OP.x1_2.^2+OP.x2_2.^2));
  OP.mw  = diag(OP.Lw'*tmp);
  OP.nw1 =-diag(OP.Lw'*(diag(OP.x2_2)./tmp));
  OP.nw2 = diag(OP.Lw'*(diag(OP.x1_2)./tmp));

  OP.me  = diag(OP.Le'*tmp);
  OP.ne1 = diag(OP.Le'*(diag(OP.x2_2)./tmp));
  OP.ne2 =-diag(OP.Le'*(diag(OP.x1_2)./tmp));

  tmp = sqrt(diag(OP.x1_1.^2+OP.x2_1.^2));
  OP.ms  = diag(OP.Ls'*tmp);
  OP.ns1 = diag(OP.Ls'*(diag(OP.x2_1)./tmp));
  OP.ns2 =-diag(OP.Ls'*(diag(OP.x1_1)./tmp));

  OP.mn  = diag(OP.Ln'*tmp);
  OP.nn1 =-diag(OP.Ln'*(diag(OP.x2_1)./tmp));
  OP.nn2 = diag(OP.Ln'*(diag(OP.x1_1)./tmp));

  OP.hmin = full(min(...
  min(min(sqrt(diag(OP.x1_1.^2+OP.x2_1.^2))))/OP.N1,...
  min(min(sqrt(diag(OP.x1_2.^2+OP.x2_2.^2))))/OP.N2));
  OP.dt = OP.hmin/OP.cp;

  % plot(OP.Lw'*OP.x1,OP.Lw'*OP.x2,'k*-')
  % hold on
  % quiver(OP.Lw'*OP.x1,OP.Lw'*OP.x2,diag(OP.nw1),diag(OP.nw2),'r')

  % plot(OP.Le'*OP.x1,OP.Le'*OP.x2,'k*-')
  % hold on
  % quiver(OP.Le'*OP.x1,OP.Le'*OP.x2,diag(OP.ne1),diag(OP.ne2),'r')

  % plot(OP.Ls'*OP.x1,OP.Ls'*OP.x2,'k*-')
  % hold on
  % quiver(OP.Ls'*OP.x1,OP.Ls'*OP.x2,diag(OP.ns1),diag(OP.ns2),'r')

  % plot(OP.Ln'*OP.x1,OP.Ln'*OP.x2,'k*-')
  % hold on
  % quiver(OP.Ln'*OP.x1,OP.Ln'*OP.x2,diag(OP.nn1),diag(OP.nn2),'r')

  % hold off
  % axis equal
  % pause

  OP.q.v1 = q.v1(OP.x1,OP.x2,0);
  OP.q.v2 = q.v2(OP.x1,OP.x2,0);
  OP.q.pr = q.pr(OP.x1,OP.x2,0);

  OP.dq.v1 = 0*OP.x1;
  OP.dq.v2 = 0*OP.x1;
  OP.dq.pr = 0*OP.x1;
end
