function B = hat_bc_sbp_block(B)
  % last argument is the upwind parameter. 0 is central, 1 is updwind
  B.f.w = hat_bc_sbp_block_side(B,B.q,B.Lw,B.nw1,B.nw2,B.alpha);
  B.f.e = hat_bc_sbp_block_side(B,B.q,B.Le,B.ne1,B.ne2,B.alpha);
  B.f.s = hat_bc_sbp_block_side(B,B.q,B.Ls,B.ns1,B.ns2,B.alpha);
  B.f.n = hat_bc_sbp_block_side(B,B.q,B.Ln,B.nn1,B.nn2,B.alpha);
end

function f = hat_bc_sbp_block_side(OP,q,L,n1,n2,a)
  Z = sqrt(OP.rho*OP.lam);
  [pr,vn] = hat_bc(L'*q.v1,L'*q.v2,L'*q.pr,n1,n2,Z,a);
  f.pr = pr;
  f.vn = vn;
end

function [pr,vn] = hat_bc(v1,v2,pr,n1,n2,Z,a)
  vn  = n1*v1 + n2*v2;
  vn  = vn + a*pr/Z;
  pr  = 0;
end
