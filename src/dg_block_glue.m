% DG_BLOCK_GLUE project the db solution to the glue
% q = dg_block_glue(q,B,G,i)
%
% inputs:
%   q: solution data structure for the block
%   B: Block data structure
%   G: Glue data structure
%   i: segment of the glue we are updating
%
%  outputs
%   dq: glue solution data structure
function dq = dg_block_glue(q,B,G,i)

  Globals2D;

  % Store modal coefficients for local glue where elements
  % are stored with increasing y

  % interpolate volume to boundary cubature
  gv1 = B.gauss.interp*q.v1;
  gv2 = B.gauss.interp*q.v2;
  gpr = B.gauss.interp*q.pr;

  pr = sqrt(G.sJ{i})*gpr(B.gmapG);
  v1 = sqrt(G.sJ{i})*gv1(B.gmapG);
  v2 = sqrt(G.sJ{i})*gv2(B.gmapG);

  % p
  dq.p = pr;

  % \vec{n} \cdot \vec{v}
  dq.v = (v1.*B.gauss.nx(B.gmapG)+v2.*B.gauss.ny(B.gmapG));
end
