% RHS_DG_BLOCK compute the RHS of an DG block
% [dq] = rhs_dg_block(OP)
%
% inputs:
%   OP:      block data structure
%
% outputs:
%   dq: rate data structure
function [dq] = rhs_dg_block(OP)

Globals2D;

assert(OP.sbp == false)

v1 = OP.q.v1;
v2 = OP.q.v2;
pr = OP.q.pr;

rho = OP.rho;
lam = OP.lam;

c = OP.cubature;
g = OP.gauss;


 mapG = OP.mapG;
gmapG = OP.gmapG;
vmapG = OP.vmapG;

rhoinv = 1/rho;
% cp = sqrt(lam/rho);
zp = sqrt(rho*lam);

cdv1dr = c.Dr*v1; cdv1ds = c.Ds*v1;
cdv2dr = c.Dr*v2; cdv2ds = c.Ds*v2;

cdv1dx = c.rx.*cdv1dr + c.sx.*cdv1ds;
cdv2dy = c.ry.*cdv2dr + c.sy.*cdv2ds;

dq.pr = -c.V'*(c.W.*lam.*(cdv1dx + cdv2dy));

% Weak Form
cpr = c.V*pr;
dq.v1 = rhoinv.*(c.Dr'*(c.W.*c.rx.*cpr) + c.Ds'*(c.W.*c.sx.*cpr));
dq.v2 = rhoinv.*(c.Dr'*(c.W.*c.ry.*cpr) + c.Ds'*(c.W.*c.sy.*cpr));

% For now we are working with constant material
glam_m = lam;

% interpolate to integration nodes
gv1 = g.interp*v1;
gv2 = g.interp*v2;
gpr = g.interp*pr;

gdv1 = gv1(g.mapM)-gv1(g.mapP);
gdv2 = gv2(g.mapM)-gv2(g.mapP);
gdpr = gpr(g.mapM)-gpr(g.mapP);

% add the pressures on both sides of interface
gapr = gpr(g.mapM)+gpr(g.mapP);

gdv1(g.mapB) = 0;
gdv2(g.mapB) = 0;
gdpr(g.mapB) = 2*gpr(g.mapB);
gapr(g.mapB) = 0;

% Compute flux with alpha term
gndotdv = g.nx.*gdv1 + g.ny.*gdv2;
galpha = (-OP.alpha.*gdpr./zp + gndotdv)/2;

% Weak Flux Term
gbeta  = -(gapr + zp.*OP.alpha.*gndotdv)/2;

gfluxv1 = gbeta.*g.nx;
gfluxv2 = gbeta.*g.ny;
gfluxpr = galpha.*glam_m;

if isfield(OP, 'pr_glue')
  % assume constant material
  flam_m = glam_m;

  gfluxv1(gmapG) = -OP.pr_glue .* g.nx(gmapG);
  gfluxv2(gmapG) = -OP.pr_glue .* g.ny(gmapG);
  gfluxpr(gmapG) =...
         (g.nx(gmapG).*gv1(gmapG)+g.ny(gmapG).*gv2(gmapG)-OP.vn_glue).*flam_m;
end

dq.v1 = dq.v1 + g.interp'*(g.W.*gfluxv1);
dq.v2 = dq.v2 + g.interp'*(g.W.*gfluxv2);
dq.pr = dq.pr + g.interp'*(g.W.*gfluxpr);

for k=1:K
  mmCHOL = cub.mmCHOL(:,:,k);
  dq.v1(:,k) = mmCHOL\(mmCHOL'\dq.v1(:,k));
  dq.v2(:,k) = mmCHOL\(mmCHOL'\dq.v2(:,k));
  dq.pr(:,k) = mmCHOL\(mmCHOL'\dq.pr(:,k));
end
