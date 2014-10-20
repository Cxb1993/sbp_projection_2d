% INIT_DG_BLOCK initialize an DG block
% [OP] = init_dg_block(rho,lam,p,mesh,curve,q,trans,alpha)
%
% inputs:
%   rho:   density
%   lam:   Lame's parameter
%   p:     DG order (polynomial order of p-1)
%   mesh:  mesh file
%   curve: How to remap the curved edges (see func_sbdg.m)
%   q:     initial condition structure (see modal_solution.m)
%   trans: how to move the grid points (see func_sbdg.m)
%   alpha: upwinding parameter
%
% output:
%   OP: block data structure for this block
function [OP] = init_dg_block(rho,lam,p,mesh,curve,q,trans,alpha)
  Globals2D;

  OP.sbp = false;
  OP.lam = lam;
  OP.rho = rho;
  OP.cp  = sqrt(OP.lam/OP.rho);
  OP.Z   = OP.rho*OP.cp;
  OP.p   = p;
  OP.alpha = alpha;

  N = p - 1;

  [Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(mesh);
  [VX,VY] = trans(VX,VY);
  StartUp2D;
  [k,f] = find(BCType==Out);
  curved = [];
  if(~isempty(k))
   outfaces = [k,f];
   curved = sort(unique(k));
   MakeCurved2D(outfaces, curve);
  end
  straight = setdiff(1:K, curved);
  BuildBCMaps2D

  intC = 2*(N+1);
  intG = 2*(N+1);

  OP.cubature = CubatureVolumeMesh2D(intC);
  OP.gauss    = GaussFaceMesh2D(intG);

  OP.gauss.f2ginterp = zeros(OP.gauss.NGauss*Nfaces,Nfp*Nfaces);

  % figure
  % PlotMesh2D()
  % figure

  % face 1
  faceR = r(Fmask(:,1));
  V1D = Vandermonde1D(N, faceR);
  I1 = InterpolationMatrix(N,OP.gauss.z,V1D);
  OP.gauss.f2ginterp(1:OP.gauss.NGauss, 1:Nfp) = I1;

  % face 2
  faceR = r(Fmask(:,2));
  V1D = Vandermonde1D(N, faceR);
  I2 = InterpolationMatrix(N,-OP.gauss.z,V1D);
  OP.gauss.f2ginterp(OP.gauss.NGauss+[1:OP.gauss.NGauss], Nfp+[1:Nfp]) = I2;

  % face 3
  faceS = s(Fmask(:,3));
  V1D = Vandermonde1D(N, faceS);
  I3 = InterpolationMatrix(N,-OP.gauss.z,V1D);
  OP.gauss.f2ginterp(2*OP.gauss.NGauss+[1:OP.gauss.NGauss], 2*Nfp+[1:Nfp]) = I3;

  rLGL = JacobiGL(0,0,N);
  rmin = abs(rLGL(1)-rLGL(2));
  dtscale = dtscale2D;

  OP.hmin = 0.75 * min(dtscale)*rmin;
  OP.dt = OP.hmin/OP.cp;

  OP.x1 = x;
  OP.x2 = y;
  % OP.MassMatrix = invV'*invV;
  OP.J = J;

  OP.q.v1 = q.v1(x,y,0);
  OP.q.v2 = q.v2(x,y,0);
  OP.q.pr = q.pr(x,y,0);

  OP.dq.v1 = 0*x;
  OP.dq.v2 = 0*x;
  OP.dq.pr = 0*x;

  % get glue grid

  vnum = [1 2; 2 3; 3 1];
  i = 1;
  for k = 1:K
    for f = 1:3
      if(BCType(k, f)==Out)
        va = EToV(k, vnum(f,1));
        vb = EToV(k, vnum(f,2));

        yave(i) = 0.5*(VY(va) + VY(vb));
        kave(i) = k;
        fave(i) = f;
        oidx = (k - 1) * Nfp * 3 + (f - 1) * Nfp + [1:Nfp];
        oave(i) = Fy(oidx(end)) > Fy(oidx(1));

        i = i + 1;
      end
    end
  end

  [ysort, I] = sort(yave);
  ksort = kave(I);
  fsort = fave(I);
  osort = oave(I);
  % [ysort' ksort' fsort' osort']

  % mapG = mapO;
  OP.g = zeros(length(I)+1,1);
  for i=1:length(I)
    k = ksort(i);
    f = fsort(i);
    o = osort(i);
    va = EToV(k, vnum(f,1));
    vb = EToV(k, vnum(f,2));

    if(VY(vb) < VY(va))
      OP.g(i)   = VY(vb);
      OP.g(i+1) = VY(va);
    else
      OP.g(i)   = VY(va);
      OP.g(i+1) = VY(vb);
    end

    oidx = (k - 1) * Nfp * 3 + (f - 1) * Nfp + [1:Nfp];
    goidx = (k-1)*OP.gauss.NGauss*3+(f-1)*OP.gauss.NGauss+[1:OP.gauss.NGauss];

    if(o)
      gidx = (i-1)*Nfp + [1:Nfp];
      ggidx = (i-1)*OP.gauss.NGauss + [1:OP.gauss.NGauss];
    else
      gidx = (i-1)*Nfp + [Nfp:-1:1];
      ggidx = (i-1)*OP.gauss.NGauss + [OP.gauss.NGauss:-1:1];
    end

    if(f==2 | f==3)
      if(f==2)
        warning('Untested fliplr for f==2');
      end
      ggidx = fliplr(ggidx);
    end

    mapG(gidx) = mapM(oidx);
    gmapG(ggidx) = OP.gauss.mapM(goidx);
  end
  OP.mapG = mapG';
  OP.gmapG = gmapG';
  OP.vmapG = vmapM(OP.mapG);

  %% size(OP.mapG)
  %% size(OP.gmapG)

  %%% Begin checks %%%
  %% fx = x(vmapM);
  %% fy = y(vmapM);

  %% fx = reshape(fx, Nfaces*Nfp, K);
  %% fy = reshape(fy, Nfaces*Nfp, K);

  %% gx = OP.gauss.f2ginterp*fx;
  %% gy = OP.gauss.f2ginterp*fy;

  %% norm(gauss.x(:)-gx(:))
  %% norm(gauss.y(:)-gy(:))

  %% fx_glue = zeros(size(nx));
  %% fy_glue = zeros(size(nx));

  %% fx_glue(OP.mapG) = fx(OP.mapG);
  %% fy_glue(OP.mapG) = fy(OP.mapG);

  %% gfx_glue = OP.gauss.f2ginterp*fx_glue;
  %% gfy_glue = OP.gauss.f2ginterp*fy_glue;

  %% norm([OP.gauss.x(gmapG)'-gfx_glue(gmapG)'])
  %% norm([OP.gauss.y(gmapG)'-gfy_glue(gmapG)'])
  %%% End checks %%%

  gK = length(OP.g) - 1;
  gN = OP.gauss.NGauss-1;
  gNp = gN + 1;
  OP.pg  = gNp;
  gr = OP.gauss.z;
  grV = Vandermonde1D(gN,gr);
  gm = [2./(2*[0:gN]+1)];

  pf2g = diag(sqrt(1./gm)) / grV;
  pg2f = grV * diag(sqrt(gm));

  OP.Pf2g = sparse(gK*gNp, gK*gNp);
  OP.Pg2f = sparse(gK*gNp, gK*gNp);

  OP.Pf2g = kron(speye(gK), pf2g);
  OP.Pg2f = kron(speye(gK), pg2f);

  OP.g = OP.g - OP.g(1);
  OP.g = 2*OP.g/OP.g(end)-1; % make it run from -1 to 1
  hg = kron(OP.g(2:end)-OP.g(1:end-1),ones(gNp,1));
  OP.gsJ = 2*OP.gauss.sJ(OP.gmapG)./hg;
  OP.gsJ = diag(sparse(OP.gsJ));
end
