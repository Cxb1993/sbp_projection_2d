% FUNC_SBPDG set up and run a two block sbp-dg test
% [B,G,energy,err,time] = func_sbpdg(N1,N2,p1,mesh,p2,AMP,doplot,tend,cfl,alpha)
%
% inputs: all optional with default in [ ]
%   N1:     N1+1 grid points perpendicular to interface on side 1 [16]
%   N2:     N2+1 grid points along the interface on side 1 [16]
%   p1:     SBP order
%   mesh:   DG mesh file to use
%   p2:     DG order (polynomial order is p2-1)
%   AMP:    size of the sine curve, should be between -1 and 1 [0.1]
%   doplot: output done every doplot timesteps [1]
%           doplot = 0 results in no output.
%           if doplot > 0 user must create directory vtk
%   tend:   number of oscilations to do in simulations [1]
%   cfl:    time step fudge factor [0.5]
%   alpha:  upwinding parameter
%
% outputs:
%   B:      blocks data structure
%   G:      glue grid data structure
%   energy: energy time history
%   err:    error time history
%   time:   time points
function [B,G,energy,err,time,berror] = func_sbpdg(varargin)
  % set defaults for optional inputs
  %          N1,N2,p1,                              mesh, p2, AMP, doplot, tend, cfl, alpha
  optargs = {16,16, 4, 'straight_12.msh',  4, 0.1,      1,   10, 0.5,     1};

  % skip any new inputs if they are empty
  newVals = cellfun(@(x) ~isempty(x), varargin);

  % now put these defaults into the valuesToUse cell array,
  % and overwrite the ones specified in varargin.
  optargs(newVals) = varargin(newVals);

  % Place optional args in memorable variable names
  [N1,N2,p1,mesh,p2,AMP,doplot,tend,cfl,alpha] = optargs{:};

  % Setup the initial conditions
  [q,~,rho,lam,w] = modal_solution();

  lx  = @(r,i,P1,P2) (1+r)*P1(i)/2 + (1-r)*P2(i)/2;
  ldx = @(r,i,P1,P2)       P1(i)/2 -       P2(i)/2;

  %%%%%%%%%%%%%
  %% BLOCK 1 %%
  %%%%%%%%%%%%%

  % coordinate transform;
  if(strncmp(mesh,'straight_',9))
    NW = [-1, 1];
    SW = [-1,-1];
    NE = [ 0, 1];
    SE = [ 0,-1];
    AMP_W = 0;
    AMP_E = AMP;
    coupling = 'E';
    dg_curve = @(x,y) deal(-AMP*sin(pi*y),y);

    LW.x1   = @(r2)  lx(r2,1,NW+[1,0],SW+[1,0])+AMP_E*sin(pi*(r2+1));
    LW.x2   = @(r2)  lx(r2,2,NW+[1,0],SW+[1,0]);
    LW.dx1  = @(r2) ldx(r2,1,NW+[1,0],SW+[1,0])+pi*AMP_E*cos(pi*(r2+1));
    LW.dx2  = @(r2) ldx(r2,2,NW+[1,0],SW+[1,0]);

    LE.x1   = @(r2)  lx(r2,1,NE+[1,0],SE+[1,0])+AMP_W*sin(pi*(r2+1));
    LE.x2   = @(r2)  lx(r2,2,NE+[1,0],SE+[1,0]);
    LE.dx1  = @(r2) ldx(r2,1,NE+[1,0],SE+[1,0])+pi*AMP_W*cos(pi*(r2+1));
    LE.dx2  = @(r2) ldx(r2,2,NE+[1,0],SE+[1,0]);

    LS.x1   = @(r1)  lx(r1,1,SE+[1,0],SW+[1,0]);
    LS.x2   = @(r1)  lx(r1,2,SE+[1,0],SW+[1,0]);
    LS.dx1  = @(r1) ldx(r1,1,SE+[1,0],SW+[1,0]);
    LS.dx2  = @(r1) ldx(r1,2,SE+[1,0],SW+[1,0]);

    LN.x1   = @(r1)  lx(r1,1,NE+[1,0],NW+[1,0]);
    LN.x2   = @(r1)  lx(r1,2,NE+[1,0],NW+[1,0]);
    LN.dx1  = @(r1) ldx(r1,1,NE+[1,0],NW+[1,0]);
    LN.dx2  = @(r1) ldx(r1,2,NE+[1,0],NW+[1,0]);

    dg_grid = transfinite(LW,LE,LS,LN);
    dg_trans = @(x,y) deal(dg_grid.x1(2*x-1, y),dg_grid.x2(2*x, y));
  else
    error(['invalid mesh ',mesh]);
  end

  LW.x1   = @(r2)  lx(r2,1,NW,SW)+AMP_W*sin(pi*(r2+1));
  LW.x2   = @(r2)  lx(r2,2,NW,SW);
  LW.dx1  = @(r2) ldx(r2,1,NW,SW)+pi*AMP_W*cos(pi*(r2+1));
  LW.dx2  = @(r2) ldx(r2,2,NW,SW);

  LE.x1   = @(r2)  lx(r2,1,NE,SE)+AMP_E*sin(pi*(r2+1));
  LE.x2   = @(r2)  lx(r2,2,NE,SE);
  LE.dx1  = @(r2) ldx(r2,1,NE,SE)+pi*AMP_E*cos(pi*(r2+1));
  LE.dx2  = @(r2) ldx(r2,2,NE,SE);

  LS.x1   = @(r1)  lx(r1,1,SE,SW);
  LS.x2   = @(r1)  lx(r1,2,SE,SW);
  LS.dx1  = @(r1) ldx(r1,1,SE,SW);
  LS.dx2  = @(r1) ldx(r1,2,SE,SW);

  LN.x1   = @(r1)  lx(r1,1,NE,NW);
  LN.x2   = @(r1)  lx(r1,2,NE,NW);
  LN.dx1  = @(r1) ldx(r1,1,NE,NW);
  LN.dx2  = @(r1) ldx(r1,2,NE,NW);

  grid = transfinite(LW,LE,LS,LN);

  [B{1}] = init_sbp_block(rho,lam,N1,N2,p1,grid,q,alpha);
  clear N1 N2 grid LE LN LS LW NE NW SE SW ldx lx


  %%%%%%%%%%%%%
  %% BLOCK 2 %%
  %%%%%%%%%%%%%

  [B{2}] = init_dg_block(rho,lam,p2,mesh,dg_curve,q,dg_trans,alpha);

  % Setup glue grid
  [G{1}] = init_glue_grid(alpha);
  if coupling == 'W'
    [G{1}] = add_block_glue(G{1},1,1,'W',false,B{1}.Lw'*B{1}.r2,B{1}.mw);
  elseif coupling == 'E'
    [G{1}] = add_block_glue(G{1},1,1,'E',false,B{1}.Le'*B{1}.r2,B{1}.me);
  else
    assert(0)
  end

  [G{1}] = add_block_glue(G{1},2,2,'X',false,B{2}.g,B{2}.gsJ);
  [G{1}] = alloc_glue_grid(G{1},B,false);

  % Time stepping
  dt = cfl*min(B{1}.dt,B{2}.dt);
  tend = tend/w;
  nsteps = ceil(tend/dt);
  dt = tend/nsteps;
  [B,G,energy,err,time,berror] = domain_update(B,G,dt,nsteps,q,doplot);
end
