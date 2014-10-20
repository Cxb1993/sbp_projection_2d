% FUNC_2BLOCK set up and run a two block projection case
% [B,G,energy,err,time] =
%                    func_2block(N1,N2,p1,p2,p3,s12,s22,s13,s23,AMP,doplot,tend,cfl,alpha)
%
% inputs: all optional with default in [ ]
%   N1: N1+1 grid points perpendicular to interface on block 1 [16]
%   N2: N2+1 grid points along the interface on block 1 [16]
%   p1: SBP order on block 1
%   p2: SBP order on block 2
%   p3: SBP order on block 3
%   s12: refinement ratio to get block 2 grids point in dimension 1 [1].
%   s13: refinement ratio to get block 2 grids point in dimension 2 [1].
%   s22: refinement ratio to get block 3 grids point in dimension 1 [1].
%   s23: refinement ratio to get block 3 grids point in dimension 2 [1].
%   AMP:    size of the sine curve, should be between -1 and 1 [0.1]
%   doplot: output done every doplot timesteps [10]
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
function [B,G,energy,err,time] = func_3block(varargin)
  % set defaults for optional inputs
  %          N1,N2,p1,p2,p3,s12,s22,s13,s23, AMP,doplot,tend,cfl,alpha
  optargs = {64,64, 8, 8, 8,  1,  1,  1,  1, 0.1,    10,   1,0.5,    1};

  MID = 0;

  % skip any new inputs if they are empty
  newVals = cellfun(@(x) ~isempty(x), varargin);

  % now put these defaults into the valuesToUse cell array,
  % and overwrite the ones specified in varargin.
  optargs(newVals) = varargin(newVals);

  % Place optional args in memorable variable names
  [N1,N2,p1,p2,p3,s12,s22,s13,s23,AMP,doplot,tend,cfl,alpha] = optargs{:};


  % material properties
  [q,~,rho,lam,w] = modal_solution();

  %%%%%%%%%%%%%
  %% BLOCK 1 %%
  %%%%%%%%%%%%%

  % coordinate transform;
  NW = [-1, 1];
  SW = [-1,-1];
  NE = [ 0, 1];
  SE = [ 0,-1];

  lx  = @(r,i,P1,P2) (r+1)*P1(i)/2 + (1-r)*P2(i)/2;
  ldx = @(r,i,P1,P2)       P1(i)/2 -       P2(i)/2;

  LW.x1   = @(r2)  lx(r2,1,NW,SW);
  LW.x2   = @(r2)  lx(r2,2,NW,SW);
  LW.dx1  = @(r2) ldx(r2,1,NW,SW);
  LW.dx2  = @(r2) ldx(r2,2,NW,SW);

  LE.x1   = @(r2)  lx(r2,1,NE,SE)+sin(pi*(r2+1))*AMP;
  LE.x2   = @(r2)  lx(r2,2,NE,SE);
  LE.dx1  = @(r2) ldx(r2,1,NE,SE)+pi*cos(pi*(r2+1))*AMP;
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

  % numerics parameters
  N1 = B{1}.N1;
  N2 = B{1}.N2;
  if(s12 > 0)
    N1 = s12*N1;
  end
  if(s22 > 0)
    N2 = s22*N2;
  end

  % coordinate transform;
  NW = [0,  1];
  SW = [AMP*sin(pi*(MID+1)),MID];
  NE = [1,  1];
  SE = [1,MID];

  lx  = @(r,i,P1,P2) (r+1)*P1(i)/2 + (1-r)*P2(i)/2;
  ldx = @(r,i,P1,P2)       P1(i)/2 -       P2(i)/2;

  LW.x1   = @(r2)  lx(r2,1,NW,SW)+sin(pi*(MID*(1-r2)/2+(1+r2)/2 + 1))*AMP;
  LW.x2   = @(r2)  lx(r2,2,NW,SW);
  LW.dx1  = @(r2) ldx(r2,1,NW,SW)+pi*(1-MID)/2*...
                                  cos(pi*(MID*(1-r2)/2+(1+r2)/2 + 1))*AMP;
  LW.dx2  = @(r2) ldx(r2,2,NW,SW);

  LE.x1   = @(r2)  lx(r2,1,NE,SE);
  LE.x2   = @(r2)  lx(r2,2,NE,SE);
  LE.dx1  = @(r2) ldx(r2,1,NE,SE);
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

  [B{2}] = init_sbp_block(rho,lam,N1,N2,p2,grid,q,alpha);
  clear N1 N2 grid LE LN LS LW NE NW SE SW ldx lx

  %%%%%%%%%%%%%
  %% BLOCK 3 %%
  %%%%%%%%%%%%%

  % numerics parameters
  N1 = B{2}.N1;
  N2 = B{2}.N2;
  if(s13 > 0)
    N1 = s13*N1;
  end
  if(s23 > 0)
    N2 = s23*N2;
  end

  % coordinate transform;
  NW = [AMP*sin(pi*(MID+1)),MID];
  SW = [0, -1];
  NE = [1,MID];
  SE = [1, -1];

  lx  = @(r,i,P1,P2) (r+1)*P1(i)/2 + (1-r)*P2(i)/2;
  ldx = @(r,i,P1,P2)       P1(i)/2 -       P2(i)/2;

  LW.x1   = @(r2)  lx(r2,1,NW,SW)+sin(pi*(-(1-r2)/2+MID*(1+r2)/2 + 1))*AMP;
  LW.x2   = @(r2)  lx(r2,2,NW,SW);
  LW.dx1  = @(r2) ldx(r2,1,NW,SW)+pi*(MID+1)/2*cos(pi*(-(1-r2)/2+MID*(1+r2)/2 + 1))*AMP;
  LW.dx2  = @(r2) ldx(r2,2,NW,SW);

  LE.x1   = @(r2)  lx(r2,1,NE,SE);
  LE.x2   = @(r2)  lx(r2,2,NE,SE);
  LE.dx1  = @(r2) ldx(r2,1,NE,SE);
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

  [B{3}] = init_sbp_block(rho,lam,N1,N2,p3,grid,q,alpha);
  clear N1 N2 grid LE LN LS LW NE NW SE SW ldx lx

  % Setup glue grid
  %% Glue Grid 1
  [G{1}] = init_glue_grid(alpha);
  [G{1}] = add_block_glue(G{1},1,1,'E',false,...
                          B{1}.Le'*B{1}.r2          ,B{1}.me);
  [G{1}] = add_block_glue(G{1},2,3,'W',false,...
                          -(1-B{3}.Lw'*B{3}.r2)/2 + MID*(1+B{3}.Lw'*B{3}.r2)/2,...
                          B{3}.mw/((MID+1)/2));
  [G{1}] = add_block_glue(G{1},2,2,'W',false,...
                          MID*(1-B{2}.Lw'*B{2}.r2)/2 + (1+B{2}.Lw'*B{2}.r2)/2,...
                          B{2}.mw/((1-MID)/2));
  [G{1}] = alloc_glue_grid(G{1},B,false);

  %% Glue Grid 2
  [G{2}] = init_glue_grid(alpha);
  [G{2}] = add_block_glue(G{2},1,2,'S',false,B{2}.Ls'*B{2}.r1,B{2}.ms);
  [G{2}] = add_block_glue(G{2},2,3,'N',false,B{3}.Ln'*B{3}.r1,B{3}.mn);
  [G{2}] = alloc_glue_grid(G{2},B,s13==0);

  % Time stepping
  t = 0;
  dt = cfl*min([B{1}.dt,B{2}.dt,B{3}.dt]);
  tend = tend/w;
  nsteps = ceil(tend/dt);
  dt = tend/nsteps;

  [B,G,energy,err,time] = domain_update(B,G,dt,nsteps,q,doplot);
end
