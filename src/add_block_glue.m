% ADD_BLOCK_GLUE add a block to the glue grid
% [OP] = add_block_glue(OP,side,block,edge,flip,r,sJ,sc)
%
% inputs:
%   OP:    glue grid structure
%   side:  which side of the glue grid (1 or 2)
%   block: which block number
%   edge:  which edge of the block (N,S,E,W)
%   flip:  flip the data order?
%   r:     glue grid
%   sJ:    surface Jacobian
%   sc:    scale factor (optional)
%
% outputs:
%   OP: updated glue grid
function [OP] = add_block_glue(OP,side,block,edge,flip,r,sJ,sc)

  OP.nB           = OP.nB+1;
  OP.block(OP.nB) = block;
  OP.edge(OP.nB)  = edge;
  OP.side(OP.nB)  = side;
  OP.r{OP.nB}     = r;
  OP.sJ{OP.nB}    = sJ;
  OP.sJI{OP.nB}   = diag(sparse(1./diag(sJ)));
  OP.flip(OP.nB)  = flip;
  OP.sc{OP.nB}    = [];
  if(nargin > 7)
    OP.sc{OP.nB} = sc;
  end
end
