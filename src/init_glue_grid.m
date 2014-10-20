% INIT_GLUE_GRID initialize a glue grid
% function [OP] = init_glue_grid(a)
%
% inputs:
% a: alpha to use in the glue grid flux. Any a >= 0 is stable with a = 0 being
%    central and a = 1 being upwind
%
% outputs:
%   OP: new glue grid
function [OP] = init_glue_grid(a)
  OP.nB = 0;     % number of blocks
  OP.alpha = a;
end
