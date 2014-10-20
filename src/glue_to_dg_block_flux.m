% GLUE_TO_DG_BLOCK_FLUX take the flux back to the grid from the glue for an DG
% block
% [B] = glue_to_dg_block_flux(B,G,i)
%
% inputs:
%   B: Block data structure
%   G: Glue data structure
%   i: the ith block of the glue
%
% outputs:
%   B: updated block data structure
function [B] = glue_to_dg_block_flux(B,G,i)
  s = G.side(i);

  % p* v*
  B.pr_glue = sqrt(G.sJI{i})*G.Pg2f{i}*G.qg{s}.ph;
  B.vn_glue = sqrt(G.sJI{i})*G.Pg2f{i}*G.qg{s}.vh;
end
