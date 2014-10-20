% GLUE_TO_SBP_BLOCK_FLUX take the flux back to the grid from the glue for an SBP
% block
% [B] = glue_to_sbp_block_flux(B,G,i)
%
% inputs:
%   B: Block data structure
%   G: Glue data structure
%   i: the ith block of the glue
%
% outputs:
%   B: updated block data structure
function [B] = glue_to_sbp_block_flux(B,G,i)
  jnd = G.jnd{i};

  s = G.side(i);

  pr = sqrt(G.sJI{i})*G.Pg2f{i}*G.qg{s}.ph(jnd,1)...
     + (1/2)*sqrt(G.sJI{i})*G.Perr{i}*G.q{i}.p;
  vn = sqrt(G.sJI{i})*G.Pg2f{i}*G.qg{s}.vh(jnd,1)...
     + (1/2)*sqrt(G.sJI{i})*G.Perr{i}*G.q{i}.v;

  if(G.flip(i))
    pr = flipud(pr);
    vn = flipud(vn);
  end

  switch(G.edge(i))
    case('W')
      B.f.w.pr = pr;
      B.f.w.vn = vn;
    case('E')
      B.f.e.pr = pr;
      B.f.e.vn = vn;
    case('S')
      B.f.s.pr = pr;
      B.f.s.vn = vn;
    case('N')
      B.f.n.pr = pr;
      B.f.n.vn = vn;
  end
end
