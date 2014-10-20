% COMPUTE_HAT compute the 'hat' variables on the glue (i.e., values going into
% the numerical flux calculation)
% [G] = compute_hat(G,B,Z)
%
% inputs:
%   G:  array of the glue data structures
%   B:  array of the block data structures
%   Z:  p-impedence
function [G] = compute_hat(G,Z)
  % projected to local glue
  for i = 1:G.nB
    jnd = G.jnd{i};
    G.qg{G.side(i)}.p(jnd,1) = G.Pf2g{i}*G.q{i}.p;
    G.qg{G.side(i)}.v(jnd,1) = G.Pf2g{i}*G.q{i}.v;
  end

  % project side 1 to glue
  G.qg{1}.pp = G.P12g*G.qg{1}.p;
  G.qg{1}.vp = G.P12g*G.qg{1}.v;

  % project side 2 to glue
  G.qg{2}.pp = G.P22g*G.qg{2}.p;
  G.qg{2}.vp = G.P22g*G.qg{2}.v;

  % compute the values
  [G.qg{1}.ph,G.qg{1}.vh,G.qg{2}.ph,G.qg{2}.vh]...
    = calculate_values(G.alpha,...
    G.qg{1}.pp,G.qg{1}.vp,...
    G.qg{2}.pp,G.qg{2}.vp,Z);

  % project glue to side 1
  G.qg{1}.ph = G.Pg21*G.qg{1}.ph;
  G.qg{1}.vh = G.Pg21*G.qg{1}.vh;

  % project glue to side 2
  G.qg{2}.ph = G.Pg22*G.qg{2}.ph;
  G.qg{2}.vh = G.Pg22*G.qg{2}.vh;
end

function [ph1,vh1,ph2,vh2] = calculate_values(a,p1,v1,p2,v2,Z)
  ph1 = (1/2)*(p1+p2) + a*(Z/2)*(v1+v2);

  ph2 = ph1;

  vh1 = (1/2)*(v1-v2) + a*(1/(2*Z))*(p1-p2);

  vh2 = -vh1;
end
