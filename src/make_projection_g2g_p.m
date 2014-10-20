function [Pa2b, Pb2a] = make_projection_g2g_p(Na, Nb)
% [Pa2b, Pb2a] = make_projection_g2g_p(Na, Nb)
%
% Generate projection operators to go from the a grid with order Na to the
% aligned b grid with order Nb:
%
%    a      b
%    o      o
%    |      |
%    |      |
%    |      |
%    |      |
%    |      |
%    |      |
%    |      |
%    |      |
%    |      |
%    |      |
%    |      |
%    o      o
%

Pa2b = eye(Nb+1, Na+1);
Pb2a = eye(Na+1, Nb+1);

return
