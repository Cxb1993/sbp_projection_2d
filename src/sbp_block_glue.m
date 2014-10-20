% SBP_BLOCK_GLUE project the sbp solution to the glue
% q = sbp_block_glue(q,B,G,i)
%
% inputs:
%   q: solution data structure for the block
%   B: Block data structure
%   G: Glue data structure
%   i: segment of the glue we are updating
%
%  outputs
%   dq: glue solution data structure
function dq = sbp_block_glue(q,B,G,i)
  switch(G.edge(i))
    case('W')
      L  = B.Lw;
      n1 = B.nw1;
      n2 = B.nw2;
    case('E')
      L  = B.Le;
      n1 = B.ne1;
      n2 = B.ne2;
    case('S')
      L  = B.Ls;
      n1 = B.ns1;
      n2 = B.ns2;
    case('N')
      L  = B.Ln;
      n1 = B.nn1;
      n2 = B.nn2;
  end

  % original
  dq.p = sqrt(G.sJ{i})*L'*q.pr;
  dq.v = sqrt(G.sJ{i})*(n1*L'*q.v1 + n2*L'*q.v2);

  if(G.flip(i))
    dq.p = flipud(dq.p);
    dq.v = flipud(dq.v);
  end
end
