% ALLOC_GLUE_GRID allocate a glue grid
% function [G] = alloc_glue_grid(G,B,idProj)
%
% inputs:
%   G:      the glue grid data structure to allocate
%   B:      the blocks data structure
%   idProj: boolean for identity projection
%
% outputs:
%   G:  newly allocated glue grid
function [G] = alloc_glue_grid(G,B,idProj)
  G.p    = 0;  % order of the glue grid
  G.r1   = []; % glue grid on side one
  G.r2   = []; % glue gird on side two
  G.ind  = {}; % where i indexes into its side glue <<< not used anymore
  G.jnd  = {}; % where i indexes into the global glue

  % figure out the glue grid order
  if(idProj)
    G.p = 0; % zero means identity projection
  else
    for i = 1:G.nB
      G.p = max(G.p,B{G.block(i)}.pg);
    end
  end

  % Create each projection matrix
  for i = 1:G.nB
    pb = B{G.block(i)}.pg;
    N = length(G.r{i})-1;
    if(idProj)
      G.Pg2f{i} = speye(N+1);
      G.Pf2g{i} = speye(N+1);
      G.Perr{i} = sparse(N+1,N+1);
    else
      if B{G.block(i)}.sbp
        [Pf2g, Pg2f]  = make_projection(N,pb);
      else
        % Given a set of nodal value make modal values
        % implement me
        [Pf2g, Pg2f]  = make_projection_dg(B{G.block(i)});
      end
      if ~isempty(G.sc{i})
        sc  = sqrt(kron(diag(sparse(   G.sc{i})),speye(pb)));
        scI = sqrt(kron(diag(sparse(1./G.sc{i})),speye(pb)));
        Pf2g = sc*Pf2g;
        Pg2f = Pg2f*scI;
      end
      [Pa2b, Pb2a] = make_projection_g2g_p(pb-1,G.p-1);
      G.Pf2g{i} = kron(speye(N),Pa2b)*Pf2g;
      G.Pg2f{i} = Pg2f*kron(speye(N),Pb2a);
      if B{G.block(i)}.sbp
        G.Perr{i} = speye(N+1) - G.Pg2f{i}*G.Pf2g{i};
      end
    end
  end

  % Create grid on either side
  j1 = 0;
  j2 = 0;
  for i = 1:G.nB
    if(G.side(i) == 1)
      G.ind{i} = length(G.r1)+1;
      if(isempty(G.r1))
        G.r1 = shiftdim(G.r{i});
      else
        G.r1 = [G.r1;shiftdim(G.r{i})];
      end
      G.ind{i} = G.ind{i}:length(G.r1);
      G.jnd{i} = j1+(1:(length(G.ind{i})-1)*max(G.p,1));
      j1 = G.jnd{i}(end);
    else
      G.ind{i} = length(G.r2)+1;
      if(isempty(G.r2))
        G.r2 = shiftdim(G.r{i});
      else
        G.r2 = [G.r2;shiftdim(G.r{i})];
      end
      G.ind{i} = G.ind{i}:length(G.r2);
      G.jnd{i} = j2+(1:(length(G.ind{i})-1)*max(G.p,1));
      j2 = G.jnd{i}(end);
    end
  end

  % Create side 1 to side 2 projection
  if(idProj)
    G.rg = G.r1;
    N = length(G.rg)-1;
    G.P12g = speye(N+1);
    G.Pg21 = speye(N+1);
    G.P22g = speye(N+1);
    G.Pg22 = speye(N+1);
    G.jnd = G.ind;
  else
    [G.rg, G.P12g, G.Pg21, G.P22g, G.Pg22] = ...
       make_projection_g2g_hr(G.p-1, unique(G.r1), unique(G.r2));
  end
  G = rmfield(G,'ind');
end
