% DOMAIN_RHS compute the domain right hand side
% [B,G] = domain_rhs(B,G,t)
%
% inputs:
%  B: block data structure
%  G: glue data structure
%  t: time of the domain
%
% outputs:
%  B: updated block data structure
%  G: updated glue data structure
function [B,G] = domain_rhs(B,G,t)

  % set the fluxes
  for i = 1:length(B)
    if B{i}.sbp
      B{i} = hat_bc_sbp_block(B{i});
    end
  end

  % do glue grid stuff
  for i = 1:length(G)
    % Send grid to glue
    for j = 1:G{i}.nB
      if B{G{i}.block(j)}.sbp
        G{i}.q{j} = sbp_block_glue(B{G{i}.block(j)}.q,B{G{i}.block(j)},G{i},j);
      else
        G{i}.q{j} = dg_block_glue(B{G{i}.block(j)}.q,B{G{i}.block(j)},G{i},j);
      end
    end

    % compute flux on the glue
    G{i} = compute_hat(G{i},B{G{i}.block(1)}.Z);

    % Send glue to grid
    for j = 1:G{i}.nB
      if B{G{i}.block(j)}.sbp
        B{G{i}.block(j)} = ...
          glue_to_sbp_block_flux(B{G{i}.block(j)},G{i},j);
      else
        B{G{i}.block(j)} = ...
          glue_to_dg_block_flux(B{G{i}.block(j)},G{i},j);
      end
    end

    % just for debugging get rid of these once we are done with them
    G{i} = rmfield(G{i},'qg');
  end

  % do the update
  for i = 1:length(B)
    if B{i}.sbp
      B{i}.rhs  = ...
        rhs_sbp_block(B{i},B{i}.q,B{i}.f.w,B{i}.f.e,B{i}.f.s,B{i}.f.n);
    else
      B{i}.rhs  = rhs_dg_block(B{i});
    end
  end
end
