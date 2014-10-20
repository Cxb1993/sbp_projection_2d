% MAKE_PROJECTION make the projection operators for the SBP to glue
% [Pf2g,Pg2f,M,H] = make_projection(N,ORDER)
%
% inputs:
%   N:     finite difference N
%   ORDER: finite difference order
%
% outputs:
%   Pf2g: projection from SBP to glue
%   Pg2f: projection from glue to SBP
%   M:    mass matrix for the glue
%   H:    H matrix for the SBP operator
function [Pf2g,Pg2f,M,H] = make_projection(N,ORDER)

  if(ORDER == 2)
    load optimal_2_q.mat
  elseif(ORDER == 4)
    load optimal_4_q.mat
  elseif(ORDER == 6)
    load optimal_6_q.mat
  elseif(ORDER == 8)
    load optimal_8_q.mat
  elseif(ORDER == 10)
    load optimal_10_q.mat
  else
    error('order must be 2, 4, 6, 8,10...');
  end

  p = ORDER-1;

  s = r-(m/2-1);

  %% Get the SBP operator
  [~,HI] = diagonal_sbp(ORDER, N);
  H = diag(sparse(1./diag(HI)));

  xf_b = 0:s+m-2; % fd grid
  xg_b = 0:s+m-2; % glue grid
  [~, ~, Mr_b] = glue_pieces(p, xf_b, xg_b);

  I = kron((s+1:N-s+1)',ones(m*(p+1),1));
  J = kron(ones(N+1-2*s,1),(1:m*(p+1))')+(p+1)*(I-1-m/2);
  qi = kron(ones(N+1-2*s,1),q(1:m*(p+1)));
  Pg2f = sparse(I,J,qi,N+1,N*(p+1));

  Qb = reshape(q(m*(p+1)+1:end),r*(p+1),s)';
  Pg2f(1:s,1:r*(p+1)) = Qb;

  Qb = reshape(diag(2*mod(1:p+1,2)-1)*flipud(reshape(rot90(Qb,2)',...
                                                     p+1,r*s)),r*(p+1),s)';
  Pg2f(end+1-s:end,end+1-r*(p+1):end) = Qb;

  M = kron(speye(N),Mr_b);
  Pf2g = kron(speye(N),diag(sparse(1./diag(Mr_b))))*Pg2f'*H;
end
