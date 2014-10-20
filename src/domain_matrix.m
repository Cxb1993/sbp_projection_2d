% DOMAIN_MATRIX create the matrix operator for the full system
% [A] = domain_matrix(B,G)
%
% inputs:
%   B: cell array of blocks data structure
%   G: cell array of glue data structure
%
% output:
%   A: matrix operator
function [A] = domain_matrix(B,G)

  % find number of dof
  N = 0;
  for b = 1:length(B)
    N = N + length(B{b}.q.v1);
    B{b}.q.v1 = zeros(size(B{b}.q.v1));
    B{b}.q.v2 = zeros(size(B{b}.q.v2));
    B{b}.q.pr = zeros(size(B{b}.q.pr));
  end

  c = 0;
  A = sparse(3*N,3*N);

  for b = 1:length(B)
    for k = 1:3
      for i = 1:length(B{b}.q.v1(:))
        c = c+1;

        % set single entry to 1
        if(k == 1)
          B{b}.q.v1(i) = 1;
        elseif(k == 2)
          B{b}.q.v2(i) = 1;
        elseif(k == 3)
          B{b}.q.pr(i) = 1;
        end

        % compute RHS
        [B,G] = domain_rhs(B,G,0);

        M = 0;
        for t = 1:length(B)
          % v1
          A(c,M+(1:length(B{t}.q.v1(:)))) = B{t}.rhs.v1(:);
          M = M + length(B{t}.q.v1(:));

          % v2
          A(c,M+(1:length(B{t}.q.v2(:)))) = B{t}.rhs.v2(:);
          M = M + length(B{t}.q.v2(:));

          % pr
          A(c,M+(1:length(B{t}.q.pr(:)))) = B{t}.rhs.pr(:);
          M = M + length(B{t}.q.pr(:));
        end

        % set single entry to 1
        if(k == 1)
          B{b}.q.v1(i) = 0;
        elseif(k == 2)
          B{b}.q.v2(i) = 0;
        elseif(k == 3)
          B{b}.q.pr(i) = 0;
        end
      end
    end
  end
end
