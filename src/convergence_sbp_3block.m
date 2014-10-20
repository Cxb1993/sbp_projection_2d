% convergence_sbp_3block driver to three block SBP convergence test
% convergence_sbp_3block(p1,p2,p3,AMP,grid_type,cfl,alpha)
%
% inputs:
%   p1:         SBP accuracy on block 1
%   p2:         SBP accuracy on block 2
%   p3:         SBP accuracy on block 3
%   AMP:        amplitude of the sine wave curve of domain
%   grid_type: 'n' nested refinement with ratio of 2
%              'u' unnested refinement with one more point than ratio of 2
%   cfl:       time step fudge factor
%   alpha:     upwinding parameter

function convergence_sbp_3block(p1,p2,p3,AMP,grid_type,cfl,alpha)
  filename = ['data/sbp_3block_tend_10',...
             '_p1_',num2str(p1),'_p2_',num2str(p2),'_p3_',num2str(p3),...
             '_AMP_',num2str(AMP),'_gt_',grid_type,...
             '_cfl_',num2str(cfl),'_alpha_',num2str(alpha),...
             '.mat'];

  format short g
  NMIN = 32;

  N = [];
  if exist(filename,'file') == 2
    load(filename);
    switch grid_type
      case 'n'
        N = N_n;
      case 'u'
        N = N_u;
    end
  end

  % loop through the number of levels we want to handle
  for i = (length(N)+1):5
    switch grid_type
      case 'n'
        disp(':: NESTED :: 3BLOCK ::')
        N_n(i,1) = NMIN*2^i;
        disp([N_n(i,1),i,p1,p2,p3])

        [~,~,energy,error,~] =...
          func_3block(N_n(i,1)/2,N_n(i,1),p1,p2,p3,2,1,0,0,AMP,false,1,cfl,alpha);
        en_n(i,1) = energy(end);
        den_n(i,1) = energy(end)-energy(1);
        er_n(i,1) = error(end);
        p_n  = log(er_n(2:end,:)./er_n(1:end-1,:))./...
               log(N_n(1:end-1,:)./N_n(2:end,:));
        disp([den_n(:,1)';er_n(:,1)';0,p_n(:,1)']);drawnow

        eval(['save ',filename,' N_n den_n en_n er_n p_n']);

      case 'u'
        disp(':: UNNESTED :: 3BLOCK ::')
        N_u(i,1) = NMIN*2^i;
        disp([N_u(i,1),i,p1,p2,p3])

        [~,~,energy,error,~] =...
          func_3block(N_u(i,1)/2,N_u(i,1),p1,p2,p3,2,1+1/N_u(i,1),0,0,...
                      AMP,false,1,cfl,alpha);
        en_u(i,1) = energy(end);
        den_u(i,1) = energy(end)-energy(1);
        er_u(i,1) = error(end);
        p_u  = log(er_u(2:end,:)./er_u(1:end-1,:))./...
               log(N_u(1:end-1,:)./N_u(2:end,:));
        disp([den_u(:,1)';er_u(:,1)';0,p_u(:,1)']);drawnow

        eval(['save ',filename,' N_u den_u en_u er_u p_u']);
    end
  end
end
