% CONVERGENCE_SBPDG driver to two block SBP convergence test
% convergence_sbpdg(p1,p2,AMP,cfl,alpha)
%
% inputs:
%   p1:        SBP accuracy
%   p2:        DG polynomial order
%   AMP:       amplitude of the sine wave curve of domain
%   cfl:       time step fudge factor
%   alpha:     upwinding parameter
function convergence_sbpdg(p1,p2,AMP,cfl,alpha)
  filename = ['data/sbpdg_tend_10',...
             '_p1_',num2str(p1),'_p2_',num2str(p2),...
             '_AMP_',num2str(AMP),...
             '_cfl_',num2str(cfl),'_alpha_',num2str(alpha),...
             '.mat'];


  format short g

  p2 = p2+1;
  NMIN = 32;

  O = ceil(NMIN/p2);

  N = [];
  if exist(filename,'file') == 2
    load(filename);
  end

  for i = (length(N)+1):5
    N(i,1) = NMIN*2^i;
    mesh = ['straight_',num2str(O*2^i),'.msh'];

    disp([N(i,1),i,p1,p2]);
    disp(mesh);
    drawnow

    [~,~,energy,error,~] =...
      func_sbpdg(N(i,1)/2,N(i,1),p1,mesh,p2,AMP,false,1,cfl,alpha);
    en(i,1) = energy(end);
    den(i,1) = energy(end)-energy(1);
    er(i,1) = error(end);

    p  = log(er(2:end,:)./er(1:end-1,:))./log(N(1:end-1,:)./N(2:end,:));
    disp([den(:,1)';er(:,1)';0,p(:,1)']);drawnow

    eval(['save ',filename,' N den en er p']);
  end
end
