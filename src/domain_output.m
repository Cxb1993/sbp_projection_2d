% DOMAIN_OUTPUT write the domain to VTK files
%
% inputs:
%   B: cell array of blocks
%   s: time step number
function domain_output(B,s)
  %%% figure(2)
  if(nargin == 1)
    for i = 1:length(B)
      if(B{i}.sbp)
        SBPWriteVTK2D(sprintf('mesh_%04d.vtk',i), ...
          B{i}.N1, B{i}.N2, B{i}.x1, B{i}.x2, ...
          {});
      else
        % Write VTK Output
        WriteVTK2D(sprintf('mesh_%04d.vtk',i), 1, {});
      end
    end
  else
    for i = 1:length(B)
      if(B{i}.sbp)
        SBPWriteVTK2D(sprintf('sol_%04d_%04d.vtk',i,s), ...
          B{i}.N1, B{i}.N2, B{i}.x1, B{i}.x2, ...
          {'v1','v2','pr'}, B{i}.q.v1, B{i}.q.v2, B{i}.q.pr);
        %%% pcolored(...
        %%%   reshape(B{i}.x1,B{i}.N1+1,B{i}.N2+1),...
        %%%   reshape(B{i}.x2,B{i}.N1+1,B{i}.N2+1),...
        %%%   reshape(B{i}.q.pr,B{i}.N1+1,B{i}.N2+1))
        %%% hold on
      else
        % Write VTK Output
        WriteVTK2D(sprintf('sol_%04d_%04d.vtk',i,s), 2*B{i}.p-1, ...
          {'v1','v2','pr'}, B{i}.q.v1, B{i}.q.v2, B{i}.q.pr);
      end
    end
  end
  %%% hold off
  %%% axis tight
  %%% caxis([-1 1])
  %%% drawnow
  %%% figure(1)
end
