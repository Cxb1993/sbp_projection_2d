function [plot] = latex_plots(filename,suf,lst)

load(filename{1});
eval(['N = N',suf,';']);
error = zeros(length(N),length(filename));
for k = 1:length(filename)
  load(filename{k});
  eval(['m = length(N',suf,');'])
  eval(['error (1:m,k) =  er',suf,';']);
end

color = {'red','blu','grn','org','prp'};
mark  = {'*','triangle*','diamond*','square*','pentagon*'};
plot = '';
for i = 1:length(filename)
  if nargin < 3
    plot = sprintf('%s\\addplot[color=%s,mark=%s] plot coordinates {%c\n',plot,color{i},mark{i},'%');
  else
    plot = sprintf('%s\\addplot[color=%s,mark=%s,%s] plot coordinates {%c\n',plot,color{i},mark{i},lst,'%');
  end
  for k = 1:length(N);
    plot = sprintf('%s(%d,%e)\n',plot,N(k),error(k,i));
  end
  plot = sprintf('%s};\n',plot);
end
