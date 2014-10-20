function [table] = latex_table(filename,suf)

load(filename{1});
eval(['N = N',suf,';']);
error = zeros(length(N),length(filename));
order = zeros(length(N),length(filename));
for k = 1:length(filename)
  load(filename{k});
  eval(['m = length(N',suf,');'])
  eval(['error (1:m,k) =  er',suf,';']);
  eval(['order  (2:m,k) =  p',suf,';']);
end

color = {'red','blu','grn','org','prp'};
table = '';
for i = 1:length(N)
  table = sprintf('%s$%3d$',table,N(i));
  for k = 1:length(filename)
    MAG = floor(log10(error(i,k)));
    VAL = round(error(i,k)/10^(MAG-2))/100;
    table = sprintf('%s & \\color{%s} $%1.1f \\times 10^{%3d}$',table,color{k},VAL,MAG);
    if i == 1
      table = sprintf('%s $\\phantom{(0.0)}$',table);
    else
      table = sprintf('%s $         (%1.1f) $',table,order(i,k));
    end
  end
  if i < length(N)
    table = sprintf('%s\\\\\n',table);
  end
end
