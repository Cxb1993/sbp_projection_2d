clear
format short g

add_matlab_bfam_paths

[Ba1,Ga1] = func_sbpdg(32,64,10,'straight_12.msh',5,0.2,false,eps,0.5,1);
[Ba0,Ga0] = func_sbpdg(32,64,10,'straight_12.msh',5,0.2,false,eps,0.5,0);

% Create the system
Aa0 = domain_matrix(Ba0,Ga0);
ea0 = eig(full(Aa0));
disp([max(real(ea0)),min(real(ea0))])

Aa1 = domain_matrix(Ba1,Ga1);
ea1 = eig(full(Aa1));
disp([max(real(ea1)),min(real(ea1))])

plot(real(ea1),imag(ea1),'b*',real(ea0),imag(ea0),'r*');
axis tight
legend('alpha = 0','alpha = 1')
