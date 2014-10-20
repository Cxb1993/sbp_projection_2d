function [q,cp,rho,lam,w1] = modal_solution()
  % material properties
  rho  = 1;
  lam  = 1;
  cp   = sqrt(lam/rho);

  a1  = 1;
  kx1 = 0.5;
  ky1 = 0.5;
  w1  = (cp/2)*sqrt(kx1^2+ky1^2);

  a2  = 1;
  kx2 = 1;
  ky2 = 1;
  w2  = (cp/2)*sqrt(kx2^2+ky2^2);
  q.pr = @(x1,x2,t)...
         a1*cos(2*pi*w1*t)*cos(kx1*pi*x1).*cos(ky1*pi*x2)...
       + a2*cos(2*pi*w2*t)*sin(kx2*pi*x1).*sin(ky2*pi*x2);

  q.v1 = @(x1,x2,t)...
          ( a1*kx1/(rho*2*w1))*sin(2*pi*w1*t)*sin(kx1*pi*x1).*cos(ky1*pi*x2)...
        + (-a2*kx2/(rho*2*w2))*sin(2*pi*w2*t)*cos(kx2*pi*x1).*sin(ky2*pi*x2);

  q.v2 = @(x1,x2,t)...
          ( a1*ky1/(rho*2*w1))*sin(2*pi*w1*t)*cos(kx1*pi*x1).*sin(ky1*pi*x2)...
        + (-a2*ky2/(rho*2*w2))*sin(2*pi*w2*t)*sin(kx2*pi*x1).*cos(ky2*pi*x2);
end
