% DOMAIN_UPDATE update a domain
% [B,G,energy,error,time,berror] = domain_update(B,G,dt,nsteps,q_e,nplot)
%
% inputs:
%   B:      blocks data structure
%   G:      glue data structure
%   dt:     time step size
%   nsteps: number of time steps
%   q_e:    exact solution (see modal_solution.m)
%   nplot:  plot every nplot time steps
%
% outputs:
%   B:      update blocks data states
%   G:      update glue data structure
%   energy: energy time history
%   error:  error time history
%   time:   time step history
%   berror: individual blocks error history

function [B,G,energy,error,time,berror] = domain_update(B,G,dt,nsteps,q_e,nplot)
  if nplot
    domain_output(B);
  end

  for s = 1:nsteps
    t = (s-1)*dt;
    [B,G] = domain_update_RK(B,G,t,dt);
    t = s*dt;

    if(mod(s,nplot)==0 && nplot)
      domain_output(B,s);
    end

    % do the update
    energy(s) = 0;
    error(s) = 0;
    berror(s) = 0;
    for i = 1:length(B)
      B{i}.e.pr  = B{i}.q.pr - q_e.pr(B{i}.x1,B{i}.x2,t);
      B{i}.e.v1  = B{i}.q.v1 - q_e.v1(B{i}.x1,B{i}.x2,t);
      B{i}.e.v2  = B{i}.q.v2 - q_e.v2(B{i}.x1,B{i}.x2,t);

      if B{i}.sbp
        energy(s) = energy(s)...
          + (B{i}.rho/2)*(B{i}.q.v1'*B{i}.J*B{i}.H1*B{i}.H2*B{i}.q.v1)...
          + (B{i}.rho/2)*(B{i}.q.v2'*B{i}.J*B{i}.H1*B{i}.H2*B{i}.q.v2)...
          + (1/(2*B{i}.lam))*(B{i}.q.pr'*B{i}.J*B{i}.H1*B{i}.H2*B{i}.q.pr);

        berror(s,i) =...
            (B{i}.rho/2)*(B{i}.e.v1'*B{i}.J*B{i}.H1*B{i}.H2*B{i}.e.v1)...
          + (B{i}.rho/2)*(B{i}.e.v2'*B{i}.J*B{i}.H1*B{i}.H2*B{i}.e.v2)...
          + (1/(2*B{i}.lam))*(B{i}.e.pr'*B{i}.J*B{i}.H1*B{i}.H2*B{i}.e.pr);

      else
        c = B{i}.cubature;

        cv1 = c.V*B{i}.q.v1;
        cv2 = c.V*B{i}.q.v2;
        cpr = c.V*B{i}.q.pr;
        energy(s) = energy(s)...
         + (B{i}.rho/2)*(cv1(:)'*(c.W(:).*cv1(:)))...
         + (B{i}.rho/2)*(cv2(:)'*(c.W(:).*cv2(:)))...
         + (1/(2*B{i}.lam))*(cpr(:)'*(c.W(:).*cpr(:)));

        cev1 = c.V*B{i}.e.v1;
        cev2 = c.V*B{i}.e.v2;
        cepr = c.V*B{i}.e.pr;
        berror(s,i) = ...
         + (B{i}.rho/2)*(cev1(:)'*(c.W(:).*cev1(:)))...
         + (B{i}.rho/2)*(cev2(:)'*(c.W(:).*cev2(:)))...
         + (1/(2*B{i}.lam))*(cepr(:)'*(c.W(:).*cepr(:)));
      end

      error(s) = error(s) + berror(s,i);
    end
    error(s) = sqrt(error(s));
    berror(s,:) = sqrt(berror(s,:));

    if(nplot)
      subplot(2,1,1)
      semilogy(dt*(1:s),energy,'*k-')
      xlim([0 nsteps*dt])
      xlabel('time')
      ylabel('energy')

      subplot(2,1,2)
      semilogy(dt*(1:s),error,'*k-')
      xlim([0 nsteps*dt])
      xlabel('time')
      ylabel('error')

      drawnow
    end
  end
  time = dt*(1:nsteps);
  disp(min(energy(1:end-1)-energy(2:end)))
end


function [B,G] = domain_update_RK(B,G,t0,dt)

  RKA = [0, -567301805773/1357537059087, -2404267990393/2016746695238,...
         -3550918686646/2091501179385, -1275806237668/842570457699];
  RKB = [1432997174477/9575080441755, 5161836677717/13612068292357,...
         1720146321549/2090206949498, 3134564353537/4481467310338,...
         2277821191437/14882151754819];
  RKC = [0, 1432997174477/9575080441755, 2526269341429/6820363962896,...
         2006345519317/3224310063776, 2802321613138/2924317926251, 1];

  for rk_i = 1:length(RKA)
    t = t0+RKC(rk_i)*dt;
    [B,G] = domain_rhs(B,G,t);

    % do the update
    for i = 1:length(B)
      B{i}.dq.v1 = RKA(rk_i)*B{i}.dq.v1 + B{i}.rhs.v1;
      B{i}.dq.v2 = RKA(rk_i)*B{i}.dq.v2 + B{i}.rhs.v2;
      B{i}.dq.pr = RKA(rk_i)*B{i}.dq.pr + B{i}.rhs.pr;
      B{i}.q.v1  = B{i}.q.v1 + dt*RKB(rk_i)*B{i}.dq.v1;
      B{i}.q.v2  = B{i}.q.v2 + dt*RKB(rk_i)*B{i}.dq.v2;
      B{i}.q.pr  = B{i}.q.pr + dt*RKB(rk_i)*B{i}.dq.pr;
    end
  end
end

function [B,G] = domain_update_sync_AB(B,G,t0,dt,s)

  numBlocks = length(B);
  numGlue   = length(G);

  switch(s)
    case 1
      AB_A = [  1;0;0];
    case 2
      AB_A = [3/2;-1/2;0];
    otherwise
      AB_A = [23/12;-16/12;5/12];
  end

  % STEP 1: update blocks
  for i = 1:numBlocks
    Bup(i) = false;
    % update block
    B{i}.q.v1  = B{i}.q.v1 + dt*...
     (AB_A(1)*B{i}.rhs{1}.v1 + AB_A(2)*B{i}.rhs{2}.v1 + AB_A(3)*B{i}.rhs{3}.v1);
    B{i}.q.v2  = B{i}.q.v2 + dt*...
     (AB_A(1)*B{i}.rhs{1}.v2 + AB_A(2)*B{i}.rhs{2}.v2 + AB_A(3)*B{i}.rhs{3}.v2);
    B{i}.q.pr  = B{i}.q.pr + dt*...
     (AB_A(1)*B{i}.rhs{1}.pr + AB_A(2)*B{i}.rhs{2}.pr + AB_A(3)*B{i}.rhs{3}.pr);

    % update block time
    % B{i}.time = B{i}.time + dt;
    % set the fluxes
    if B{i}.sbp
      B{i} = hat_bc_sbp_block(B{i});
    end
  end

  % STEP 2: update fluxes
  for i = 1:numGlue

    for j = 1:G{i}.nB
      AB_G = AB_A;
      G{i}.q{j}.v  = G{i}.q{j}.v + dt*...
        (AB_G(1)*G{i}.rhs{j,1}.v + AB_G(2)*G{i}.rhs{j,2}.v...
         + AB_G(3)*G{i}.rhs{j,3}.v);
      G{i}.q{j}.p  = G{i}.q{j}.p + dt*...
        (AB_G(1)*G{i}.rhs{j,1}.p + AB_G(2)*G{i}.rhs{j,2}.p...
         + AB_G(3)*G{i}.rhs{j,3}.p);
    end

    % compute flux on the glue
    G{i} = compute_hat(G{i},B{G{i}.block(1)}.Z,B{G{i}.block(1)}.Z);

    % update the states
    for j = 1:G{i}.nB
      if B{G{i}.block(j)}.sbp
        B{G{i}.block(j)} = ...
          glue_to_sbp_block_flux(B{G{i}.block(j)},G{i},j);
      else
        B{G{i}.block(j)} = ...
          glue_to_dg_block_flux(B{G{i}.block(j)},G{i},j);
      end
    end
  end

  % STEP 3: update RHS blocks
  for i = 1:numBlocks
    % 2 to 3
    B{i}.rhs{3} = B{i}.rhs{2};

    % 1 to 2
    B{i}.rhs{2} = B{i}.rhs{1};

    % new RHS
    if B{i}.sbp
      B{i}.rhs{1}  = ...
        rhs_sbp_block(B{i},B{i}.q,B{i}.f.w,B{i}.f.e,B{i}.f.s,B{i}.f.n);
    else
      B{i}.rhs{1}  = rhs_dg_block(B{i});
    end
  end

  % STEP 4: update RHS glue
  for i = 1:numGlue
    for j = 1:G{i}.nB
      % 2 to 3
      G{i}.rhs{j,3} = G{i}.rhs{j,2};

      % 1 to 2
      G{i}.rhs{j,2} = G{i}.rhs{j,1};

      % pull RHS from block j
      if B{G{i}.block(j)}.sbp
        G{i}.rhs{j,1} = sbp_block_glue(B{G{i}.block(j)}.rhs{1},...
                                       B{G{i}.block(j)},G{i},j);
      else
        G{i}.rhs{j,1} = dg_block_glue(B{G{i}.block(j)}.rhs{1},...
                                      B{G{i}.block(j)},G{i},j);
      end
    end
  end
end
