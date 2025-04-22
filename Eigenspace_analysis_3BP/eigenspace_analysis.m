clc;clear

% Solving lienarized dimensionless equation of motion of CR3BP to find the
% collinear points L1, L2 and L3  (y=0 & z=0)

mu = 0.0121; % Cislunar regime
tStep = 1*10^(-3);  % 5*10^-4 works for very small perturbations =< 10^(-2)

% L1 (between m1 and m2)
f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 + mu*(x+mu)^2;
g = @(x) f(x, mu);
l1_pos = fzero(g, 0); % x-position in rotating frame of L1

% Find eigenvalues and eigenvector of dynamics around L1 to study manifolds
omega = (1-mu)/abs(l1_pos+mu)^3 + mu/abs(l1_pos-(1-mu))^3;
a = 1+2*omega;
b = 1-omega;
A = [0 0 1 0
     0 0 0 1
     a 0 0 2
     0 b -2 0];   % Characteristic matrix at L1 position (2D xy plane)
 
 [l1_eigvec, l1_eigv] = eig(A);
 l1_eigvec = real(l1_eigvec); % Cause these eigenvectors will be used in real space
 
 
 %% Planar Lyapunov Orbit (Periodic orbit with imaginary eigenvalues)
 
 T_p = 3*pi/abs(l1_eigv(3,3)); % Third eigenvalue of the system in the form +-i*nu
 %del_x_fam = [3*10^(-3) 8*10^(-2) 5*10^(-1) 6*10^(-1) 7*10^(-1)];
 del_x_fam = [8*10^(-2)];
 tspan = 0:0.001:T_p;
 
 STM_fam = [];
 X_dc_fam = [];
 path_len = [];
 x0 = [l1_pos;0;0;0] + del_x_fam(1)*l1_eigvec(:,3);
 % differential correction for family of periodic orbits
 for i=1:length(del_x_fam)
    x0 = [l1_pos;0;0;0] + del_x_fam(i)*l1_eigvec(:,3); 
    [a,b, STM_b] = diffcorr_l1(x0, mu, tStep, del_x_fam(i));
    STM_fam = [STM_fam STM_b];
    [t_dc, X_dc] = ode45(@(t,x) PRTBP_dyn(t,x,mu), 0:0.001:2*b, a);
    plot(X_dc(:,1), X_dc(:,2), '-k', l1_pos, 0, 'ob') % Plots periodic orbit
    hold on;   
    pause(1)
    X_dc_fam = [X_dc_fam; X_dc];
    path_len = [path_len, length(t_dc)];
 end

%  [a,b, STM_b] = diffcorr_l1(x0, mu, tStep);
%  [t_dc, X_dc] = ode45(@(t,x) PRTBP_dyn(t,x,mu), 0:0.001:2*b, a);

%  load('L1_periodic')
%  load('L1_STM')
 
 plot(X_dc(:,1), X_dc(:,2), '-k', l1_pos, 0, 'ob') % Plots periodic orbit
 hold on;
 
 %%%%%%% Prpagating over the stable/unstable manifold of STM matrix on corrected initial conditions from DC
 disp = -10^(-3);
 [vec, val] = eig(STM_b);
 x0_smani = X_dc(1,:)' + disp*vec(:,2)/norm(vec(:,2)); % Taking eigenvector direction for second eigenvalue since that's the stable direction 
 x0_unsmani = X_dc(1,:)' + disp*vec(:,1)/norm(vec(:,1)); % Taking eigenvector direction for first eigenvalue since that's the unstable direction  
 %%%%% Forward integration in time %%%%
%  tSpan_mani_forw = 0:0.001:(length(X_dc)*tStep);
%  [t_sm_f, X_sm_f] = ode45(@(t,x) PRTBP_dyn(t,x,mu), tSpan_mani_forw, x0_smani);
%  [t_unsm_f, X_unsm_f] = ode45(@(t,x) PRTBP_dyn(t,x,mu), tSpan_mani_forw, x0_unsmani);
 %%%%% Backward integration in time %%%%
 tSpan_mani_backw = 0:-0.001:-(length(X_dc)*tStep);
 tSpan_mani_forw = 0:0.001:(length(X_dc)*tStep);
 % Invariant manifold calculation for specific orbits 
 num_points = 64;
 ps_points = [];
 for i=1:num_points
     x0_smani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,2)/norm(vec(:,2));
     x0_unsmani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,1)/norm(vec(:,1));
     [t_sm_b, X_sm_b] = ode45(@(t,x) PRTBP_dyn(t,x,mu), tSpan_mani_backw, x0_smani);
     [t_unsm_b, X_unsm_b] = ode45(@(t,x) PRTBP_dyn(t,x,mu), tSpan_mani_forw, x0_unsmani);
     % Poincare preprocessing
     [ps_v, ps_index] = min(abs(X_unsm_b(:,1)-(1-mu)));
     ps_points = [ps_points; [X_unsm_b(ps_index,2) X_unsm_b(ps_index,4)]];
     %%%%%%%%%%%%%%%%%%%%%%%%
     plot(X_sm_b(:,1), X_sm_b(:,2), '-r', 1-mu, 0, 'ok'); hold on;   %Only stable manifold plot
     plot(X_unsm_b(:,1), X_unsm_b(:,2), '-b', 1-mu, 0, 'ok')
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %  plot(X_sm_f(:,1), X_sm_f(:,2), '-g', 1-mu, 0, 'ob', X_unsm_f(:,1), X_unsm_f(:,2), '-r', X_sm_b(:,1), X_sm_b(:,2), '--g', X_unsm_b(:,1), X_unsm_b(:,2), '--r')
% % % % %  plot(X_sm_f(:,1), X_sm_f(:,2), '-g', X_sm_b(:,1), X_sm_b(:,2), '--g', 1-mu, 0, 'ob'); hold on; pause(4)  %Only stable manifold plot
% % % % %  plot(X_unsm_f(:,1), X_unsm_f(:,2), '-r', X_unsm_b(:,1), X_unsm_b(:,2), '--r', 1-mu, 0, 'ob')  %Only unstable manifold plot

%% Chaos Analysis
%%%%%%% Poincare Section %%%%%%%%%%%%%%%%%%%
figure(2)
plot(ps_points(:,1), ps_points(:,2), 'ob');

 
