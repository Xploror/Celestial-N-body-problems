clc;clear

dt = 0.01;
transfer = 0;

G = 2*10^(-6);
M = 19000000; % Jupiter  1.9*10^27
m_s1 = 893;  % Io  8.93*10^22
m_s2 = 480;  % Europa   4.8*10^22
m_s3 = 1480;  % Ganymede  14.8*10^22
m_s4 = 1075.9;  % Callisto 10.759*10^22
m_sat = 0.0001;

resonance_ratio = [2,1];
e_s1 = 0.0041;  % Io
e_s2 = 0.009;   % Europa
e_s3 = 0.0013;  % Ganymede
e_s4 = 0.0074;  % Callisto
e_sat = e_s2;
a_s1 = 1;  % Io 421700
a_s2 = 670900/421700;   % Europa  670900
a_s3 = 1070400/421700;   % Ganymede 1070400
a_s4 = 1882700/421700;   % Callisto 1882700
T_s1 = ((4*pi^2/(G*M))*a_s1^3)^(1/2);
T_s2 = ((4*pi^2/(G*M))*a_s2^3)^(1/2);
T_s3 = ((4*pi^2/(G*M))*a_s3^3)^(1/2);
T_s4 = ((4*pi^2/(G*M))*a_s4^3)^(1/2);

T_sat = resonance_ratio(2)*T_s2/resonance_ratio(1);
a_sat_target = a_s2*(T_sat/T_s2)^(2/3);

r_s1_init = [a_s1*(1+e_s1),0];
r_s2_init = [a_s2*(1+e_s2),0];
r_s3_init = [a_s3*(1+e_s3),0];
r_s4_init = [a_s4*(1+e_s4),0];
v1_init = [0,sqrt(G*M*((2/r_s1_init(1))-(1/a_s1)))];
v2_init = [0,sqrt(G*M*((2/r_s2_init(1))-(1/a_s2)))];
v3_init = [0,sqrt(G*M*((2/r_s3_init(1))-(1/a_s3)))];
v4_init = [0,sqrt(G*M*((2/r_s4_init(1))-(1/a_s4)))];

%% Resonance Plotting
% %%%%% Transfer mechanism (Starts with an arbitrary parking orbit around primary 1)
% e1_sat = 0; % Circular orbit (Hyperparameter)
% a1_sat = a_sat_target*(1-e_sat); % (Hyperparameter)
% T1_sat = ((4*pi^2/(G*M))*a1_sat^3)^(1/2);
% r_sat_init = [a1_sat*(1+e1_sat),0];
% v_sat_init = sqrt(G*M*((2/norm(r_sat_init))-(1/a1_sat)))*cross([0,0,1], [r_sat_init,0])*[1 0 0; 0 1 0]'/norm(r_sat_init);
% t_trans = 40*T1_sat/2;  % (Hyperparameter)
% 
% %Dynamic propogation
% r_s1 = r_s1_init;
% r_s2 = r_s2_init;
% r_s3 = r_s3_init;
% v_s1 = v1_init;
% v_s2 = v2_init;
% v_s3 = v3_init;
% r_sat = r_sat_init;
% v_sat = v_sat_init;
% t = 0;
% 
% for i=0:dt:100
%     acc_s1 = -((G*M/norm(r_s1)^3)*r_s1 + (G*m_sat/norm(r_sat-r_s1)^3)*(r_s1 - r_sat) + (G*m_s2/norm(r_s2-r_s1)^3)*(r_s1 - r_s2) + (G*m_s3/norm(r_s3-r_s1)^3)*(r_s1 - r_s3));
%     acc_s2 = -((G*M/norm(r_s2)^3)*r_s2 + (G*m_sat/norm(r_sat-r_s2)^3)*(r_s2 - r_sat) + (G*m_s1/norm(r_s2-r_s1)^3)*(r_s2 - r_s1) + (G*m_s3/norm(r_s3-r_s2)^3)*(r_s2 - r_s3));
%     acc_s3 = -((G*M/norm(r_s3)^3)*r_s3 + (G*m_sat/norm(r_sat-r_s3)^3)*(r_s3 - r_sat) + (G*m_s2/norm(r_s2-r_s3)^3)*(r_s3 - r_s2) + (G*m_s1/norm(r_s3-r_s1)^3)*(r_s3 - r_s1));
%     acc_sat = -((G*M/norm(r_sat)^3)*r_sat + (G*m_s1/norm(r_sat-r_s1)^3)*(r_sat - r_s1) + (G*m_s2/norm(r_sat-r_s2)^3)*(r_sat - r_s2) + (G*m_s3/norm(r_sat-r_s3)^3)*(r_sat - r_s3));
%     v_s1 = v_s1 + acc_s1*dt;
%     v_s2 = v_s2 + acc_s2*dt;
%     v_s3 = v_s3 + acc_s3*dt;
%     v_sat = v_sat + acc_sat*dt;
%     r_s1 = r_s1 + v_s1*dt;
%     r_s2 = r_s2 + v_s2*dt;
%     r_s3 = r_s3 + v_s3*dt;
%     r_sat = r_sat + v_sat*dt;
%     
%     %%%%%%%%%%%%% Transfer mechanism %%%%%%%%%%%%%%%
% %     if t >= t_trans && transfer == 0
% %         a_sat = (norm(r_sat) + a_sat_target*(1+e_sat))/2;
% %         v_sat = sqrt(G*M*((2/norm(r_sat))-(1/a_sat_target)))*cross([0,0,1], [r_sat,0])*[1 0 0; 0 1 0]'/norm(r_sat_init);  % This line of code is analogous to delta_v thrust point
% %         transfer = 1;
% %     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     figure(1)
%     plot(r_s1(1), r_s1(2), '.r', r_s2(1), r_s2(2), '.b', r_s3(1), r_s3(2), '.g', 0,0, '.k', r_sat(1), r_sat(2), '.y', 'MarkerSize', 20);
%     xlim([-3,3])
%     ylim([-3,3])
%     
%     %%%%%% This is the plot to visualize the celestial crossing of the
%     %%%%%% satellite 
% %     if abs((r_s1/norm(r_s1))*(r_sat/norm(r_sat))' - 1) <= 0.0005
% %         figure(2)
% %         plot(r_s1(1), r_s1(2), '.b', r_sat(1), r_sat(2), '.g')
% %         xlim([-5,5])
% %         ylim([-5,5])
% %         hold on;
% %     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % Contour plot
%     x = -3:0.05:3;
%     y = -3:0.05:3;
%     [X, Y] = meshgrid(x,y);
%     len_x = length(X);
%     len_y = length(Y);
%     U = -(G*M*10^(-4)./sqrt(X.^2 + Y.^2) + G*m_s1./sqrt((X-r_s1(1)*ones(len_x,len_x)).^2 + (Y-r_s1(2)*ones(len_y,len_y)).^2) + G*m_s2./sqrt((X-r_s2(1)*ones(len_x,len_x)).^2 + (Y-r_s2(2)*ones(len_y,len_y)).^2) + G*m_s3./sqrt((X-r_s3(1)*ones(len_x,len_x)).^2 + (Y-r_s3(2)*ones(len_y,len_y)).^2));
%     figure(2)
%     contourf(X, Y, U);
%     hold on;
%  
%     R = sqrt(X.^2+Y.^2);
%     acc_obj_x = -((G*10^(-3)*M./R^3).*X + (G*m_s1./((X-r_s1(1)*ones(len_x,len_x)).^2 + (Y-r_s1(2)*ones(len_y,len_y)).^2)^(3/2)).*(X - r_s1(1)*ones(len_x,len_x)) + (G*m_s2./((X-r_s2(1)*ones(len_x,len_x)).^2 + (Y-r_s2(2)*ones(len_y,len_y)).^2)^(3/2)).*(X - r_s2(1)*ones(len_x,len_x)) + (G*m_s3./((X-r_s3(1)*ones(len_x,len_x)).^2 + (Y-r_s3(2)*ones(len_y,len_y)).^2)^(3/2)).*(X - r_s3(1)*ones(len_x,len_x)));
%     acc_obj_y = -((G*10^(-3)*M./R^3).*Y + (G*m_s1./((X-r_s1(1)*ones(len_x,len_x)).^2 + (Y-r_s1(2)*ones(len_y,len_y)).^2)^(3/2)).*(Y - r_s1(1)*ones(len_x,len_x)) + (G*m_s2./((X-r_s2(1)*ones(len_x,len_x)).^2 + (Y-r_s2(2)*ones(len_y,len_y)).^2)^(3/2)).*(Y - r_s2(1)*ones(len_x,len_x)) + (G*m_s3./((X-r_s3(1)*ones(len_x,len_x)).^2 + (Y-r_s3(2)*ones(len_y,len_y)).^2)^(3/2)).*(Y - r_s3(1)*ones(len_x,len_x)));
%     quiver(X, Y, acc_obj_x, acc_obj_y);
%     hold off;
%     pause(0.0000001)
%     
%     t = t + dt;
% end

%% Jupiter-Io/Europa/Ganymede/Callisto 3BP in x-y plane analysis

% mu1 = m_s3/M;
% tStep = 1*10^(-3);
% 
% % L1 (between Jupiter and Io)
% f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 + mu*(x+mu)^2;
% g = @(x) f(x, mu1);
% l1_pos_io = fzero(g, 0);
% 
% % Eigenspace construction for L1 
% omega = (1-mu1)/abs(l1_pos_io+mu1)^3 + mu1/abs(l1_pos_io-(1-mu1))^3;
% a = 1+2*omega;
% b = 1-omega;
% A = [0 0 1 0
%      0 0 0 1
%      a 0 0 2
%      0 b -2 0];   % Characteristic matrix at L1 position (2D xy plane)
%  
% [l1_eigvec_io, l1_eigv_io] = eig(A);
% l1_eigvec_io = real(l1_eigvec_io); % Cause these eigenvectors will be used in real space
% 
% T_p = 3*pi/abs(l1_eigv_io(3,3)); % Third eigenvalue of the system in the form +-i*nu
% del_x = 8*10^(-3);
% tspan = 0:0.001:T_p;
% x0 = [l1_pos_io;0;0;0] + del_x*l1_eigvec_io(:,3);
% 
% [a,b, STM_b] = diffcorr_l1(x0, mu1, tStep);
% [t_dc, X_dc] = ode45(@(t,x) PRTBP_dyn(t,x,mu1), 0:0.001:2*b, a);
% 
% % load('L1_J-Io_8e-2.mat');
% % load('L1_STM_J-Io_8e-2.mat');
% 
% plot(X_dc(:,1), X_dc(:,2), '-k'); hold on; pause(0.0004)
% %%%%%%% Propagating over the stable/unstable manifold of STM matrix on corrected initial conditions from DC
% % Invariant manifold calculation for specific orbits (Inputs --> num_points, X_dc, STM_b, disp)
%  disp = 10^(-4);
%  [vec, val] = eig(STM_b);
%  num_points = 256;
%  ps_points = [];
%  tSpan_mani_backw = 0:-0.001:-(length(X_dc)*tStep);
%  tSpan_mani_forw = 0:0.001:(length(X_dc)*tStep);
%  for i=1:num_points
%      x0_smani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,2)/norm(vec(:,2));
%      x0_unsmani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,2)/norm(vec(:,2));
%      [t_sm_b, X_sm_b] = ode45(@(t,x) PRTBP_dyn(t,x,mu1), tSpan_mani_backw, x0_smani);
%      [t_unsm_b, X_unsm_b] = ode45(@(t,x) PRTBP_dyn(t,x,mu1), tSpan_mani_forw, x0_unsmani);
%      % Poincare preprocessing
%      [ps_v, ps_index] = min(abs(X_sm_b(:,1)-(1-mu1)));
%      ps_points = [ps_points; [X_sm_b(ps_index,2) X_sm_b(ps_index,4)]];
%      %%%%%%%%%%%%%%%%%%%%%%%%
%      plot(X_sm_b(:,1), X_sm_b(:,2), '-r', 1-mu1, 0, 'ok'); hold on;   %Only stable manifold plot
%      plot(X_unsm_b(:,1), X_unsm_b(:,2), '-b', 1-mu1, 0, 'ok')   %Only unstable manifold plot
%  end

%% Jupiter-Io/Europa/Ganymede/Callisto 3BP in x-y-z region analysis

mu1 = m_s1/M;
tStep = 1*10^(-3);

% L1 (between Jupiter and moons)
f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 + mu*(x+mu)^2;
g = @(x) f(x, mu1);
l1_pos_io = fzero(g, 0);

% L2 (beyond Jupiter and moons)
f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 - mu*(x+mu)^2;
g = @(x) f(x, mu1);
l2_pos_io = fzero(g, 0);

%% Eigenspace construction for L1 in 3D spatial coordinates (6 dimensional state-space)
A = crtbp_J(l1_pos_io, 0, 0, mu1);
 
[l1_eigvec_io, l1_eigv_io] = eig(A);
l1_eigvec_io = real(l1_eigvec_io); % Cause these eigenvectors will be used in real space

T_p = 3*pi/abs(l1_eigv_io(3,3)); % Third eigenvalue of the system in the form +-i*nu
del_x = -5.48*10^(-2);
tspan = 0:0.001:T_p;
x0 = [l1_pos_io;0;0;0;0;0] + del_x*l1_eigvec_io(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eigenspace construction for L2 in 3D spatial coordinates (6 dimensional state-space)
% A = crtbp_J(l2_pos_io, 0, 0, mu1);
%  
% [l2_eigvec_io, l2_eigv_io] = eig(A);
% l2_eigvec_io = real(l2_eigvec_io); % Cause these eigenvectors will be used in real space
% 
% T_p = 3*pi/abs(l2_eigv_io(3,3)); % Third eigenvalue of the system in the form +-i*nu
% del_x = 1*10^(-3);
% tspan = 0:0.001:T_p;
% x0 = [l2_pos_io;0;0;0;0;0] + del_x*l2_eigvec_io(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%% Comment this subseciton if differential correction is to be ignored once the stored X_dc and STM matrrix is saved for perticular configuration
zz = [-0.001];
for i=1:1
    x0 = [l1_pos_io;0;zz(i);0;0;0] + del_x*l1_eigvec_io(:,3);
[a,b, STM_b] = halo_DC(x0, mu1, tStep);
[t_dc, X_dc] = ode45(@(t,x) statediffeq(t,x,mu1), 0:0.001:2*b, a);

plot3(X_dc(:,1), X_dc(:,2), X_dc(:,3), '-k'); hold on;
pause(0.003)
hold on;
end

%% Periodic Orbits and Invariant manifolds

% load('L2_J-G_1e-1_3d.mat');
% load('L1_J-G_8e-3_3d.mat');


%%%%%%% Propagating over the stable/unstable manifold of STM matrix on corrected initial conditions from DC
% Invariant manifold calculation for specific orbits (Inputs --> num_points, X_dc, STM_b, disp)
%  disp = 1*10^(-3);
%  %[vec, val] = eig(STM_b);
%  num_points = 64;
%  ps_points = [];
%  tSpan_mani_backw = 0:-0.001:-2*(length(X_dc)*tStep);
%  tSpan_mani_forw = 0:0.001:2*(length(X_dc)*tStep);
%  l2_orb = load('L2_J-G_1e-1_3d.mat');
%  l1_orb = load('L1_J-G_1e-1_3d.mat');
%  l2_STM = load('L2_STM_J-G_1e-1_3d.mat');
%  l1_STM = load('L1_STM_J-G_1e-1_3d.mat');
%  for k=1:2
%      if k==1
%          X_dc = l1_orb.X_dc;
%          [vec, val] = eig(l1_STM.STM_b);
%      else
%          X_dc = l2_orb.X_dc;
%          [vec, val] = eig(l2_STM.STM_b);
%      end
%  for i=1:num_points
%      x0_smani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,1)/norm(vec(:,3));
%      x0_unsmani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,1)/norm(vec(:,3));
%      [t_sm_b, X_sm_b] = ode45(@(t,x) statediffeq(t,x,mu1), tSpan_mani_backw, x0_smani);
%      [t_unsm_b, X_unsm_b] = ode45(@(t,x) statediffeq(t,x,mu1), tSpan_mani_forw, x0_unsmani);
%      % Poincare preprocessing
% %      [ps_v, ps_index] = min(abs(X_unsm_b(:,1)-(1-mu1)));
% %      ps_points = [ps_points; [X_unsm_b(ps_index,2) X_unsm_b(ps_index,4)]];
%      %%%%%%%%%%%%%%%%%%%%%%%%
%      plot3(X_sm_b(:,1), X_sm_b(:,2), X_sm_b(:,3), '-r', 1-mu1, 0, 0, 'ok'); hold on;   %Only stable manifold plot
%      xlabel('x');ylabel('y');zlabel('z')
%      plot3(X_unsm_b(:,1), X_unsm_b(:,2), X_unsm_b(:,3), '-b', 1-mu1, 0, 0, 'ok')   %Only unstable manifold plot
%      xlim([0.8,1.2])
%  end
%  pause(0.03)
%  end
