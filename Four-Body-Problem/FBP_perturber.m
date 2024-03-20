clc;clear;

dt = 0.01;

M = 19000000;
m_1 = 893;
m_2 = 480;
m_3 = 1480;
m_4 = 1075.9;
e_1 = 0.004;
e_2 = 0.009;
e_3 = 0.0013;
e_4 = 0.0074;
a_1 = 1.97; 
a_2 = 3.135; 
a_3 = 5;
a_4 = 8.79;

tau = 3*10^(2); % Scaling of time --> tau sec/sec
dt = tau*0.01;
G = 2*10^(-6);
%G = 6.7*10^(-11)*(5/(1.07*10^9))^3*(10^20)*(tau)^2;   % Taking care of rescaling in length, mass and time

laplace_resonance = [4,2,1]; % Io-Europa-Ganymede
T_1 = ((4*pi^2/(G*M))*a_1^3)^(1/2);
T_2 = ((4*pi^2/(G*M))*a_2^3)^(1/2);
T_3 = ((4*pi^2/(G*M))*a_3^3)^(1/2);
T_4 = ((4*pi^2/(G*M))*a_4^3)^(1/2);

% Rotating to Inertial %%%%%
n_1 = sqrt(G*(M+m_1)/(a_1^3));
n_2 = sqrt(G*(M+m_2)/(a_2^3));
n_3 = sqrt(G*(M+m_3)/(a_3^3));
n_4 = sqrt(G*(M+m_4)/(a_4^3));
Rot_mat = @(n,t) [cos(n*t) -sin(n*t); sin(n*t) cos(n*t)];

mu1 = m_2/(M+m_2);
mu2 = m_3/(M+m_2);
%mu2 = 0;
% [t1, X1] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[0,a_3/a_2],n_3,3),0:0.001:60,[0.8,0.1,-0.01,0.6,0]);
% [t2, X2] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[0,a_3/a_2],n_3,3),0:0.001:60,[0.8,0.1,-0.01,0.6,0]);
% X_1 = zeros([length(X1(:,1)),2]);
% X_2 = zeros([length(X2(:,1)),2]);
% for i=1:length(t1)
%     X_1(i,:) = a_2*(Rot_mat(n_2,t1(i))*X1(i,1:2)')';
%     X_2(i,:) = a_2*(Rot_mat(n_2,t2(i))*X2(i,1:2)')';
% end
% plot(X_1(:,1),X_1(:,2),'k');hold on; plot(X_2(:,1),X_2(:,2),'r')

tStep = 1*10^(-3);

% L1 (between Jupiter and moons)
f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 + mu*(x+mu)^2;
g = @(x) f(x, mu1);
l1_pos_io = fzero(g, 0);

% L2 (beyond Jupiter and moons)
f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 - mu*(x+mu)^2;
g = @(x) f(x, mu1);
l2_pos_io = fzero(g, 0);

A_3d = crtbp_J(l2_pos_io, 0, 0, mu1);
A = [A_3d(1:2,1:2) A_3d(1:2,4:5); A_3d(4:5,1:2) A_3d(4:5,4:5)];
 
[l2_eigvec_io, l2_eigv_io] = eig(A);
l2_eigvec_io = real(l2_eigvec_io); % Cause these eigenvectors will be used in real space

T_p = 3*pi/abs(l2_eigv_io(3,3)); % Third eigenvalue of the system in the form +-i*nu
del_x = -1*10^(-3);
tspan = 0:0.001:T_p;
x0 = [l2_pos_io;0;0;0] + del_x*l2_eigvec_io(:,3);
x0 = [x0; pi/2]; % The fifth state is the angle of third primary with first 
[a,b, STM_b] = diffcorr_l2_FBP(x0, mu1, 0, [(a_3/a_2)+mu1,0], n_3, tStep);
[t_dc, X_dc] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_2)+mu1,0], n_3, 3), 0:0.001:2*b, a);
X_dc(:,5) = x0(5)*ones([length(X_dc(:,5)),1]);  % third primary angle for all the entries of periodic orbit should have same value as the initial condition x0 

figure(1)
hold on;
plot(X_dc(:,1), X_dc(:,2),'-k', 1-mu1, 0, 'ok'); hold on;
pause(0.000000001)
hold on;

 disp = 1*10^(-3);
 [vec, val] = eig(STM_b);
 num_points = 1;
 %ps_points = [];
 %tSpan_mani_backw = 0:-0.001:-1.5*(length(X_dc)*tStep);
 tSpan_mani_forw = 0:0.001:1*(length(X_dc)*tStep);
 for i=1:num_points
     %x0_smani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,2)/norm(vec(:,2));
     x0_unsmani = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)' + disp*vec(:,1)/norm(vec(:,1));
     %[t_sm_b, X_sm_b] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_2)+mu1,0], n_3, 3), tSpan_mani_backw, x0_smani);
     [t_unsm_b, X_unsm_b] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[(a_3/a_2)+mu1,0], n_3, 3), tSpan_mani_forw, x0_unsmani);
     % Poincare preprocessing
%      [ps_v, ps_index] = min(abs(X_unsm_b(:,1)-(1-mu1)));
%      ps_points = [ps_points; [X_unsm_b(ps_index,2) X_unsm_b(ps_index,4)]];
     %%%%%%%%%%%%%%%%%%%%%%%%
%      figure(1)
%      plot(X_sm_b(:,1), X_sm_b(:,2), '-r'); hold on;   %Only stable manifold plot
%      xlabel('x');ylabel('y');
%      plot(X_unsm_b(:,1), X_unsm_b(:,2), '-b')   %Only unstable manifold plot
%      pause(0.00001)
 end
 
 % Visualizing third primary in rotating frame (only for one manifold for one specific initial condition from periodic orbit)
 for i=1:length(X_unsm_b(:,1))
     X_3 = [mu1+((a_3/a_2)+mu1)*cos(X_unsm_b(i,5)) ((a_3/a_2)+mu1)*sin(X_unsm_b(i,5))];
     figure(1)
     plot(X_3(1), X_3(2), '.r', X_unsm_b(i,1), X_unsm_b(i,2), '.b', 1-mu1, 0, 'ok'); 
     xlabel('x');ylabel('y');
     pause(0.000000001)
 end
 
 % 