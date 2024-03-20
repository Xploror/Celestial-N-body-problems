% Make a high fidelity model example 4 and 5 body problem
% Idea: Look into Saturn low energy transfers into Titan, include ring perturbations too

clc;clear;

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

n_1 = sqrt(G*(M+m_s1)/(a_s1^3));
n_2 = sqrt(G*(M+m_s2)/(a_s2^3));
n_3 = sqrt(G*(M+m_s3)/(a_s3^3));
%n_3 = 0.5514;
n_4 = sqrt(G*(M+m_s4)/(a_s4^3));
Rot_mat = @(n,t) [cos(n*t) -sin(n*t); sin(n*t) cos(n*t)];

r_s1_init = [a_s1*(1+e_s1),0];
r_s2_init = [a_s2*(1+e_s2),0];
r_s3_init = [a_s3*(1+e_s3),0];
r_s4_init = [a_s4*(1+e_s4),0];
v1_init = [0,sqrt(G*M*((2/r_s1_init(1))-(1/a_s1)))];
v2_init = [0,sqrt(G*M*((2/r_s2_init(1))-(1/a_s2)))];
v3_init = [0,sqrt(G*M*((2/r_s3_init(1))-(1/a_s3)))];
v4_init = [0,sqrt(G*M*((2/r_s4_init(1))-(1/a_s4)))];

Hamilt = @(x,y,z,vx,vy,vz,mu,E) 0.5*vx^2 + 0.5*vy^2 + 0.5*vz^2 - 0.5*x^2 - 0.5*y^2 - 0.5*z^2 - (1-mu)/sqrt((x+mu)^2 + y^2 + z^2) - (mu)/sqrt((x+mu-1)^2 + y^2 + z^2) - E;

%% Poincare Section
% Poincare Section for orbits around the primary within the region between primary and secondary

% num_iter = 1000;
% mu = m_s1/(M+m_s1);
% 
% tspan = 0:0.01:20;
% opt = odeset('Events',@ps_crossing);
% 
% for i=1:num_iter
%     x_init = 0.4 + rand([1,1])*(0.5-0.4); % [0.25, 0.6]
%     E = -3 + rand([1,1])*(-2.5+3);  % [-2,-3]
%     find_vy = @(vy) Hamilt(x_init,0,0.1,0.01,vy,0,mu,E);
%     syms vy
%     vy_init = solve(find_vy, vy);
%     x_sat_init = [x_init,0,0.1,0.01,double(abs(vy_init(1))),0];
%     [t, X, te, Xe, ie] = ode45(@(t,x) statediffeq(t,x,mu), tspan, x_sat_init, opt);
%     
%     plot(Xe(:,3),Xe(:,6),'.b')
%     %plot3(X(:,1),X(:,2),X(:,3),'.b')
%     %xlim([0.3,0.4])
%     pause(0.0003)
%     hold on;
% end

%% Stroboscopic map for PCC4BP

num_points = 100;

tStep = 5*10^-3;
lce_tspan = 200 ;

mu1 = m_s2/(M+m_s2);
mu2 = m_s3/(M+m_s2);
p = [(a_s3/a_s2)+mu1,0];

func = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 - mu*(x+mu)^2;
g = @(x) func(x, mu1);
l2_pos = fzero(g, 0);
del_x = -1*10^(-3);
[t_dc, X_dc, STM] = libration_orbit_calculator(l2_pos, del_x, mu1, 0); % gives the periodic orbit values
[vec, val] = eig(STM);
%--------------- Calculation of periodic orbit done ------------------%

tspan = 0:0.01*tStep:lce_tspan;
opt = odeset('Events',@theta_crossing);
disp = 10^(-4);

f = @(x) PCC4BP_eqn(1, x, mu1, mu2, p, n_3, 3);
figure(2)
for i=1:num_points
    x_init = X_dc(int32((i-1)*length(X_dc)/num_points+1),:)'; %+ disp*vec(:,1)/norm(vec(:,1));
    [t, X, te, Xe, ie] = ode45(@(t,x) PCC4BP_eqn(t, x, mu1, mu2, p, n_3, 3), tspan, x_init, opt);
    Xe = X(find(rem(X(:,5),2*pi)<2.5*10^(-4)==1),:);
    
    plot(Xe(:,1),Xe(:,2),'.b')
    pause(0.0000003)
    hold on;
end

% X_ = X_dc(1,:)';
% for i=1:length(X_dc)
%    [t_, X_] = ode45(@(t,x) PCC4BP_eqn(t, x, mu1, 0, p, n_3, 3), [0,tStep], X_);
%    X_ = X_(end,:)';
%    X_3 = [-mu1+p(1)*cos(X_(5)) p(1)*sin(X_(5))];
%    
%    figure(1)
%    plot(X_(1), X_(2), '.r')
%    hold on;
%    
%    figure(2)
%    plot(X_3(1), X_3(2), '.b')
%    hold on;
%     
% end


%% FLI map
% Taking L2 periodic orbit manifolds analysis using Planar 4BP model --- Finding each LCE for initial condition x0(1:2) and tangent vector v0 (stable and unstable direction) 

% tStep = 5*10^-3;
% lce_tspan = 50;
% 
% mu1 = m_s2/(M+m_s2);
% mu2 = m_s3/(M+m_s2);
% p = [(a_s3/a_s2)+mu1,0];
% 
% func = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 - mu*(x+mu)^2;
% g = @(x) func(x, mu1);
% l2_pos = fzero(g, 0);
% del_x = -1*10^(-3);
% [t_dc, X_dc, STM] = libration_orbit_calculator(l2_pos, del_x, mu1, 0); % gives the periodic orbit values
% [vec, val] = eig(STM);
% 
% x0 = X_dc(end,:)' + 4*10^(-3)*vec(:,1); % Initial state on the STM stable eigenvector i.e on the stable manifold 
% % v_s0 = [-vec(2,2);vec(1,2);-vec(4,2);vec(3,2);vec(5,2)];
% % v_uns0 = [-vec(2,1);vec(1,1);-vec(4,1);vec(3,1);vec(5,2)];
% v_s0 = [-vec(2,2);vec(1,2);vec(3:4,2);vec(5,2)];
% v_uns0 = [-vec(2,1);vec(1,1);vec(3:4,1);vec(5,1)];
% x = x0;
% v_s = v_s0;
% v_uns = v_uns0;
% f = @(x) PCC4BP_eqn(1, x, mu1, mu2, p, n_3, 3);
% g = @(J,v) J*v;
% LCE_s = []; % Stable tangent vector LCE 
% LCE_uns = []; % Unstable tangent vector LCE
% 
% for t = tStep:tStep:lce_tspan
%    
%    J = PCC4BP_J(x(1), x(2), x(5), mu1, mu2, p);
%    
% %    v_s = J*v_s;
% %    v_uns = J*v_uns;
%    % RK 6th Order Numerical Integration  
%    k1 = tStep*g(J,v_s);
%    k2 = tStep*g(J,v_s+k1);
%    k3 = tStep*g(J,v_s+(3*k1+k2)/8);
%    k4 = tStep*g(J,v_s+(8*k1+2*k2+8*k3)/27);
%    k5 = tStep*g(J,v_s+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
%    k6 = tStep*g(J,v_s+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
%    k7 = tStep*g(J,v_s+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
%    v_s = v_s + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
%    
%    k1 = tStep*g(J,v_uns);
%    k2 = tStep*g(J,v_uns+k1);
%    k3 = tStep*g(J,v_uns+(3*k1+k2)/8);
%    k4 = tStep*g(J,v_uns+(8*k1+2*k2+8*k3)/27);
%    k5 = tStep*g(J,v_uns+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
%    k6 = tStep*g(J,v_uns+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
%    k7 = tStep*g(J,v_uns+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
%    v_uns = v_uns + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
%    
%    LCE_s = [LCE_s (log(norm(v_s)/norm(v_s0)))];
%    LCE_uns = [LCE_uns (log(norm(v_uns)/norm(v_uns0)))]; 
%    % normalize perturbation tangent vectors
%    v_s = v_s/norm(v_s);
%    v_uns = v_uns/norm(v_uns);
%     
%    k1 = tStep*f(x);
%    k2 = tStep*f(x+k1);
%    k3 = tStep*f(x+(3*k1+k2)/8);
%    k4 = tStep*f(x+(8*k1+2*k2+8*k3)/27);
%    k5 = tStep*f(x+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
%    k6 = tStep*f(x+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
%    k7 = tStep*f(x+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
%    x = x + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
% end

%% Twist map (Keplerian) of PCR3BP 




%% Twist map of 4 Body Problem (PCCR4BP)









%% Defined functions 

function [position,isterminal,direction] = ps_crossing(t,y)
  position = y(2); % The value that we want to be zero --- Poincare section existing in the y=0 line and x>0
  isterminal = 0;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction -- moving towards positive y direction
end

function [position,isterminal,direction] = theta_crossing(t,y)
 position = (rem(y(5),2*pi)>2.5*10^(-5)); % Want theta=0 crossing
 isterminal = 0;
 direction = 0;
end