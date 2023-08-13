% CisLunar Numerical Simulation and Plotting
clc;clear;

%% Variables

epoch = 5000;
mu = 1.215*10^(-2);  % mass ratio
major1 = [-mu 0 0];
major2 = [(1-mu) 0 0];
t = 0;
dt = 2*pi*10^(-6)*5*60/2.361;  % taking 5 min time interval for Earth-Moon system
da = 2*pi*10^(-6)*5*60/2.361; 
%dt = 0.001;  % Time interval 0.1 was working fine
%da = 0.001; % angle change
NB = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

%% Initial Conditions (Rotating frame)
pos = [0 -0.1 0];
x = pos(1);
y = pos(2);
z = pos(3);
vel = 4.24*[cos(-10*pi/180) sin(-10*pi/180) 0]/1.025;  % 26000 kmph ---> 7.22        2.62 giving elliptical orbit for satellite
velx = vel(1);
vely = vel(2);
velz = vel(3);

X = [major1(1) major2(1)];
Y = [major1(2) major2(2)];
Z = [major1(3) major2(3)];

% figure(1);
% %plot3(X, Y, Z, '.r','MarkerSize', 15);
% plot(X, Y, '.r','MarkerSize', 10);
% xlabel('X')
% ylabel('Y')
% %zlabel('z')
% title('Rotating frame')
% hold on;

%% Numerical integration
pos_i = NB*pos';
for i=1:epoch
    
    %%%%%%%%%%%%%%%%%%%%%%% Rotating Frame %%%%%%%%%%%%%%%%%%%
    r1 = pos - major1;
    r1 = sqrt(sum(r1.*r1));
    r2 = pos - major2;
    r2 = sqrt(sum(r2.*r2));
    
    state0 = [x, y, z, velx, vely, velz];
    [time, state] = ode45(@statediffeq, [0,dt], state0);
    x = state(end,1);
    y = state(end,2);
    z = state(end,3);
    velx = state(end,4);
    vely = state(end,5);
    velz = state(end,6);
    pos = [x y z];
    
    NB = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%% Inertial Frame %%%%%%%%%%%%%%%%%%%
    major1_i = NB*major1';
    major2_i = NB*major2';
    
    pos_i = NB*pos';
%     x_i = pos_i(1);
%     y_i = pos_i(2);
%     z_i = pos_i(3);
    vel_i = NB*[velx-y; vely+x; velz];
    velx_i = vel_i(1);
    vely_i = vel_i(2);
    velz_i = vel_i(3);
    sqrt(sum(vel_i.*vel_i))
    %pos_i = pos_i + vel_i*dt;
    x_i = pos_i(1);
    y_i = pos_i(2);
    z_i = pos_i(3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Inertial frame 3D Plotting %%%%%%%%%%%%%%%%%%%
    
    figure(2)
    plot3([major1_i(1) major2_i(1)], [major1_i(2) major2_i(2)], [0 0], '.r','MarkerSize', 15)
    hold on;
    plot3(x_i, y_i, z_i, '.g');
    hold on;
    title('Inertial Frame')
    xlim([-1 1]);
    ylim([-1 1]);
    pause(0.0008);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Rotation frame 2D Plotting %%%%%%%%%%%%%%%%%%%
    
    figure(1)
    %%% 3D plot
    %plot3(x, y, z, '.b');
    %%% 2D plot
    plot(x, y, '.b');
    hold on;
    title('Rotating Frame')
    xlim([-5 5]);
    ylim([-20 20]);
    %zlim([-1 1]);
    pause(0.0008);
    
    t = t + da;
    
end

% state0 = [x, y, z, velx, vely, velz];
% [time, state] = ode45(@statediffeq, [0,1000*dt], state0);
% plot(state(:,1), state(:,2))
% Wierd results for epoch = 2000, pos = [-0.25 0 0], mu = 1.215*10^(-3), dt = 2.661*10^(-6)*3600, vel = 2.552*[cos(-90*pi/180) sin(-90*pi/180) 0]/1.025