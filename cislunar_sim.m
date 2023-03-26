% CisLunar Numerical Simulation and Plotting
clc;clear;

%% Variables

epoch = 2000;
mu = 1.215*10^(-3);  % mass ratio
major1 = [-mu 0 0];
major2 = [(1-mu) 0 0];
t = 0;
dt = 2.661*10^(-6)*3600;  % taking 1 hr time interval for Earth-Moon system
da = 2.661*10^(-6)*3600;
%dt = 0.001;  % Time interval 0.1 was working fine
%da = 0.001; % angle change
NB = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

%% Initial Conditions (Rotating frame)
pos = [-0.25 0 0];
x = pos(1);
y = pos(2);
z = pos(3);
vel = 2.52*[cos(-90*pi/180) sin(-90*pi/180) 0]/1.025;  % 26000 kmph ---> 7.22        2.62 giving elliptical orbit for satellite
velx = vel(1);
vely = vel(2);
velz = vel(3);

X = [major1(1) major2(1)];
Y = [major1(2) major2(2)];
Z = [major1(3) major2(3)];

figure(1);
%plot3(X, Y, Z, '.r','MarkerSize', 15);
plot(X, Y, '.r','MarkerSize', 10);
xlabel('X')
ylabel('Y')
%zlabel('z')
title('Rotating frame')
hold on;

%% Numerical integration
pos_i = NB*pos';
for i=1:epoch
    
    %%%%%%%%%%%%%%%%%%%%%%% Rotating Frame %%%%%%%%%%%%%%%%%%%
    r1 = pos - major1;
    r1 = sqrt(sum(r1.*r1));
    r2 = pos - major2;
    r2 = sqrt(sum(r2.*r2));
    
    dvx = (2*vely + x - ((1-mu)*(x + mu)/(r1^3)) - mu*(x - (1-mu))/(r2^3))*dt;
    dvy = (y - 2*velx - ((1-mu)*y/(r1^3)) - mu*y/(r2^3))*dt;
    dvz = - (((1-mu)*z/r1^3) - mu*z/r2^3)*dt;
    
    velx = velx + dvx;
    vely = vely + dvy;
    velz = velz + dvz;
    
    x = x + velx*dt;
    y = y + vely*dt;
    z = z + velz*dt;
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
    %pos_i = pos_i + vel_i*dt;
    x_i = pos_i(1);
    y_i = pos_i(2);
    z_i = pos_i(3);
    
    figure(2)
    plot([major1_i(1) major2_i(1)], [major1_i(2) major2_i(2)], '.r','MarkerSize', 15)
    hold on;
    plot(x_i, y_i, '.g');
    hold on;
    title('Inertial Frame')
    %xlim([-15 15]);
    %ylim([-5 5]);
    %pause(0.0008);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(1)
    %plot3(x, y, z, '.b');
    plot(x, y, '.b');
    hold on;
    xlim([-5 5]);
    ylim([-20 20]);
    %zlim([-1 1]);
    %pause(0.0000000001)
    
    t = t + da;
    
end

% Wierd results for epoch = 2000, pos = [-0.25 0 0], mu = 1.215*10^(-3), dt = 2.661*10^(-6)*3600, vel = 2.552*[cos(-90*pi/180) sin(-90*pi/180) 0]/1.025