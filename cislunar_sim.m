% CisLunar Numerical Simulation and Plotting
clc;clear;

%% Variables

epoch = 1300;
mu = 1.215*10^(-2);  % Earth-Moon mass ratio 
major1 = [-mu 0 0];
major2 = [(1-mu) 0 0];
t = 0;
dt = 2*pi*10^(-6)*60*60/2.361;  % taking 5 min time interval for Earth-Moon system
da = 2*pi*10^(-6)*60*60/2.361; 
%dt = 0.001;  % Time interval 0.1 was working fine
%da = 0.001; % angle change
NB = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

%% Initial Conditions (Rotating frame)
pos = [-0.2 -0.1 0];
x = pos(1);
y = pos(2);
z = pos(3);
vel = 2.62*[cos(-90*pi/180) sin(-90*pi/180) 0]/1.025;  % 26000 kmph ---> 7.22        2.62 giving elliptical orbit for satellite    4.24 is a flyby (-10 deg and [0 -0.1 0])
velx = vel(1);
vely = vel(2);
velz = vel(3);

X = [major1(1) major2(1)];
Y = [major1(2) major2(2)];
Z = [major1(3) major2(3)];

vidobj_iner = VideoWriter("out_iner.jpg", "MPEG-4");
vidobj_rot = VideoWriter("out_rot.jpg", "MPEG-4");
vidobj_iner.FrameRate = 30;
vidobj_rot.FrameRate = 30;
open(vidobj_iner);
open(vidobj_rot);

figure(1);
ang=0:0.01:2*pi; 
xp=1*cos(ang);
yp=1*sin(ang);

figure(1);
plot(X, Y, '.r','MarkerSize', 15);hold on;
plot((-mu)+xp,yp,'--r');
xlabel('X')
ylabel('Y')
title('Rotating frame')
grid on;hold on;
f1 = getframe(gcf);
writeVideo(vidobj_rot, f1);

figure(2);
plot3((-mu)+xp,yp,zeros(length(xp),1),'--r');
grid on;hold on;
f2 = getframe(gcf);
writeVideo(vidobj_iner, f2);

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
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    grid on;
    f2 = getframe(gcf);
    writeVideo(vidobj_iner, f2);
    pause(0.0000001);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Rotation frame 2D Plotting %%%%%%%%%%%%%%%%%%%
    
    figure(1)
    %%% 3D plot
    %plot3(x, y, z, '.b');
    %%% 2D plot
    plot(x, y, '.b');
    hold on;
    title('Rotating Frame')
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    %zlim([-1 1]);
    pause(0.0000001);
    f1 = getframe(gcf);
    writeVideo(vidobj_rot, f1);
    pause(0.0000001);
    
    t = t + da;
    
end
vidobj_iner.close();
vidobj_rot.close();