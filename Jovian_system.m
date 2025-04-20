clc;clear

dt = 10/1000;
transfer = 0;

%G = 2*10^(-6);
t_scaling = 24*3600;  %1 day in secs
G = 6.6743*(t_scaling^2)/(421700)^3;
M = 19000000; % Jupiter  1.9*10^27
m_s1 = 893;  % Io  8.93*10^22
m_s2 = 480;  % Europa   4.8*10^22
m_s3 = 1480;  % Ganymede  14.8*10^22
m_sat = 200*1e-20;  %200kg

scaling_jupiter_plot = 3e-4;

resonance_ratio = [2,1];
e_s1 = 0.0041;  % Io
e_s2 = 0.009;   % Europa
e_s3 = 0.0013;  % Ganymede
e_sat = 0.1;
a_s1 = 1;  % Io 421700
a_s2 = 670900/421700;   % Europa  670900
a_s3 = 1070400/421700;   % Ganymede 1070400
T_s1 = ((4*pi^2/(G*M))*a_s1^3)^(1/2);
T_s2 = ((4*pi^2/(G*M))*a_s2^3)^(1/2);
T_s3 = ((4*pi^2/(G*M))*a_s3^3)^(1/2);

% Transfer target
T_sat = resonance_ratio(2)*T_s3/resonance_ratio(1);
a_sat_target = a_s3*(T_sat/T_s3)^(2/3);

%%%%% Moon orbit pre-plotting
figure(1);
ang=0:0.01:2*pi;
io_orbit = a_s1*[cos(ang); sin(ang)];
europa_orbit = a_s2*[cos(ang); sin(ang)];
ganymede_orbit = a_s3*[cos(ang); sin(ang)];


r_s1_init = [a_s1*(1+e_s1),0];
r_s2_init = [a_s2*(1+e_s2),0];
r_s3_init = [a_s3*(1+e_s3),0];
v1_init = [0,sqrt(G*M*((2/r_s1_init(1))-(1/a_s1)))];
v2_init = [0,sqrt(G*M*((2/r_s2_init(1))-(1/a_s2)))];
v3_init = [0,sqrt(G*M*((2/r_s3_init(1))-(1/a_s3)))];

%%%%% Transfer mechanism (Starts with an arbitrary parking orbit around primary 1)
e1_sat = 0; % Circular orbit (Hyperparameter)
a1_sat = a_s1*(1-e_sat); % (Hyperparameter)
T1_sat = ((4*pi^2/(G*M))*a1_sat^3)^(1/2);
r_sat_init = [a1_sat*(1+e1_sat),0];
v_sat_init = sqrt(G*M*((2/norm(r_sat_init))-(1/a1_sat)))*cross([0,0,1], [r_sat_init,0])*[1 0 0; 0 1 0]'/norm(r_sat_init);
t_trans = 1*T1_sat;  % (Hyperparameter)

%Dynamic propogation
r_s1 = r_s1_init;
r_s2 = r_s2_init;
r_s3 = r_s3_init;
v_s1 = v1_init;
v_s2 = v2_init;
v_s3 = v3_init;
r_sat = r_sat_init;
v_sat = v_sat_init;
t = 0;

% Initial Plotting
figure(1);
plot(io_orbit(1,:),io_orbit(2,:),'--r', europa_orbit(1,:),europa_orbit(2,:),'--b', ganymede_orbit(1,:),ganymede_orbit(2,:),'--g');
title('2D Jovian System')
xlabel('X (dimensionless)');ylabel('Y (dimensionless)');
grid on; hold on;
sat_traj = [];

for i=0:dt:10
    acc_s1 = -((G*M/norm(r_s1)^3)*r_s1 + (G*m_sat/norm(r_sat-r_s1)^3)*(r_s1 - r_sat) + (G*m_s2/norm(r_s2-r_s1)^3)*(r_s1 - r_s2) + (G*m_s3/norm(r_s3-r_s1)^3)*(r_s1 - r_s3));
    acc_s2 = -((G*M/norm(r_s2)^3)*r_s2 + (G*m_sat/norm(r_sat-r_s2)^3)*(r_s2 - r_sat) + (G*m_s1/norm(r_s2-r_s1)^3)*(r_s2 - r_s1) + (G*m_s3/norm(r_s3-r_s2)^3)*(r_s2 - r_s3));
    acc_s3 = -((G*M/norm(r_s3)^3)*r_s3 + (G*m_sat/norm(r_sat-r_s3)^3)*(r_s3 - r_sat) + (G*m_s2/norm(r_s2-r_s3)^3)*(r_s3 - r_s2) + (G*m_s1/norm(r_s3-r_s1)^3)*(r_s3 - r_s1));
    acc_sat = -((G*M/norm(r_sat)^3)*r_sat + (G*m_s1/norm(r_sat-r_s1)^3)*(r_sat - r_s1) + (G*m_s2/norm(r_sat-r_s2)^3)*(r_sat - r_s2) + (G*m_s3/norm(r_sat-r_s3)^3)*(r_sat - r_s3));
    v_s1 = v_s1 + acc_s1*dt;
    v_s2 = v_s2 + acc_s2*dt;
    v_s3 = v_s3 + acc_s3*dt;
    v_sat = v_sat + acc_sat*dt;
    r_s1 = r_s1 + v_s1*dt;
    r_s2 = r_s2 + v_s2*dt;
    r_s3 = r_s3 + v_s3*dt;
    r_sat = r_sat + v_sat*dt;
    
    %%%%%%%%%%%%% Transfer mechanism %%%%%%%%%%%%%%%
    if t >= t_trans && transfer == 0
        a_sat = (norm(r_sat) + a_sat_target*(1+e_sat))/2;
        v_sat = sqrt(G*M*((2/norm(r_sat))-(1/a_sat_target)))*cross([0,0,1], [r_sat,0])*[1 0 0; 0 1 0]'/norm(r_sat_init);  % This line of code is analogous to delta_v thrust point
        transfer = 1;
        fprintf('Transfer completed!')
        figure(1)
        plot(r_sat(1), r_sat(2), '*k');hold on;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sat_traj = [sat_traj; r_sat];
    figure(1);
    ppp = plot(r_s1(1), r_s1(2), '.r', r_s2(1), r_s2(2), '.b', r_s3(1), r_s3(2), '.g', 0,0, '.k', r_sat(1), r_sat(2), '.y', 'MarkerSize', 20);hold on;
    plot(sat_traj(:,1), sat_traj(:,2), '--y');
    xlim([-3,3])
    ylim([-3,3])

    % Contour plot
    x = -3:0.05:3;
    y = -3:0.05:3;
    [X, Y] = meshgrid(x,y);
    len_x = length(X);
    len_y = length(Y);
    U = -(scaling_jupiter_plot*G*M./sqrt(X.^2 + Y.^2) + G*m_s1./sqrt((X-r_s1(1)*ones(len_x,len_x)).^2 + (Y-r_s1(2)*ones(len_y,len_y)).^2) + G*m_s2./sqrt((X-r_s2(1)*ones(len_x,len_x)).^2 + (Y-r_s2(2)*ones(len_y,len_y)).^2) + G*m_s3./sqrt((X-r_s3(1)*ones(len_x,len_x)).^2 + (Y-r_s3(2)*ones(len_y,len_y)).^2));
    figure(2)
    contourf(X, Y, U);
    hold on;
    
    R = sqrt(X.^2+Y.^2);
    acc_obj_x = -((scaling_jupiter_plot*G*M./R^3).*X + (G*m_s1./((X-r_s1(1)*ones(len_x,len_x)).^2 + (Y-r_s1(2)*ones(len_y,len_y)).^2)^(3/2)).*(X - r_s1(1)*ones(len_x,len_x)) + (G*m_s2./((X-r_s2(1)*ones(len_x,len_x)).^2 + (Y-r_s2(2)*ones(len_y,len_y)).^2)^(3/2)).*(X - r_s2(1)*ones(len_x,len_x)) + (G*m_s3./((X-r_s3(1)*ones(len_x,len_x)).^2 + (Y-r_s3(2)*ones(len_y,len_y)).^2)^(3/2)).*(X - r_s3(1)*ones(len_x,len_x)));
    acc_obj_y = -((scaling_jupiter_plot*G*M./R^3).*Y + (G*m_s1./((X-r_s1(1)*ones(len_x,len_x)).^2 + (Y-r_s1(2)*ones(len_y,len_y)).^2)^(3/2)).*(Y - r_s1(2)*ones(len_y,len_y)) + (G*m_s2./((X-r_s2(1)*ones(len_x,len_x)).^2 + (Y-r_s2(2)*ones(len_y,len_y)).^2)^(3/2)).*(Y - r_s2(2)*ones(len_y,len_y)) + (G*m_s3./((X-r_s3(1)*ones(len_x,len_x)).^2 + (Y-r_s3(2)*ones(len_y,len_y)).^2)^(3/2)).*(Y - r_s3(2)*ones(len_y,len_y)));
    quiver(X, Y, real(acc_obj_x), real(acc_obj_y));
    title('Contour Plot (Potential Field)');hold on;
    xlabel('X (dimensionless)');ylabel('Y (dimensionless)');
    hold off;
    pause(0.0000001)
    
    t = t + dt;
    
    delete(ppp)
end

