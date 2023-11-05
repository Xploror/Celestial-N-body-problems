clc;clear

dt = 0.01;

G = 2*10^(-6);
M = 1000000;
m_s = 800;
m_sat = 1;
transfer = 0;  % Bool

resonance_ratio = [2,1];
e_s = 0.3;    %(Hyperparameter)
e_sat = 0.3;  %(Hyperparameter)
a_s = 3;      %(Hyperparameter)
T_s = ((4*pi^2/(G*M))*a_s^3)^(1/2);
T_sat = resonance_ratio(2)*T_s/resonance_ratio(1);
a_sat_target = a_s*(T_sat/T_s)^(2/3);

r_s_init = [a_s*(1+e_s),0];
v_init = [0,sqrt(G*M*((2/r_s_init(1))-(1/a_s)))];

%%%%% Transfer mechanism (Starts with an arbitrary parking orbit around primary 1)
e1_sat = 0; % Circular orbit (Hyperparameter)
a1_sat = a_sat_target*(1-e_sat); % (Hyperparameter)
T1_sat = ((4*pi^2/(G*M))*a1_sat^3)^(1/2);
r_sat_init = [a1_sat*(1+e1_sat),0];
v_sat_init = sqrt(G*M*((2/r_sat_init(1))-(1/a1_sat)))*cross([0,0,1], [r_sat_init,0])*[1 0 0; 0 1 0]'/norm(r_sat_init);
t_trans = 7*T1_sat/2;  % (Hyperparameter)

%Dynamic propogation
r_s = r_s_init;
v_s = v_init;
r_sat = r_sat_init;
v_sat = v_sat_init;
t = 0;

for i=0:dt:500
    acc_s = -((G*M/norm(r_s)^3)*r_s + (G*m_sat/norm(r_sat-r_s)^3)*(r_s - r_sat));
    acc_sat = -((G*M/norm(r_sat)^3)*r_sat + (G*m_s/norm(r_sat-r_s)^3)*(r_sat - r_s));
    v_s = v_s + acc_s*dt;
    v_sat = v_sat + acc_sat*dt;
    r_s = r_s + v_s*dt;
    r_sat = r_sat + v_sat*dt;
    
    if t >= t_trans && transfer == 0
        a_sat = (norm(r_sat) + a_sat_target*(1+e_sat))/2;
        v_sat = sqrt(G*M*((2/norm(r_sat))-(1/a_sat_target)))*cross([0,0,1], [r_sat,0])*[1 0 0; 0 1 0]'/norm(r_sat_init);  % This line of code is analogous to delta_v thrust point
        transfer = 1;
    end
    
    figure(1)
    plot(r_s(1), r_s(2), '.b', r_sat(1), r_sat(2), '.g', 0,0, '.b', 'MarkerSize', 20);
    xlim([-5,5])
    ylim([-5,5])
    
    if abs((r_s/norm(r_s))*(r_sat/norm(r_sat))' - 1) <= 0.0005
        %%%%% Plotting the spatial positions where body2 and satellite has
        %%%%% same direction of its radial vector
        figure(2)
        plot(r_s(1), r_s(2), '.b', r_sat(1), r_sat(2), '.g')
        xlim([-5,5])
        ylim([-5,5])
        hold on;
    end
    pause(0.0000001)
    
    t = t + dt;
end
