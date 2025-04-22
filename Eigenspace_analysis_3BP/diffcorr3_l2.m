function [x0_corr, t_corr, phi] = diffcorr3_l2(x0, mu, tStep)

% Jacobian of dynamics x_dot = f(x) (2D xy plane) for L1
t = 0; % time-count
phi = eye(length(x0));
x0_corr = x0;
x = x0_corr;
x_prev = x;
f = @(x) statediffeq(1, x, mu);
phidot = @(A, phi) A*phi;
count = 0;  % Counter used to avoiding the initialized small values of x(2) and x(3)
k = 0; % Counter used for checking stagnation in x(3) value
x4_prev = 0;

% Integrate states and find monodromy matrix from t0 to t1 taking original x0
while abs(x(4))>1*10^(-6) || count<2
    x(4)
    % preparing init conditions for next iteration by updating with STM matrix...
    t = 0; % time-count
    del_xt = [0;0-x(2);0;0-x(4);0;0];
    del_x0_req = inv(phi)*del_xt;
    if x(1)>(1-mu)
        x0_corr = x0_corr + [0;0;0;0;0.0001*del_x0_req(5);0];
        %x0_corr(4) = 1.00005*x0_corr(4);
    else
        x0_corr(5) = 0.95*x0_corr(5);
    end
    x = x0_corr;
    X = [x0];
    phi = eye(length(x0)); % State Transition Matrix (Monodromy matrix)
    while sign(x(2))*sign(x_prev(2))==1 || t<100*tStep  % tolerance (position crossing x-axis perpendicularly) & dont stop before 1 iterations of tStep
        x_prev = x; %(previous iteration states)
        
        % Runge Kutta 6th Order Numerical Integration for x & phi individually
        A = crtbp_J(x(1),x(2),x(3),mu);
        
        %%% 6th Order RK method for integration
        k1 = tStep*phidot(A, phi);
        k2 = tStep*phidot(A, phi+k1);
        k3 = tStep*phidot(A, phi+(3*k1+k2)/8);
        k4 = tStep*phidot(A, phi+(8*k1+2*k2+8*k3)/27);
        k5 = tStep*phidot(A, phi+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
        k6 = tStep*phidot(A, phi+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
        k7 = tStep*phidot(A, phi+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
        phi = phi + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
        
        k1 = tStep*f(x);
        k2 = tStep*f(x+k1);
        k3 = tStep*f(x+(3*k1+k2)/8);
        k4 = tStep*f(x+(8*k1+2*k2+8*k3)/27);
        k5 = tStep*f(x+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
        k6 = tStep*f(x+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
        k7 = tStep*f(x+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
        x = x + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;  
        X = [X x];

        t = t + tStep;
    end
    % Plotting for debugging
%     plot3(X(1,:),X(2,:), X(3,:), '*r', 1-mu, 0, 0, 'og')
%      %hold on;
%     yline([0])
%     pause(0.1)
   
    if abs(x4_prev-x(4))<10^(-8) && x(4)<10^(-4)
        k = k + 1;
    else
        x4_prev = x(4);
        k = 0;
    end
    if k>30
        break;
    end
    
    t_corr = t;
    count = count + 1;
    
end
end