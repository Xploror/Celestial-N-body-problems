function y = PCC4BP_eqn(t, x, mu, mu_bar, x1t, omega, perturber)
% Planar Concentric circular resticted 4 Body Problem
% Important x1t is the distance vector of third primary wrt first primary (the massive one)
p = [-mu+x1t(1)*cos(x(5)) x1t(1)*sin(x(5))];  % Third primary's radial vector with origin of rotating frame

if perturber == 3
    r1 = x(1:2)' - [-mu 0];
    r1 = sqrt(sum(r1.*r1));
    r2 = x(1:2)' - [(1-mu) 0];
    r2 = sqrt(sum(r2.*r2));
    r3 = x(1:2)' - p;
    r3 = sqrt(sum(r3.*r3));
    
    r13 = norm(x1t);

    y = zeros(length(x),1);
    y(1:2) = x(3:4);
    y(3) = 2*x(4) + x(1) - (1-mu)*(x(1)+mu)/r1^3 - mu*(x(1)-1+mu)/r2^3 - mu_bar*(x(1)-p(1))/r3^3 - mu_bar*cos(x(5))/r13^2;
    y(4) = x(2) - 2*x(3) - (1-mu)*x(2)/r1^3 - mu*x(2)/r2^3 - mu_bar*(x(2)-p(2))/r3^3 - mu_bar*sin(x(5))/r13^2;
    y(5) = omega-1;
    
else
    r1 = x(1:2)' - [-mu 0];
    r1 = sqrt(sum(r1.*r1));
    r3 = x(1:2)' - [(1-mu) 0];
    r3 = sqrt(sum(r3.*r3));
    r2 = x(1:2)' - p';
    r2 = sqrt(sum(r2.*r2));
    
    r12 = norm(x1t);

    y = zeros(length(x),1);
    y(1:2) = x(3:4);
    y(3) = 2*x(4) + x(1) - (1-mu)*(x(1)+mu)/r1^3 - mu*(x(1)-1+mu)/r3^3 - mu_bar*(x(1)-p(1))/r2^3 - mu_bar*cos(x(5))/r12^2;
    y(4) = x(2) - 2*x(3) - (1-mu)*x(2)/r1^3 - mu*x(2)/r3^3 - mu_bar*(x(2)-p(2))/r2^3 - mu_bar*sin(x(5))/r12^2;
    y(5) = omega-1;
end


end