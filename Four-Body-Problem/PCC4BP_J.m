function A = PCC4BP_J(x, y, theta, mu, mu1, x1t)
% Note: x1t is the radial vector between first and third primary!!!

r1t = norm(x1t);
p = [-mu+x1t(1)*cos(theta) x1t(1)*sin(theta)]; % Radial distance of perturbing primary wrt origin of the rotating frame 

%%%%%%% Jacobian elements %%%%%%%
f41 = 1 - (1-mu)/sqrt((x+mu)^2+y^2)^3 - mu/sqrt((x-(1-mu))^2+y^2)^3 -mu1/norm([x y]-p)^3 + 3*(1-mu)*(x+mu)^2/sqrt((x+mu)^2+y^2)^5 + 3*mu*(x-(1-mu))^2/sqrt((x-(1-mu))^2+y^2)^5 + 3*mu1*((x-p(1))^2)/norm([x y]-p)^5;
f52 = 1 - (1-mu)/sqrt((x+mu)^2+y^2)^3 - mu/sqrt((x-(1-mu))^2+y^2)^3 -mu1/norm([x y]-p)^3 + 3*(1-mu)*y^2/sqrt((x+mu)^2+y^2)^5 + 3*mu*y^2/sqrt((x-(1-mu))^2+y^2)^5 + 3*mu1*((y-p(2))^2)/norm([x y]-p)^5;
f42 = 3*(1-mu)*(x+mu)*y/sqrt((x+mu)^2+y^2)^5 + 3*mu*(x-(1-mu))*y/sqrt((x-(1-mu))^2+y^2)^5 + 3*mu1*((x-p(1))*(y-p(2)))/norm([x y]-p)^5;
f51 = 3*(1-mu)*(x+mu)*y/sqrt((x+mu)^2+y^2)^5 + 3*mu*(x-(1-mu))*y/sqrt((x-(1-mu))^2+y^2)^5 + 3*mu1*((y-p(2))*(x-p(1)))/norm([x y]-p)^5;
f43 = mu1*sin(theta)/r1t^2; %- mu1*norm(x1t)*sin(theta)/norm([x y]-p)^3 + 3*mu1*norm(x1t)*(-(x-p(1))*sin(theta) + (y-p(2))*cos(theta))/norm([x y]-p)^5;
f53 = -mu1*cos(theta)/r1t^2; %+ mu1*norm(x1t)*cos(theta)/norm([x y]-p)^3 + 3*mu1*norm(x1t)*(-(x-p(1))*sin(theta) + (y-p(2))*cos(theta))/norm([x y]-p)^5;

A = [0 0 1 0 0;
    0 0 0 1 0;
    f41 f42 0 2 f43; 
    f51 f52 -2 0 f53;
    0 0 0 0 0];

end