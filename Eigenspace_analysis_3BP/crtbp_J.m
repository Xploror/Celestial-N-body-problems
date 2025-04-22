function A = crtbp_J(x, y, z, mu)

%mu = 1.215*10^(-2);

%%%%%%% Jacobian elements %%%%%%%
f41 = 1 - (1-mu)/sqrt((x+mu)^2+y^2+z^2)^3 - mu/sqrt((x-(1-mu))^2+y^2+z^2)^3 + 3*(1-mu)*(x+mu)^2/sqrt((x+mu)^2+y^2+z^2)^5 + 3*mu*(x-(1-mu))^2/sqrt((x-(1-mu))^2+y^2+z^2)^5;
f52 = 1 - (1-mu)/sqrt((x+mu)^2+y^2+z^2)^3 - mu/sqrt((x-(1-mu))^2+y^2+z^2)^3 + 3*(1-mu)*y^2/sqrt((x+mu)^2+y^2+z^2)^5 + 3*mu*y^2/sqrt((x-(1-mu))^2+y^2+z^2)^5;
f63 = (1-mu)/sqrt((x+mu)^2+y^2+z^2)^3 - mu/sqrt((x-(1-mu))^2+y^2+z^2)^3 + 3*(1-mu)*z^2/sqrt((x+mu)^2+y^2+z^2)^5 + 3*mu*z^2/sqrt((x-(1-mu))^2+y^2+z^2)^5;
f42 = 3*(1-mu)*(x+mu)*y/sqrt((x+mu)^2+y^2+z^2)^5 + 3*mu*(x-(1-mu))*y/sqrt((x-(1-mu))^2+y^2+z^2)^5;
f51 = f42;
f43 = 3*(1-mu)*(x+mu)*z/sqrt((x+mu)^2+y^2+z^2)^5 + 3*mu*(x-(1-mu))*z/sqrt((x-(1-mu))^2+y^2+z^2)^5;
f61 = f43;
f53 = 3*(1-mu)*z*y/sqrt((x+mu)^2+y^2+z^2)^5 + 3*mu*y*z/sqrt((x-(1-mu))^2+y^2+z^2)^5;
f62 = f53;

A = [0 0 0 1 0 0;
    0 0 0 0 1 0
    0 0 0 0 0 1
    f41 f42 f43 0 2 0
    f51 f52 f53 -2 0 0
    f61 f62 f63 0 0 0];

end