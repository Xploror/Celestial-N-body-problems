function y = statediffeq(t, st0)
%%% st0(1) --> x | st0(2) --> y | st0(3) --> z | st0(4) --> vx | st0(5) --> vy | st0(6) --> vz
mu = 1.215*10^(-2);
r1 = st0(1:3)' - [-mu 0 0];
r1 = sqrt(sum(r1.*r1));
r2 = st0(1:3)' - [(1-mu) 0 0];
r2 = sqrt(sum(r2.*r2));

y = zeros(length(st0),1);
y(1:3) = st0(4:end);
y(4) = 2*st0(5) + st0(1) - (1-mu)*(st0(1)+mu)/r1^3 - mu*(st0(1)-1+mu)/r2^3;
y(5) = st0(2) - 2*st0(4) - (1-mu)*st0(2)/r1^3 - mu*st0(2)/r2^3;
y(6) = -(1-mu)*st0(3)/r1^3 - mu*st0(3)/r2^3;
 
end