function [t,X,STM_b] = libration_orbit_calculator(l2, del_x, mu1, theta)

tStep = 1*10^(-3);

A_3d = crtbp_J(l2, 0, 0, mu1);
A = [A_3d(1:2,1:2) A_3d(1:2,4:5); A_3d(4:5,1:2) A_3d(4:5,4:5)];
 
[l2_eigvec, l2_eigv] = eig(A);
l2_eigvec = real(l2_eigvec); % Cause these eigenvectors will be used in real space

T_p = 3*pi/abs(l2_eigv(3,3)); % Third eigenvalue of the system in the form +-i*nu
x0 = [l2;0;0;0] + del_x*l2_eigvec(:,3); 
x0 = [x0; theta];

[a,b, STM_b] = diffcorr_l2_FBP(x0, mu1, 0, [1.5,0], 0, tStep);
[t, X] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[1.5,0], 0, 3), 0:0.001:2*b, a);
X(:,5) = x0(5)*ones([length(X(:,5)),1]);


end