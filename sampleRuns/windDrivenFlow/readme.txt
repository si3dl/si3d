Model run of a rectangular basin with a constant wind forcing. 
This model run is for the purpose of validating si3d against theoretical solutions for the velocity profile induced by wind in a recatangular basin. 

Details of the test case can be read in the si3d_inp.txt file used for the model run. 

The analytical solution is:

%% Analytical solution Huang 1993
k = 0.005;      % [m/s] Linearized friction coefficient
dl_dx = (3/2)*(tau/(rho_w*g*h))*((2*K_v+k*h)/(3*K_v+k*h));
c1 = (1/(6*K_v))*g*h^2;
c2 = (tau*h/(2*rho_w*K_v));
u_Huang = 100*(c1.*dl_dx.*(3.*(dh-1).^2-1)+c2.*(2.*dh-1));