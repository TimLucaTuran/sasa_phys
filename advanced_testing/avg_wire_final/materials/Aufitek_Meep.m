function efkt = Aufitek_Meep(lambda)
% Meep Au-fit from Kay
% Use ONLY in FMM code (for single wavelength in the loop)

epsinf = 5.53;
delta_1 = 2178.43;
a_1 = 1;
b_1 = 0.30978;
c_1 = 0;
delta_2 = 465.79;
a_2 = 1;
b_2 = 2.94869;
c_2 = 228.713;

% lambda=1./linspace(1/2.5,1/0.5,101);
omega = 2*pi./lambda;

%%% Der Fit ohne die hohen Resonanzfrequenzen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_kay = epsinf+delta_1./(-a_1*omega.^2-1i*b_1*omega+c_1)+...
                delta_2./(-a_2*omega.^2-1i*b_2*omega+c_2);
efkt = epsilon_kay;