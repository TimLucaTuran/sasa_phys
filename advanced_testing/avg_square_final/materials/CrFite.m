function eps_fun = CrFite(LAMBDA_)

% Calculate epsilon of chrome for a given lambda using a fit function
% supplied by Kay Dietrich.

% LAMBDA_=1./linspace(1/3,1/0.2,101);

epsinf = 1;
delta_1 = 497.9688;
a_1 = 1;
b_1 = 0.238032;
c_1 = 0;
delta_2 = 447.5791;
a_2 = 1;
b_2 = 16.07984;
c_2 = 0.375532;
delta_3 = 444.615;
a_3 = 1;
b_3 = 6.609194;
c_3 = 7.562677;
delta_4 = 3405.751;
a_4 = 1;
b_4 = 13.55265;
c_4 = 99.54246;
% delta_5 = 2445.382*0;
% a_5 = 1;
% b_5 = 6.761129;
% c_5 = 1975.014;

omega=2*pi./LAMBDA_;

%%% Der Fit ohne die hohen Resonanzfrequenzen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measured by Kay
eps_fun = epsinf + delta_1./(-a_1*omega.^2 - 1i*b_1*omega+c_1)+...
    delta_2./(-a_2*omega.^2 - 1i*b_2*omega + c_2)+...
    delta_3./(-a_3*omega.^2 - 1i*b_3*omega + c_3)+...
    delta_4./(-a_4*omega.^2 - 1i*b_4*omega + c_4);%+...
%     delta_5./(-a_5*omega.^2-1i*b_5*omega+c_5);



