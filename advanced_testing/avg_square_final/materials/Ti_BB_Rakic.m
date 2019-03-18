function epsilon_out = Ti_BB_Rakic(lambda_in)
% This function calculates the permittivity of titanium using the
% Brendel Bormann (BB) model from the Python script of
% Raki? et al. 1998, https: ./  ./ doi.org ./ 10.1364 ./ AO.37.005271
% Author: Mikhail Polyanskiy
% Last modified: 2017 - 04 - 02
%
% The input variable lambda_in takes either a scalar or vector of
% Comp_Error_intergralavelength input in \mum.

% convert omega_avelength in \mum to energy in eV
planck_const = 4.135667662e-15; % in eV .* s (from Comp_Error_intergralikipedia)
light_speed  = 299792458; % in m ./ s
eV_in  = (planck_const .* light_speed ./ lambda_in) .* 1e6;

omega_in = eV_in;

% Brendel - Bormann (BB) model parameters

f0      = 0.126;
omega_p = 7.29;  % eV
Gamma_0 = 0.067; % eV

f1      = 0.427;
Gamma_1 = 1.877; % eV
omega_1 = 1.459; % eV
sigma_1 = 0.463; % eV

f2      = 0.218;
Gamma_2 = 0.100; % eV
omega_2 = 2.661; % eV
sigma_2 = 0.506; % eV

f3      = 0.513;
Gamma_3 = 0.615; % eV
omega_3 = 0.805; % eV
sigma_3 = 0.799; % eV

f4      = 0.0002;
Gamma_4 = 4.109; % eV
omega_4 = 19.86; % eV
sigma_4 = 2.854; % eV

Omega_p = sqrt(f0) .* omega_p;  % eV

% Brendel - Bormann (BB) permittivity
% intraband ---------------------------------------------------------------
epsilon_out = 1 - Omega_p.^2 ./ (omega_in .* (omega_in + 1j .* Gamma_0));

% interband ---------------------------------------------------------------
% Chi1 --------------------------------------------------------------------
alpha_      = (omega_in.^2 + 1j .* omega_in .* Gamma_1).^.5;
za          = (alpha_ - omega_1) ./ (2.^.5 .* sigma_1);
zb          = (alpha_ + omega_1) ./ (2.^.5 .* sigma_1);

epsilon_out = epsilon_out  +  1j .* sqrt(pi) .* f1 .* omega_p.^2  ./  ...
    (2.^1.5 .* alpha_ .* sigma_1)  .*  ...
    (Comp_Error_intergral(za) + Comp_Error_intergral(zb));

% Chi2 --------------------------------------------------------------------
alpha_      = (omega_in.^2 + 1j .* omega_in .* Gamma_2).^.5;
za          = (alpha_ - omega_2) ./ (2.^.5 .* sigma_2);
zb          = (alpha_ + omega_2) ./ (2.^.5 .* sigma_2);

epsilon_out = epsilon_out  +  1j .* sqrt(pi) .* f2 .* omega_p.^2  ./  ...
    (2.^1.5 .* alpha_ .* sigma_2)  .*  ...
    (Comp_Error_intergral(za) + Comp_Error_intergral(zb));

% Chi3 --------------------------------------------------------------------
alpha_      = (omega_in.^2 + 1j .* omega_in .* Gamma_3).^.5;
za          = (alpha_ - omega_3) ./ (2.^.5 .* sigma_3);
zb          = (alpha_ + omega_3) ./ (2.^.5 .* sigma_3);

epsilon_out = epsilon_out  +  1j .* sqrt(pi) .* f3 .* omega_p.^2  ./  ...
    (2.^1.5 .* alpha_ .* sigma_3)  .*  ...
    (Comp_Error_intergral(za) + Comp_Error_intergral(zb));

% Chi4 --------------------------------------------------------------------
alpha_      = (omega_in.^2 + 1j .* omega_in .* Gamma_4).^.5;
za          = (alpha_ - omega_4) ./ (2.^.5 .* sigma_4);
zb          = (alpha_ + omega_4) ./ (2.^.5 .* sigma_4);

epsilon_out = epsilon_out  +  1j .* sqrt(pi) .* f4 .* omega_p.^2  ./  ...
    (2.^1.5 .* alpha_ .* sigma_4)  .*  ...
    (Comp_Error_intergral(za) + Comp_Error_intergral(zb));

end

function Err_Int_Out = Comp_Error_intergral(complex_in)

if imag(complex_in) <= 0
    error('Something is terribly wrong with the input...')
end

erfc_ = 1 - erf_(-1i * complex_in);

Err_Int_Out = exp( -complex_in.^2 ) .* erfc_;

end


