function epsilon_out = Ti_LD_Rakic(lambda_in)
% This function calculates the permittivity of titanium using the
% Lorentz-Drude (LD) model from the Python script of
% Raki? et al. 1998, https://doi.org/10.1364/AO.37.005271
% Author: Mikhail Polyanskiy
% Last modified: 2017-04-02
%
% The input variable lambda_in takes either a scalar or vector of
% wavelength input in \mum.

% convert omega_avelength in \mum to energy in eV
planck_const = 4.135667662e-15; % in eV*s (from Wikipedia)
light_speed  = 299792458; % in m/s
eV_in  = (planck_const * light_speed ./ lambda_in) * 1e6;

omega_in = eV_in;

% Lorentz-Drude (LD) model parameters
f0      = 0.148;
Gamma_0 = 0.082; % eV
omega_p = 7.29;  % eV

f1      = 0.899;
Gamma_1 = 2.276; % eV
omega_1 = 0.777; % eV

f2      = 0.393;
Gamma_2 = 2.518; % eV
omega_2 = 1.545; % eV

f3      = 0.187;
Gamma_3 = 1.663; % eV
omega_3 = 2.509; % eV

f4      = 0.001;
Gamma_4 = 1.762; % eV
omega_4 = 19.43; % eV

Omega_p = sqrt(f0) * omega_p; % eV

% Lorentz-Drude permittivity
epsilon_out = 1 - Omega_p^2 ./ (omega_in .* (omega_in + 1j * Gamma_0) ) ...
    + f1 * omega_p^2 ./ ( (omega_1^2 - omega_in^2) - 1j * omega_in .* Gamma_1) ...
    + f2 * omega_p^2 ./ ( (omega_2^2 - omega_in^2) - 1j * omega_in .* Gamma_2) ...
    + f3 * omega_p^2 ./ ( (omega_3^2 - omega_in^2) - 1j * omega_in .* Gamma_3) ...
    + f4 * omega_p^2 ./ ( (omega_4^2 - omega_in^2) - 1j * omega_in .* Gamma_4);