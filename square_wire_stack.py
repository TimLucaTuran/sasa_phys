import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from stack import MetaLayer, NonMetaLayer, Stack

# load refractive indices
def refractive_ind(wav_vec, material,fit_interpret):
    if fit_interpret == 1:
        if material == 1:
            # saphire
            index = np.sqrt(1+1.4313493/(1-(0.0726631/wav_vec)**2)
                            + 0.65054713/(1-(0.1193242/wav_vec)**2)
                            + 5.3414021/(1-(18.028251/wav_vec)**2))
        if material == 2:
            # silica
            index = np.sqrt(1+0.6961663/(1-(0.0684043/wav_vec)**2)
                            + 0.4079426/(1-(0.1162414/wav_vec)**2)
                            + 0.8974794/(1-(9.896161/wav_vec)**2))
        if material == 7:
            # Si4N4
            index = np.sqrt(1 + 3.0249/(1 - (0.1353406/wav_vec)**2)
                    + 40314/(1 - (1239.842/wav_vec)**2))
        if material == 12:
            # Ti from Rakic et al.
            eps_ti = Ti_LD_Rakic(wav_vec)
            index = np.sqrt(eps_ti)
    return index

def Ti_LD_Rakic(wav_vec):
    planck_const = 4.135667662e-15
    light_speed  = 299792458
    eV_in  = (planck_const * light_speed / wav_vec) * 10**6
    omega_in = eV_in
    # Lorentz-Drude (LD) model parameters
    f0      = 0.148
    Gamma_0 = 0.082
    omega_p = 7.29
    f1      = 0.899
    Gamma_1 = 2.276
    omega_1 = 0.777
    f2      = 0.393
    Gamma_2 = 2.518
    omega_2 = 1.545
    f3      = 0.187
    Gamma_3 = 1.663
    omega_3 = 2.509
    f4      = 0.001
    Gamma_4 = 1.762
    omega_4 = 19.43

    Omega_p = np.sqrt(f0) * omega_p

    # Lorentz-Drude permittivity
    epsilon_out = (1 - Omega_p**2 / (omega_in * (omega_in + 1j * Gamma_0) )
        + f1 * omega_p**2 / ( (omega_1**2 - omega_in**2) - 1j * omega_in * Gamma_1)
        + f2 * omega_p**2 / ( (omega_2**2 - omega_in**2) - 1j * omega_in * Gamma_2)
        + f3 * omega_p**2 / ( (omega_3**2 - omega_in**2) - 1j * omega_in * Gamma_3)
        + f4 * omega_p**2 / ( (omega_4**2 - omega_in**2) - 1j * omega_in * Gamma_4))
    return epsilon_out
"""
def Ti_BB_Rakic(wav_vec):
    planck_const = 4.135667662e-15
    light_speed  = 299792458
    eV_in  = (planck_const * light_speed / wav_vec) * 1e6
    omega_in = eV_in
    # Brendel - Bormann (BB) model parameters
    f0      = 0.126
    omega_p = 7.29
    Gamma_0 = 0.067
    f1      = 0.427
    Gamma_1 = 1.877
    omega_1 = 1.459
    sigma_1 = 0.463
    f2      = 0.218
    Gamma_2 = 0.100
    omega_2 = 2.661
    sigma_2 = 0.506
    f3      = 0.513
    Gamma_3 = 0.615
    omega_3 = 0.805
    sigma_3 = 0.799
    f4      = 0.0002
    Gamma_4 = 4.109
    omega_4 = 19.86
    sigma_4 = 2.854
    Omega_p = np.sqrt(f0) * omega_p

    # Brendel - Bormann (BB) permittivity
    epsilon_out = 1 - Omega_p**2 / (omega_in * (omega_in + 1j * Gamma_0))
    alpha_      = np.sqrt((omega_in**2 + 1j * omega_in * Gamma_1))
    za          = (alpha_ - omega_1) / (np.sqrt(2) * sigma_1)
    zb          = (alpha_ + omega_1) / (np.sqrt(2) * sigma_1)
    epsilon_out = epsilon_out  +  1j * np.sqrt(pi) * f1 * omega_p**2
                    / (2**1.5 * alpha_ * sigma_1)  .*  ...
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
"""


path1 = 'advanced_testing/SFA030b/SFA030b.mat'
data1 = loadmat(path1)

# load FMM data1
path2 = 'advanced_testing/avg_wire_final/avg_wire_final_Daten_gesamt.mat'
data2 = loadmat(path2)
smat_wire = data2["SMAT_"]

path3 = 'advanced_testing/avg_square_final/avg_square_final_Daten_gesamt.mat'
data3 = loadmat(path3)
smat_square = data3["SMAT_"]

# wavelegth in \mu m
wav_vec = np.linspace(0.47,1.2,512)

# FMM-based data1
n_Futu = np.zeros((1,len(wav_vec)))
n_Sili = np.zeros((1,len(wav_vec)))
print(refractive_ind(wav_vec,12,1))
