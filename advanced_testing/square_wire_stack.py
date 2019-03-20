import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.io import loadmat
from stack import MetaLayer, NonMetaLayer, Stack
from scipy.interpolate import PchipInterpolator

path1 = 'advanced_testing/SFA030b/SFA030b.mat'
data1 = loadmat(path1)

# load FMM data1
path2 = 'advanced_testing/avg_wire_final/avg_wire_final_Daten_gesamt.mat'
data2 = loadmat(path2)
smat_wire = data2["SMAT_"]

path3 = 'advanced_testing/avg_square_final/avg_square_final_Daten_gesamt.mat'
data3 = loadmat(path3)
smat_square = data3["SMAT_"]
# print(smat_wire)


# wavelegth in \mu m
wav_vec = np.linspace(0.47,1.2,512)

# FMM-based data
n_Futu = loadmat('Speichern.mat')["n_Futu"].squeeze()
n_Sili = loadmat('Speichern.mat')["n_Sili"].squeeze()


print("1. Futu",n_Futu.shape)

n_vac = np.ones((1,len(wav_vec)))

# SASA
h_spacer      = 450 / 1000
h_cladding    = 585 / 1000

meta_2 = MetaLayer(s_mat = smat_square,
                  cladding = n_Futu,
                  substrate = n_Sili)
meta_1 = MetaLayer(s_mat = smat_wire,
                    cladding= n_Sili,
                    substrate = n_Futu)
spacer_1 = NonMetaLayer(n_Futu,
                        height = h_spacer)
spacer_2 = NonMetaLayer(n_Sili,
                        height = h_cladding)

stack = Stack(layer_list = [meta_1, spacer_1, meta_2, spacer_2],
              wav_vec = wav_vec,
              cladding = n_Sili,
              substrate = n_vac)

s_out = stack.build()
print(s_out[0,:,:])

wl = 1

s_matlab = loadmat('Ergebnis.mat')["SMAT_SASA"][wl-1,:,:]
print(s_out[wl-1,:,:]-s_matlab)



# print(refractive_ind(wav_vec,1,2))
# data_test = h5py.File('advanced_testing/materials/Malitson-o.mat')
# print(data_test["IN_FILE"][:,1])
