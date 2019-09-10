import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from stack import MetaLayer, NonMetaLayer, Stack


# Load S-matrices of externally simulated/measured Metasurface
#------------------------------------------------------------------------
# S-matrices are shaped Lx4x4 where L is the number of simulated/measured
# wavelengths

data = loadmat('SASA_example_data.mat')

# First Meta-Layer parameters
H1_ind = 0  # height [50, 55] nm
W1_ind = 0  # width 75:5:90 nm
L1_ind = 3  # length 220:5:240 nm
s_mat_1 = data["SMAT_1"]
s_mat_1 = np.squeeze(s_mat_1[H1_ind, W1_ind, L1_ind, :, :, :])

# Second Meta-Layer parameters
H2_ind = 2  # height [55 65 75] nm
RAD_ind = 2  # corner radius [20 22.5 25] nm
W2_ind = 3  # width 60:5:80 nm
L2_ind = 0  # length 220:5:250 nm
s_mat_2 = data["SMAT_2"]
s_mat_2 = np.squeeze(s_mat_2[H2_ind, RAD_ind, W2_ind, L2_ind, :, :, :])

print(np.matrix(s_mat_1[0,:,:]))
# print(np.matrix(s_mat_2[0,:,:]))

# simulated wavelengths
wavleghts = np.linspace(0.6, 1.62, 64)

# smat_1 and smat_2 are now 64x4x4
# define the refractiv indeces of further materials in the stack at the
# simulated wavelengths
n_SiO2 = data["n_SiO2"]
n_vac = np.ones(wavleghts.size)

# Define the layers in the stack the target looks like this:
# ---------------------------------------------------------------------
# | Cladding | Meta-layer 1 | Spacer | Meta-layer 2 | substrate | Vacuum | <--Light
# ---------------------------------------------------------------------

meta1 = MetaLayer(s_mat = s_mat_1,
                  cladding = n_SiO2,
                  substrate = n_SiO2)

sp_h = 350 / 1000
spacer = NonMetaLayer(n_SiO2, # at this point you can also imput two vectors for
                      height = sp_h) # anisotropic behavior



meta2 = MetaLayer(s_mat = s_mat_2,
                  cladding = n_SiO2,
                  substrate = n_SiO2)

subs_h = 910 / 1000
substrate = NonMetaLayer(n_SiO2,
                         height=subs_h)

# Define the stack
stack = Stack(layer_list = [meta1, spacer, meta2,meta2,meta2, substrate],
              wav_vec = wavleghts,
              cladding = n_SiO2,
              substrate = n_vac) # The "substrate" of the Stack is in this case vacuum

# Optionally apply symmetry opperations to the layers (mirror, flip and rotate)
meta1.rotate(80) #in deg

# calculate the s-matrix describing the whole stack
s_out = stack.build()
my_smat =s_out[0,:,:]
# my_smat_h = np.conj(my_smat).transpose()
# print(my_smat)
# plot the results
intensity = np.abs( s_out[:, 2, 2] )**2 / n_SiO2
plt.plot(wavleghts, np.squeeze(intensity))
plt.show()

# Further fuctionality
# ------------------------------------------------------------------------
#
