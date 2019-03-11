import numpy as np
import scipy.io

data = scipy.io.loadmat('SASA_example_data.mat')

#First Meta-Layer parameters
H1_ind  = 0 # height [50, 55] nm
W1_ind  = 0 # width 75:5:90 nm
L1_ind  = 3 # length 220:5:240 nm

#Second Meta-Layer parameters
H2_ind  = 2 # height [55 65 75] nm
RAD_ind = 2 # corner radius [20 22.5 25] nm
W2_ind  = 3 # width 60:5:80 nm
L2_ind  = 0 # length 220:5:250 nm

SMAT_1 = data["SMAT_1"]
SMAT_1 = np.squeeze(SMAT_1[H1_ind, W1_ind, L1_ind, :,:,:])

SMAT_2 = data["SMAT_2"]
SMAT_2 = np.squeeze(SMAT_2[H2_ind, RAD_ind, W2_ind, L2_ind, :,:,:])

print(SMAT_2)
