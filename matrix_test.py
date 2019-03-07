import numpy as np

A = np.ones((2,4,4))
A[0,:,:] = np.arange(16).reshape((1,4,4))
A[1,:,:] = np.arange(16,32).reshape((1,4,4))
print(A)