# Thinly wrapped numpy
import autograd.numpy as np
# Basically everything you need
from autograd import grad
# Define a function like normal with Python and Numpy
def tanh(x):
    y = np.exp(-x)
    return (1.0 - y) / (1.0 + y)
# Create a function to compute the gradient
grad_tanh = grad(np.sin)
# Evaluate the gradient at x = 1.0
print(grad_tanh(0.0))
