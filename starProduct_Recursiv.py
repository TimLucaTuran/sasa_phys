import numpy as np
from starProductanalyt import starProductanalyt


def starProduct_Recursiv(SMAT_LIST):

    if not type(SMAT_LIST) is list:
        raise TypeError("Input has to be a list")
    elif len(SMAT_LIST) <= 1:
        raise ValueError("List has to be length 2 or larger")

    currentStarMat = SMAT_LIST[0]

    if len(SMAT_LIST) >= 3:
        currentStarMat = starProductanalyt(currentStarMat, starProduct_Recursiv(SMAT_LIST[1:]))
    else:
        currentStarMat = starProductanalyt(currentStarMat, SMAT_LIST[1])
    return currentStarMat

a = np.ones((1,4,4))*4
b = np.ones((1,4,4))*3
c = np.ones((1,4,4))*2
l= [a,b,c]
print(starProduct_Recursiv(l))
