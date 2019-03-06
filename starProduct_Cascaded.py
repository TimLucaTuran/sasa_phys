import numpy as np
from starProductanalyt import starProductanalyt

def starProduct_Cascaded(SMAT_LIST):
    """
     Iteratively calculates the starproduct (Li, 1996) of N S-matrices, where
     N >= 2. The iteration goes the through the starproduct pair-wise, so
     that: S = ((((((S1 * S2) * S3) * S4) * ... ) * Sn-1) * Sn).

     Usage: StarMat = starProduct_Cascaded( SMAT_LIST )

     SMAT_LISt: An 1-by-N list containing N L-by-4-by-4 S-matrices,
                whith N >= 2 and where L denotes the number of sampled
                wavelengths. Make sure that L is exactly the same in each
                S-matrix, otherwise starProductanalyt throws an error.

     Output: An L-by-4-by-4 S-matrix.

     Note: On some machines (those with the capability of parallel
           computation) a recursive approach that computes the starproduct
           pair-wise might be more efficient.
    """
    if not type(SMAT_LIST) is list:
        raise TypeError("Input has to be a list")
    elif len(SMAT_LIST) <= 1:
        raise ValueError("List has to be length 2 or larger")

    currentStarMat = SMAT_LIST[0]
    for i in range(1, len(SMAT_LIST)):
        currentStarMat = starProductanalyt(currentStarMat, SMAT_LIST[i])
    return currentStarMat

a = np.ones((1,4,4))*4
b = np.ones((1,4,4))*3
c = np.ones((1,4,4))*2
l= [a,b,c]
print(starProduct_Cascaded(l))
