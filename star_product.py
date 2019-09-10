import numpy as np


def star_product_analyt(SIN_1,SIN_2):
    """
    Calculate Lifeng Li's starproduct for two S-matrices SIN_1 and SIN_2,
    such that S = S1 * S2. The starproduct between two arbitrary S-matrices
    was precalculated analytically with Mathematica.

    Parameters
    ----------
    SIN_1 : HxLx4x4 numpy array
            H is height_vec_len, the dimension of the height vector
            given to the layer object. (Most of the time equal to 1)
            L is wav_vec_len the number of measured wavelengths
    SIN_2 : HxLx4x4 numpy array
            H is height_vec_len, the dimension of the height vector
            given to the layer object. (Most of the time equal to 1)
            L is wav_vec_len the number of measured wavelengths

    Returns
    -------
    s_out : HxLx4x4 numpy array
    """
    height_vec_len = max(SIN_1.shape[0], SIN_2.shape[0])

    # S-matrix 1
    TF_1 = SIN_1[:,:,0:2,0:2]
    TB_1 = SIN_1[:,:,2:4,2:4]
    RF_1 = SIN_1[:,:,2:4,0:2]
    RB_1 = SIN_1[:,:,0:2,2:4]
    # S-matrix 2
    TF_2 = SIN_2[:,:,0:2,0:2]
    TB_2 = SIN_2[:,:,2:4,2:4]
    RF_2 = SIN_2[:,:,2:4,0:2]
    RB_2 = SIN_2[:,:,0:2,2:4]
    # number of wavelengths
    wav_vec_len = TF_1.shape[1]
    # declare output matrix
    s_out = np.zeros((height_vec_len,wav_vec_len,4,4)).astype(complex)

    left_kernel = np.linalg.inv(np.eye(2) - RB_1 @ RF_2)
    right_kernel = np.linalg.inv(np.eye(2) - RF_2 @ RB_1)

    TF = TF_2 @ left_kernel @ TF_1
    TB = TB_1 @ right_kernel @ TB_2
    RF = RF_1 + TB_1 @ RF_2 @ left_kernel @ TF_1
    RB = RB_2 + TF_2 @ RB_1 @ right_kernel @ TB_2
    # Assemble the resulting s-matrix using the elements from above
    s_out[:,:,0:2,0:2] = TF
    s_out[:,:,2:4,2:4] = TB
    s_out[:,:,2:4,0:2] = RF
    s_out[:,:,0:2,2:4] = RB
    return s_out

def star_product_geometric(SIN_1, SIN_2, order):
    """
    A version of star_product where the [I - a @ b]**-1 term is developed as
    a geometric series to the nth order.

    Parameters
    ----------
    SIN_1: HxLx4x4 numpy array
        H is height_vec_len, the dimension of the height vector
        given to the layer object. (Most of the time equal to 1)
        L is wav_vec_len the number of measured wavelengths
    SIN_2: HxLx4x4 numpy array
        H is height_vec_len, the dimension of the height vector
        given to the layer object. (Most of the time equal to 1)
        L is wav_vec_len the number of measured wavelengths

    order : int

    Returns
    -------
    s_out : HxLx4x4 numpy array
    """
    TF_1 = SIN_1[:,:,0:2,0:2]
    TB_1 = SIN_1[:,:,2:4,2:4]
    RF_1 = SIN_1[:,:,2:4,0:2]
    RB_1 = SIN_1[:,:,0:2,2:4]

    TF_2 = SIN_2[:,:,0:2,0:2]
    TB_2 = SIN_2[:,:,2:4,2:4]
    RF_2 = SIN_2[:,:,2:4,0:2]
    RB_2 = SIN_2[:,:,0:2,2:4]

    height_vec_len = max(SIN_1.shape[0], SIN_2.shape[0])
    wav_vec_len = TF_1.shape[1]
    left_kernel = np.zeros((height_vec_len, wav_vec_len, 2, 2)).astype(complex)
    right_kernel = np.zeros((height_vec_len, wav_vec_len, 2, 2)).astype(complex)

    for n in range(1,order+1):
        left_kernel = left_kernel + np.linalg.matrix_power(RB_1 @ RF_2, n)
        right_kernel = right_kernel + np.linalg.matrix_power(RF_2 @ RB_1, n)

    TF = TF_2 @ TF_1 + TF_2 @ left_kernel @ TF_1
    TB = TB_1 @ TB_2 + TB_1 @ right_kernel @ TB_2
    RF = RF_1 + TB_1 @ RF_2 @ TF_1 + TB_1 @ RF_2 @ left_kernel @ TF_1
    RB = RB_2 + TF_2 @ RB_1 @ TB_2 + TF_2 @ RB_1 @ right_kernel @ TB_2

    s_out = np.zeros((height_vec_len, wav_vec_len, 4, 4)).astype(complex)
    s_out[:,:,0:2,0:2] = TF
    s_out[:,:,2:4,2:4] = TB
    s_out[:,:,2:4,0:2] = RF
    s_out[:,:,0:2,2:4] = RB
    return s_out


def star_product_cascaded(smat_list):
    """
    Iteratively calculates the starproduct (Li, 1996) of N S-matrices, where
    N >= 2. The iteration goes the through the starproduct pair-wise, so
    that: S = ((((((S1 * S2) * S3) * S4) * ... ) * Sn-1) * Sn).

    Parameters
    ----------
    smat_list : list
                A list containing N HxLx4x4 S-matrices

    Returns
    -------
    smat : HxLx4x4 numpy array
    """
    if not type(smat_list) is list:
        raise TypeError("Input has to be a list")
    elif len(smat_list) <= 1:
        raise ValueError("List has to be length 2 or larger")

    smat = smat_list[0]
    for i in range(1, len(smat_list)):
        smat = star_product_analyt(smat, smat_list[i])

    return smat

def star_product_cascaded_geo(smat_list, order):

    """
    A version of star_product_cascaded unsing star_product_geometric.

    Parameters
    ----------
    smat_list : list
                A list containing N HxLx4x4 S-matrices

    order : int

    Returns
    -------
    smat : An L-by-4-by-4 S-matrix.
    """


    if not type(smat_list) is list:
        raise TypeError("Input has to be a list")
    elif len(smat_list) <= 1:
        raise ValueError("List has to be length 2 or larger")

    smat = smat_list[0]

    for i in range(1, len(smat_list)):
        smat = star_product_geometric(smat, smat_list[i], order)

    return smat
