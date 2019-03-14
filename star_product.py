import numpy as np


def star_product_analyt(SIN_1,SIN_2):
    """
    Calculate Lifeng Li's starproduct for two S-matrices SIN_1 and SIN_2,
    such that S = S1 * S2. The starproduct between two arbitrary S-matrices
    was precalculated analytically with Mathematica.

    Parameters
    ----------
    SIN_1, SIN_2 : HxLx4x4 numpy array
                   H is height_vec_len, the dimension of the height vector
                   given to the layer object. (Most of the time equal to 1)
                   L is wav_vec_len the number of measured wavelengths

    Returns
    -------
    s_out : HxLx4x4 numpy array
    """
    height_vec_len = max(SIN_1.shape[0], SIN_2.shape[0])

    # S-matrix 1
    TF1XX = SIN_1[:,:,0,0]
    TF1XY = SIN_1[:,:,0,1]
    RB1XX = SIN_1[:,:,0,2]
    RB1XY = SIN_1[:,:,0,3]
    TF1YX = SIN_1[:,:,1,0]
    TF1YY = SIN_1[:,:,1,1]
    RB1YX = SIN_1[:,:,1,2]
    RB1YY = SIN_1[:,:,1,3]
    RF1XX = SIN_1[:,:,2,0]
    RF1XY = SIN_1[:,:,2,1]
    TB1XX = SIN_1[:,:,2,2]
    TB1XY = SIN_1[:,:,2,3]
    RF1YX = SIN_1[:,:,3,0]
    RF1YY = SIN_1[:,:,3,1]
    TB1YX = SIN_1[:,:,3,2]
    TB1YY = SIN_1[:,:,3,3]

    # S-matrix 2
    TF2XX = SIN_2[:,:,0,0]
    TF2XY = SIN_2[:,:,0,1]
    RB2XX = SIN_2[:,:,0,2]
    RB2XY = SIN_2[:,:,0,3]
    TF2YX = SIN_2[:,:,1,0]
    TF2YY = SIN_2[:,:,1,1]
    RB2YX = SIN_2[:,:,1,2]
    RB2YY = SIN_2[:,:,1,3]
    RF2XX = SIN_2[:,:,2,0]
    RF2XY = SIN_2[:,:,2,1]
    TB2XX = SIN_2[:,:,2,2]
    TB2XY = SIN_2[:,:,2,3]
    RF2YX = SIN_2[:,:,3,0]
    RF2YY = SIN_2[:,:,3,1]
    TB2YX = SIN_2[:,:,3,2]
    TB2YY = SIN_2[:,:,3,3]
    # number of wavelengths
    wav_vec_len = TF1XX.shape[1]
    # declare output matrix
    SOUT = np.zeros((height_vec_len,wav_vec_len,4,4)).astype(complex)

    # Plain analytic form of the staproduct
    TFXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YX*(RB1XX*RF2XY*TF2XX+
        RB1XY*RF2YY*TF2XX+TF2XY+(-1)*RB1XX*RF2XX*TF2XY+(-1)*RB1XY*
        RF2YX*TF2XY)+TF1XX*(TF2XX+(-1)*RB1YX*RF2XY*TF2XX+(-1)*
        RB1YY*RF2YY*TF2XX+RB1YX*RF2XX*TF2XY+RB1YY*RF2YX*TF2XY))
    # -------------------------------------------------------------------------
    TFXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YY*(RB1XX*RF2XY*TF2XX+
        RB1XY*RF2YY*TF2XX+TF2XY+(-1)*RB1XX*RF2XX*TF2XY+(-1)*RB1XY*
        RF2YX*TF2XY)+TF1XY*(TF2XX+(-1)*RB1YX*RF2XY*TF2XX+(-1)*
        RB1YY*RF2YY*TF2XX+RB1YX*RF2XX*TF2XY+RB1YY*RF2YX*TF2XY))
    # -------------------------------------------------------------------------
    TFYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YX*(RB1XX*RF2XY*TF2YX+
        RB1XY*RF2YY*TF2YX+TF2YY+(-1)*RB1XX*RF2XX*TF2YY+(-1)*RB1XY*
        RF2YX*TF2YY)+TF1XX*(TF2YX+(-1)*RB1YX*RF2XY*TF2YX+(-1)*
        RB1YY*RF2YY*TF2YX+RB1YX*RF2XX*TF2YY+RB1YY*RF2YX*TF2YY))
    # -------------------------------------------------------------------------
    TFYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YY*(RB1XX*RF2XY*TF2YX+
        RB1XY*RF2YY*TF2YX+TF2YY+(-1)*RB1XX*RF2XX*TF2YY+(-1)*RB1XY*
        RF2YX*TF2YY)+TF1XY*(TF2YX+(-1)*RB1YX*RF2XY*TF2YX+(-1)*
        RB1YY*RF2YY*TF2YX+RB1YX*RF2XX*TF2YY+RB1YY*RF2YX*TF2YY))
    # -------------------------------------------------------------------------
    RBXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2XX*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1XY*RB1YX*RF2YY*TB2XX*TF2XX+RB1XY*TB2YX*TF2XX+(-1)
        *RB1XY*RB1YX*RF2XY*TB2YX*TF2XX+RB1YX*TB2XX*TF2XY+(-1)*
        RB1XY*RB1YX*RF2YX*TB2XX*TF2XY+RB1YY*TB2YX*TF2XY+RB1XY*
        RB1YX*RF2XX*TB2YX*TF2XY+RB1XX*(RB1YY*TB2YX*(RF2XY*TF2XX+(
        -1)*RF2XX*TF2XY)+TB2XX*(TF2XX+(-1)*RB1YY*RF2YY*TF2XX+RB1YY*
        RF2YX*TF2XY)))
    # -------------------------------------------------------------------------
    RBXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2XY*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1XY*RB1YX*RF2YY*TB2XY*TF2XX+RB1XY*TB2YY*TF2XX+(-1)
        *RB1XY*RB1YX*RF2XY*TB2YY*TF2XX+RB1YX*TB2XY*TF2XY+(-1)*
        RB1XY*RB1YX*RF2YX*TB2XY*TF2XY+RB1YY*TB2YY*TF2XY+RB1XY*
        RB1YX*RF2XX*TB2YY*TF2XY+RB1XX*(RB1YY*TB2YY*(RF2XY*TF2XX+(
        -1)*RF2XX*TF2XY)+TB2XY*(TF2XX+(-1)*RB1YY*RF2YY*TF2XX+RB1YY*
        RF2YX*TF2XY)))
    # -------------------------------------------------------------------------
    RBYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2YX*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1XY*RB1YX*RF2YY*TB2XX*TF2YX+RB1XY*TB2YX*TF2YX+(-1)
        *RB1XY*RB1YX*RF2XY*TB2YX*TF2YX+RB1YX*TB2XX*TF2YY+(-1)*
        RB1XY*RB1YX*RF2YX*TB2XX*TF2YY+RB1YY*TB2YX*TF2YY+RB1XY*
        RB1YX*RF2XX*TB2YX*TF2YY+RB1XX*(RB1YY*TB2YX*(RF2XY*TF2YX+(
        -1)*RF2XX*TF2YY)+TB2XX*(TF2YX+(-1)*RB1YY*RF2YY*TF2YX+RB1YY*
        RF2YX*TF2YY)))
    # -------------------------------------------------------------------------
    RBYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2YY*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1XY*RB1YX*RF2YY*TB2XY*TF2YX+RB1XY*TB2YY*TF2YX+(-1)
        *RB1XY*RB1YX*RF2XY*TB2YY*TF2YX+RB1YX*TB2XY*TF2YY+(-1)*
        RB1XY*RB1YX*RF2YX*TB2XY*TF2YY+RB1YY*TB2YY*TF2YY+RB1XY*
        RB1YX*RF2XX*TB2YY*TF2YY+RB1XX*(RB1YY*TB2YY*(RF2XY*TF2YX+(
        -1)*RF2XX*TF2YY)+TB2XY*(TF2YX+(-1)*RB1YY*RF2YY*TF2YX+RB1YY*
        RF2YX*TF2YY)))
    # -------------------------------------------------------------------------
    RFXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1XX*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1YY*RF2XY*RF2YX*TB1XX*TF1XX+RF2YX*TB1XY*TF1XX+(-1)
        *RB1YX*RF2XY*RF2YX*TB1XY*TF1XX+RF2XY*TB1XX*TF1YX+(-1)*
        RB1XY*RF2XY*RF2YX*TB1XX*TF1YX+RB1XX*RF2XY*RF2YX*TB1XY*
        TF1YX+RF2YY*TB1XY*TF1YX+RF2XX*(RF2YY*TB1XY*(RB1YX*TF1XX+(-1)
        *RB1XX*TF1YX)+TB1XX*(TF1XX+(-1)*RB1YY*RF2YY*TF1XX+RB1XY*
        RF2YY*TF1YX)))
    # -------------------------------------------------------------------------
    RFXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1XY*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1YY*RF2XY*RF2YX*TB1XX*TF1XY+RF2YX*TB1XY*TF1XY+(-1)
        *RB1YX*RF2XY*RF2YX*TB1XY*TF1XY+RF2XY*TB1XX*TF1YY+(-1)*
        RB1XY*RF2XY*RF2YX*TB1XX*TF1YY+RB1XX*RF2XY*RF2YX*TB1XY*
        TF1YY+RF2YY*TB1XY*TF1YY+RF2XX*(RF2YY*TB1XY*(RB1YX*TF1XY+(-1)
        *RB1XX*TF1YY)+TB1XX*(TF1XY+(-1)*RB1YY*RF2YY*TF1XY+RB1XY*
        RF2YY*TF1YY)))
    # -------------------------------------------------------------------------
    RFYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1YX*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1YY*RF2XY*RF2YX*TB1YX*TF1XX+RF2YX*TB1YY*TF1XX+(-1)
        *RB1YX*RF2XY*RF2YX*TB1YY*TF1XX+RF2XY*TB1YX*TF1YX+(-1)*
        RB1XY*RF2XY*RF2YX*TB1YX*TF1YX+RB1XX*RF2XY*RF2YX*TB1YY*
        TF1YX+RF2YY*TB1YY*TF1YX+RF2XX*(RF2YY*TB1YY*(RB1YX*TF1XX+(-1)
        *RB1XX*TF1YX)+TB1YX*(TF1XX+(-1)*RB1YY*RF2YY*TF1XX+RB1XY*
        RF2YY*TF1YX)))
    # -------------------------------------------------------------------------
    RFYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1YY*(((-1)+RB1YX*RF2XY)*
        ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
        RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
        RF2YY)+RB1YY*RF2XY*RF2YX*TB1YX*TF1XY+RF2YX*TB1YY*TF1XY+(-1)
        *RB1YX*RF2XY*RF2YX*TB1YY*TF1XY+RF2XY*TB1YX*TF1YY+(-1)*
        RB1XY*RF2XY*RF2YX*TB1YX*TF1YY+RB1XX*RF2XY*RF2YX*TB1YY*
        TF1YY+RF2YY*TB1YY*TF1YY+RF2XX*(RF2YY*TB1YY*(RB1YX*TF1XY+(-1)
        *RB1XX*TF1YY)+TB1YX*(TF1XY+(-1)*RB1YY*RF2YY*TF1XY+RB1XY*
        RF2YY*TF1YY)))
    # -------------------------------------------------------------------------
    TBXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1XY*(RB1XX*RF2YX*TB2XX+
        RB1YX*RF2YY*TB2XX+TB2YX+(-1)*RB1XX*RF2XX*TB2YX+(-1)*RB1YX*
        RF2XY*TB2YX)+TB1XX*(TB2XX+(-1)*RB1XY*RF2YX*TB2XX+(-1)*
        RB1YY*RF2YY*TB2XX+RB1XY*RF2XX*TB2YX+RB1YY*RF2XY*TB2YX))
    # -------------------------------------------------------------------------
    TBXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1XY*(RB1XX*RF2YX*TB2XY+
        RB1YX*RF2YY*TB2XY+TB2YY+(-1)*RB1XX*RF2XX*TB2YY+(-1)*RB1YX*
        RF2XY*TB2YY)+TB1XX*(TB2XY+(-1)*RB1XY*RF2YX*TB2XY+(-1)*
        RB1YY*RF2YY*TB2XY+RB1XY*RF2XX*TB2YY+RB1YY*RF2XY*TB2YY))
    # -------------------------------------------------------------------------
    TBYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1YY*(RB1XX*RF2YX*TB2XX+
        RB1YX*RF2YY*TB2XX+TB2YX+(-1)*RB1XX*RF2XX*TB2YX+(-1)*RB1YX*
        RF2XY*TB2YX)+TB1YX*(TB2XX+(-1)*RB1XY*RF2YX*TB2XX+(-1)*
        RB1YY*RF2YY*TB2XX+RB1XY*RF2XX*TB2YX+RB1YY*RF2XY*TB2YX))
    # -------------------------------------------------------------------------
    TBYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
        RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
        RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1YY*(RB1XX*RF2YX*TB2XY+
        RB1YX*RF2YY*TB2XY+TB2YY+(-1)*RB1XX*RF2XX*TB2YY+(-1)*RB1YX*
        RF2XY*TB2YY)+TB1YX*(TB2XY+(-1)*RB1XY*RF2YX*TB2XY+(-1)*
        RB1YY*RF2YY*TB2XY+RB1XY*RF2XX*TB2YY+RB1YY*RF2XY*TB2YY))


    # Assemble the resulting s-matrix using the elements from above

    SOUT[:,:,0,0] = TFXX
    SOUT[:,:,0,1] = TFXY
    SOUT[:,:,0,2] = RBXX
    SOUT[:,:,0,3] = RBXY
    SOUT[:,:,1,0] = TFYX
    SOUT[:,:,1,1] = TFYY
    SOUT[:,:,1,2] = RBYX
    SOUT[:,:,1,3] = RBYY
    SOUT[:,:,2,0] = RFXX
    SOUT[:,:,2,1] = RFXY
    SOUT[:,:,2,2] = TBXX
    SOUT[:,:,2,3] = TBXY
    SOUT[:,:,3,0] = RFYX
    SOUT[:,:,3,1] = RFYY
    SOUT[:,:,3,2] = TBYX
    SOUT[:,:,3,3] = TBYY

    return SOUT

def star_product_geometric(SIN_1, SIN_2, order):
    """
    A version star_product where the [I - a @ b]**-1 term in the star product
    is developt as a geometric series to the nth order.

    Parameters
    ----------
    SIN_1, SIN_2 : HxLx4x4 numpy array
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
