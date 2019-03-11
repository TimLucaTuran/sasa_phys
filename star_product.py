import numpy as np


def starProductanalyt(SIN_1,SIN_2):
    SIN_1 = np.squeeze(SIN_1)
    SIN_2 = np.squeeze(SIN_2)
    # S-matrix 1
    TF1XX = np.squeeze(SIN_1[:,0,0])
    TF1XY = np.squeeze(SIN_1[:,0,1])
    RB1XX = np.squeeze(SIN_1[:,0,2])
    RB1XY = np.squeeze(SIN_1[:,0,3])
    TF1YX = np.squeeze(SIN_1[:,1,0])
    TF1YY = np.squeeze(SIN_1[:,1,1])
    RB1YX = np.squeeze(SIN_1[:,1,2])
    RB1YY = np.squeeze(SIN_1[:,1,3])
    RF1XX = np.squeeze(SIN_1[:,2,0])
    RF1XY = np.squeeze(SIN_1[:,2,1])
    TB1XX = np.squeeze(SIN_1[:,2,2])
    TB1XY = np.squeeze(SIN_1[:,2,3])
    RF1YX = np.squeeze(SIN_1[:,3,0])
    RF1YY = np.squeeze(SIN_1[:,3,1])
    TB1YX = np.squeeze(SIN_1[:,3,2])
    TB1YY = np.squeeze(SIN_1[:,3,3])

    # S-matrix 2
    TF2XX = np.squeeze(SIN_2[:,0,0])
    TF2XY = np.squeeze(SIN_2[:,0,1])
    RB2XX = np.squeeze(SIN_2[:,0,2])
    RB2XY = np.squeeze(SIN_2[:,0,3])
    TF2YX = np.squeeze(SIN_2[:,1,0])
    TF2YY = np.squeeze(SIN_2[:,1,1])
    RB2YX = np.squeeze(SIN_2[:,1,2])
    RB2YY = np.squeeze(SIN_2[:,1,3])
    RF2XX = np.squeeze(SIN_2[:,2,0])
    RF2XY = np.squeeze(SIN_2[:,2,1])
    TB2XX = np.squeeze(SIN_2[:,2,2])
    TB2XY = np.squeeze(SIN_2[:,2,3])
    RF2YX = np.squeeze(SIN_2[:,3,0])
    RF2YY = np.squeeze(SIN_2[:,3,1])
    TB2YX = np.squeeze(SIN_2[:,3,2])
    TB2YY = np.squeeze(SIN_2[:,3,3])
    # number of wavelengths
    numwavel = TF1XX.size


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

    SOUT = np.zeros((numwavel,4,4)).astype(complex)

    SOUT[:,0,0] = TFXX
    SOUT[:,0,1] = TFXY
    SOUT[:,0,2] = RBXX
    SOUT[:,0,3] = RBXY
    SOUT[:,1,0] = TFYX
    SOUT[:,1,1] = TFYY
    SOUT[:,1,2] = RBYX
    SOUT[:,1,3] = RBYY
    SOUT[:,2,0] = RFXX
    SOUT[:,2,1] = RFXY
    SOUT[:,2,2] = TBXX
    SOUT[:,2,3] = TBXY
    SOUT[:,3,0] = RFYX
    SOUT[:,3,1] = RFYY
    SOUT[:,3,2] = TBYX
    SOUT[:,3,3] = TBYY

    return SOUT



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

    StarMat = SMAT_LIST[0]
    print(SMAT_LIST[0].shape)
    for i in range(1, len(SMAT_LIST)):
        currentStarMat = starProductanalyt(StarMat, SMAT_LIST[i])
        
    return StarMat



    """
def starProduct_Recursiv(SMAT_LIST):
    A recursiv variant of starProduct_Cascaded

    if not type(SMAT_LIST) is list:
        raise TypeError("Input has to be a list")
    elif len(SMAT_LIST) <= 1:
        raise ValueError("List has to be length 2 or larger")

    StarMat = SMAT_LIST[0]

    if len(SMAT_LIST) >= 3:
        StarMat = starProductanalyt(StarMat, starProduct_Recursiv(SMAT_LIST[1:]))
    else:
        StarMat = starProductanalyt(StarMat, SMAT_LIST[1])
    return StarMat
    """
