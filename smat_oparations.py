import numpy as np

def array_mirror_smat(SMAT):
    """
    INPUT: L x 4 x 4 Array SMAT
    OUTPUT: mirrored SMAT """
    mask = np.array([[1, -1, 1,-1]
                    ,[-1,1,-1,1]
                    ,[1,-1,1,-1]
                    ,[-1,1,-1,1]])
    return SMAT*mask

def array_flip_smat(SMAT):
    """
    Array-Flipping
    INPUT: L x 4 x 4 Array SMAT
    OUTPUT: flipped SMAT
    """
    Sout = np.block([
                    [SMAT[:,2:4,2:4],   SMAT[:,2:4,0:2]],
                    [SMAT[:,0:2,2:4],   SMAT[:,0:2,0:2]]])
    return array_mirror_smat(Sout)

def array_rot_smat(SMAT,ANG):
    """
    Array-Rotation
    INPUT 1: L x 4 x 4 Array SMAT
    INPUT 2: rotationangle in degrees
    OUTPUT: rotated SMAT by angle ANG
    """
    #convert rotation angle in rad
    phi = ANG * np.pi/180
    #number of wavelengths
    numel_wav = SMAT.shape[0]
    # Determine 2x2 rotation matrix to be applied on the matrix blocks
    R = np.array([   [np.cos(phi), -np.sin(phi) ],
            [np.sin(phi), np.cos(phi) ] ])
    #Define right sided rotation operator of size 4x4
    Rot_op = np.block([     [R,                np.zeros((2,2))],
                            [np.zeros((2,2)),  R]])

    #rotate SMAT
    Sout = Rot_op.T @ SMAT @ Rot_op
    return Sout

def array_phase_shift(SMAT,ANG):
    """
    Shifting the Phase of a Matrix
    INPUT 1: L x 4 x 4 Array SMAT
    INPUT 2: shifting angle
    OUTPUT: shifted Matrix
    """
    smat_arg = np.angle(SMAT)
    smat_abs = np.abs(SMAT)
    return smat_abs*np.exp(1j*(smat_arg+ANG))
