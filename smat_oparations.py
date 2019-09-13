import numpy as np

def mirror_smat(s_mat):
    """
    Mirror a given S-Matrix

    Parameters
    ----------
    s_mat: L x 4 x 4 numpy Array
        S-Matrix

    Returns
    -------
    s_out: L x 4 x 4 numpy Array
        mirrored S-Matrix
    """
    mask = np.array([[1, -1, 1,-1]
                    ,[-1,1,-1,1]
                    ,[1,-1,1,-1]
                    ,[-1,1,-1,1]])
    s_out = s_mat * mask
    return s_out

def flip_smat(s_mat):
    """
    Flip a given S-Matrix

    Parameters
    ----------
    s_mat: L x 4 x 4 numpy Array
        S-Matrix

    Returns
    -------
    s_out: L x 4 x 4 numpy Array
        flipped S-Matrix
    """
    s_out = np.block([
                    [SMAT[:,2:4,2:4],   SMAT[:,2:4,0:2]],
                    [SMAT[:,0:2,2:4],   SMAT[:,0:2,0:2]]])
    s_out = array_mirror_smat(Sout)
    return s_out

def rot_smat(s_mat,ang):
    """
    Rotate a given S-Matrix by a given angle

    Parameters
    ----------
    s_mat : Lx4x4 Array
    ang   : float
            rotaion angle in rad

    Returns
    -------
    s_out: Lx4x4 Array
           rotated S-Matrix
    """
    #convert rotation angle in rad
    phi = ang * np.pi/180
    #number of wavelengths
    numel_wav = s_mat.shape[0]
    # Determine 2x2 rotation matrix to be applied on the matrix blocks
    R = np.array([[np.cos(phi), -np.sin(phi) ],
                  [np.sin(phi),  np.cos(phi) ] ])
    #Define right sided rotation operator of size 4x4
    Rot_op = np.block([     [R,                np.zeros((2,2))],
                            [np.zeros((2,2)),  R]])

    #rotate SMAT
    s_out = Rot_op.T @ s_mat @ Rot_op
    return s_out

def phase_shift(smat,ang):
    """
    Shifting the phase of a given S-Matrix by a given angle

    Parameters
    ----------
    s_mat: L x 4 x 4 numpy Array
        S-Matrix
    ang: float
        rotaion angle in rad

    Returns
    -------
    s_out: L x 4 x 4 numpy Array
        shifted S-Matrix
    """
    smat_arg = np.angle(smat)
    smat_abs = np.abs(smat)
    s_out = smat_abs*np.exp(1j*(smat_arg+ang))
    return s_out
