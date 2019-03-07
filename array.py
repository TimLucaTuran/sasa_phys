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
    """
    #rotation angle in rad
    phi = ANG * np.pi/180
    #number of wavelengths
    numel_wav = SMAT.shape[0]
    # Determine 2x2 rotation matrix to be applied on the matrix blocks
    R = np.array([   [np.cos(phi), -np.sin(phi) ],
            [np.sin(phi), np.cos(phi) ] ])
    #Define right sided rotation operator of size 4x4
    Rot_op = np.block([     [R,                np.zeros((2,2))],
                            [np.zeros((2,2)),  R]])
    #Define left sided rotation operator of size 4x4
    Rot_op_t = Rot_op.T
    #rotate SMAT 
    Sout = Rot_op_t @ SMAT @ Rot_op
    return Sout
        
A = np.ones((2,4,4))
A[0,:,:] = np.arange(16).reshape((1,4,4))
A[1,:,:] = np.arange(16,32).reshape((1,4,4))
print(array_rot_smat(A,90))