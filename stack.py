import numpy as np

class Layer:
    def __init__(self, wav_vec ):
        self.wav_vec = wav_vec
        self.wav_vec_len = self.wav_vec.size
        self.mirror = False
        self.flip = False
        self.angle = 0

    def set_options(self, mirror=False, flip=False, angle=0):
        self.mirror = mirror
        self.flip = flip
        self.angle = angle

    def flip(self):
        self.flip = True
        return

    def mirror(self):
        self.mirror = True

    def rotate(self, angle):
        self.angle = angle



class MetaLayer(Layer):
    def __init__(self, wav_vec, s_mat, n_embed):
        Layer.__init__(self, wav_vec)
        self.s_mat = s_mat
        self.n_embed = n_embed

class NonMetaLayer(Layer):
    def __init__(self, wav_vec, height, *n_material):
        Layer.__init__(self, wav_vec)
        self.height = height
        self.height_len = np.size(self.height)
        self.n_x = n_material[0]
        #isotropic material
        if len(n_material) == 1:
            self.n_y = self.n_x
        #anisotropic material
        elif len(n_material) == 2:
            self.n_y = n_material[1]
        else:
            raise ValueError("input 1 or 2 refrectiv index vectors")



class Stack:

    def __init__(self, layer_list, cladding, substrate):
        """

        Parameters
        ----------
        layer_list : list of Layer objects
        cladding : float / vector
                   The refrectiv Index of the material on top of the stack
                   if the input is a single float n_i wavelength independent
                   behavior will be assumed.
        substrate : float / vectors
                    The refractiv index of the material below the stack

        """
        self.layer_list = layer_list
        self.cladding = cladding
        self.substrate = substrate

    def create_propagator(self, nml):
        """
        Creates the propergator S-Matrix of a Non-Meta-Layers

        Parameters
        ----------
        nml: NonMetaLayer object
        """

        #Height is a scalar
        if nml.height_len == 1:
            nml.height = np.array([nml.height])

        s_mat_list = np.zeros((nml.height_len, nml.wav_vec_len,4,4)).astype(complex)
        for i in range(nml.height_len):
            prop_x = np.exp(1j * nml.n_x * nml.height[i] * 2*np.pi /nml.wav_vec)
            prop_y = np.exp(1j * nml.n_y * nml.height[i] * 2*np.pi /nml.wav_vec)
            s_mat_list[i,:,0,0] = prop_x
            s_mat_list[i,:,1,1] = prop_y
            s_mat_list[i,:,2,2] = prop_x
            s_mat_list[i,:,3,3] = prop_y

        return s_mat_list

    def create_interface(self, nml_1, nml_2):
        """
        Creates the interface S-Matrix for the transmission between 2 Non-Meta-Layers

        Parameters
        ----------
        nml_1 , nml_2:  NonMetaLayer Objects
        """
        #transmission and reflection in x and y directions
        T_x = 2*nml1.n_x/(nml_1.n_x + nml_2.n_x)
        T_y = 2*nml1.n_y/(nml_1.n_y + nml_2.n_y)
        R_x = (nml_1.n_x - nml_2.n_x)/(nml_1.n_x + nml_2.n_x)
        R_y = (nml_1.n_y - nml_2.n_y)/(nml_1.n_y + nml_2.n_y)
        return np.array([[ T_x  , 0    , R_x,    0],
                         [ 0    , T_y  ,   0,  R_y],
                         [-1*R_x, 0    , T_x,  0  ],
                         [ 0    ,-1*R_y, 0  , T_y ]
                        ])


    def build(self):
        """


        """



layer = NonMetaLayer(1/np.arange(1,5), [5.3, 4.3], np.arange(4), 2*np.arange(4))
stack = Stack([layer])
s_mat = stack.create_propagator(stack.layer_list[0])
for i in range(4):
    print(s_mat[1,2,i,i])
