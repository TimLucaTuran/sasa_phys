import numpy as np
from star_product import *

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
    def __init__(self, wav_vec, s_mat, cladding, substrate):
        Layer.__init__(self, wav_vec)
        self.s_mat = s_mat
        self.cladding = cladding
        self.substrate = substrate

class NonMetaLayer(Layer):
    """
    Parameters
    ----------
    wav_vec : vector of the measured wavelengths
    height : height in (μm)
    n_vec : one ``
    """
    def __init__(self, wav_vec, height, *n_vec):
        Layer.__init__(self, wav_vec)
        self.height = height
        self.height_len = np.size(self.height)
        self.n_x = n_vec[0]
        #isotropic material
        if len(n_vec) == 1:
            self.n_y = self.n_x
        #anisotropic material
        elif len(n_vec) == 2:
            self.n_y = n_vec[1]
        else:
            raise ValueError("input 1 or 2 refrectiv index vectors")



class Stack:
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
    def __init__(self, layer_list, cladding,
                 clad_height, substrate, subs_height):

        self.layer_list = layer_list
        self.cladding = cladding
        self.clad_height = clad_height
        self.substrate = substrate
        self.subs_height = subs_height
        self.wav_vec = self.layer_list[0].wav_vec
        self.wav_vec_len = len(self.wav_vec)

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

        return np.squeeze(s_mat_list)

    def create_interface(self, l_1, l_2):
        """
        Creates the interface S-Matrix for the transmission between 2 Non-Meta-Layers

        Parameters
        ----------
        l_1 , l_2:  NonMetaLayer or MetaLayer Objects
        """

        #load n_* from the Layers
        if (type(l_1) is NonMetaLayer):
            n1_x = l_1.n_x
            n1_y = l_1.n_y
        else:
            n1_x = l_1.substrate
            n1_y = l_1.substrate

        if(type(l_2) is NonMetaLayer) :
            n2_x = l_2.n_x
            n2_y = l_2.n_y
        else:
            n2_x = l_2.cladding
            n2_y = l_2.cladding

        #transmission and reflection in x and y direction

        T_x = 2*n1_x/(n1_x + n2_x)
        T_y = 2*n1_y/(n1_y + n2_y)
        R_x = (n1_x - n2_x)/(n1_x + n2_x)
        R_y = (n1_y - n2_y)/(n1_y + n2_y)
        s_mat_list = np.zeros((self.wav_vec_len,4,4)).astype(complex)
        s_mat_list[:,0,0] = T_x
        s_mat_list[:,1,1] = T_y
        s_mat_list[:,2,2] = T_x
        s_mat_list[:,3,3] = T_y
        s_mat_list[:,0,2] = R_x
        s_mat_list[:,1,3] = R_y
        s_mat_list[:,2,0] = -1*R_x
        s_mat_list[:,3,1] = -1 * R_y
        """
        S_out =  np.array([[ T_x  , zer    , R_x,    zer],
                         [ zer   , T_y  ,   zer,  R_y],
                         [-1*R_x, zer   , T_x,  zer  ],
                         [ zer    ,-1*R_y, zer  , T_y ]
                        ])
    """
        return s_mat_list


    def build(self):
        """
        Build all the propagation and interface matrices and multiplies them.

        Returns
        -------
        s_mat : Lx4x4 S-matrix describing the whole stack
        """

        #Create Layer-Object for the cladding
        clad_layer = NonMetaLayer(self.wav_vec,
                                  self.clad_height,
                                  self.cladding
                                  )
        #Create Layer-Object for the substrate
        subs_layer = NonMetaLayer(self.wav_vec,
                                  self.subs_height,
                                  self.substrate
                                  )
        #add the substrate layer to the back
        self.layer_list.append(subs_layer)
        #create interface between the cladding and the first layer
        inter = self.create_interface(clad_layer, self.layer_list[0])
        s_mat_list = [inter]
        for i in range(len(self.layer_list) - 1):
            current_layer = self.layer_list[i]
            next_layer = self.layer_list[i+1]

            if type(current_layer) is NonMetaLayer:
                prop = self.create_propagator(current_layer)

            elif type(current_layer) is MetaLayer:
                prop = current_layer.s_mat

            else:
                raise ValueError("Stack has to consist of Mata and \
                                NonMetaLayers")

            inter = self.create_interface(current_layer, next_layer)
            s_mat_list.append(prop)
            print(i)
            print("prop", prop[0,:,:])
            s_mat_list.append(inter)
            print("inter", inter[0,:,:])
        #end building loop


        s_mat_list.append(self.create_propagator(subs_layer))
        #print(s_mat_list)

        return starProduct_Cascaded(s_mat_list)




"""
layer = NonMetaLayer(1/np.arange(1,5), [5.3, 4.3], np.arange(4), 2*np.arange(4))

stack = Stack([layer])
s_mat = stack.create_propagator(stack.layer_list[0])
layer2 = NonMetaLayer()
for i in range(4):
    print(s_mat[1,2,i,i])
"""
