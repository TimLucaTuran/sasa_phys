import numpy as np
from .star_product import *
from .smat_oparations import *


class Layer:
    """
    Parrent class of Meta- and NonMetaLayer, contains information about
    which symmetry opperations will be applied.
    """

    def __init__(self):
        self.mirror_bool = False
        self.flip_bool = False
        self.angle = 0

    def flip(self):
        self.flip_bool = True
        return

    def mirror(self):
        self.mirror_bool = True

    def rotate(self, angle):
        self.angle = angle


class MetaLayer(Layer):
    """
    Class to describe a Meta-Surface in the Stack.

    Parameters
    ----------
    s_mat : L x 4 x 4 numpy Array
            the Lx4x4 S-Matrix of the Meta-Layer, externally simulated/measured
    cladding : vector
               containing the refraction indices of the cladding.

    substrate : vector
                containing the refraction indices of the substrate.
    """

    def __init__(self, s_mat, cladding, substrate):
        Layer.__init__(self)
        self.s_mat = s_mat
        self.cladding = cladding
        self.substrate = substrate


class NonMetaLayer(Layer):
    """
    Class to describe a homogenous isotropic or anisotropic Layer.

    Parameters
    ----------
    height : height in (μm)
    n_vec : one or two vactors containing the diffraction indeces.
            If only one vector is given homogenous behavior will be assumed.
    """

    def __init__(self, *n_vec, height):
        Layer.__init__(self)
        self.height = height
        self.height_len = np.size(self.height)
        self.n_x = n_vec[0]
        # isotropic material
        if len(n_vec) == 1:
            self.n_y = self.n_x
        # anisotropic material
        elif len(n_vec) == 2:
            self.n_y = n_vec[1]
        else:
            raise ValueError("input 1 or 2 refrectiv index vectors")


class Stack:
    """
    Class to describe the whole Stack, contains information about the layers,
    cladding, substrate and further options.

    Parameters
    ----------
    layer_list : list of Layer objects
    wav_vec : vector
              The target wavelengths where the Meta-Surfaces were simulated/
              measured
    cladding : vector
               The refrectiv indeces of the cladding.
    substrate : vector
                The refractiv indeces of the substrate. The first material to be
                hit by light.

    """

    def __init__(self, layer_list, wav_vec, cladding, substrate):

        self.layer_list = layer_list
        self.cladding = cladding
        self.substrate = substrate
        self.wav_vec = wav_vec
        self.wav_vec_len = len(self.wav_vec)
        self.__geo_bool = False
        self.__geo_order = 5

    def create_propagator(self, layer):
        """
        Creates the propagator S-Matrix

        Parameters
        ----------
        layer : NonMetaLayer or MetaLayer object

        Returns
        -------
        s_mat : H x L x 4 x 4 numpy array
                propagation S-Matrix

        """
        if type(layer) is NonMetaLayer:

            s_mat = np.zeros((layer.height_len, self.wav_vec_len, 4, 4)).astype(complex)
            prop_x = np.exp(2j*np.pi * np.outer(layer.height, layer.n_x/self.wav_vec).squeeze())
            prop_y = np.exp(2j*np.pi * np.outer(layer.height, layer.n_y/self.wav_vec).squeeze())
            s_mat[:, :, 0, 0] = prop_x
            s_mat[:, :, 1, 1] = prop_y
            s_mat[:, :, 2, 2] = prop_x
            s_mat[:, :, 3, 3] = prop_y

        elif type(layer) is MetaLayer:
            s_mat = layer.s_mat.reshape((1, self.wav_vec_len, 4, 4))
        else:
            raise ValueError("Stack has to consist of Mata and \
                            NonMetaLayers")
        # apply symmetry opperations
        if layer.mirror_bool:
            s_mat = mirror_smat(s_mat)
        if layer.flip_bool:
            s_mat = flip_smat(s_mat)
        if layer.angle != 0:
            s_mat = rot_smat(s_mat, layer.angle)

        return s_mat

    def create_interface(self, l_2, l_1):
        """
        Creates the interface S-Matrix for the transmission between two Layers

        Parameters
        ----------
        l_1 :  NonMetaLayer or MetaLayer Objects
        l_2 :  NonMetaLayer or MetaLayer Objects

        Returns
        -------
        s_mat : L x 4 x 4 numpy array
                interface S-Matrix
        """

        # load n_* from the Layers
        if (type(l_1) is NonMetaLayer):
            n1_x = l_1.n_x
            n1_y = l_1.n_y
        else:
            n1_x = l_1.cladding
            n1_y = l_1.cladding

        if(type(l_2) is NonMetaLayer):
            n2_x = l_2.n_x
            n2_y = l_2.n_y
        else:
            n2_x = l_2.substrate
            n2_y = l_2.substrate

        # transmission and reflection in x and y direction

        s_mat_list = np.zeros((self.wav_vec_len, 4, 4)).astype(complex)
        # Transmission
        s_mat_list[:, 0, 0] = 2*n1_x/(n1_x + n2_x)
        s_mat_list[:, 1, 1] = 2*n1_y/(n1_y + n2_y)
        s_mat_list[:, 2, 2] = 2*n2_x/(n1_x + n2_x)
        s_mat_list[:, 3, 3] = 2*n2_y/(n1_y + n2_y)
        # Reflection
        R_x = (n1_x - n2_x)/(n1_x + n2_x)
        R_y = (n1_y - n2_y)/(n1_y + n2_y)
        s_mat_list[:, 0, 2] = R_x
        s_mat_list[:, 1, 3] = R_y
        s_mat_list[:, 2, 0] = -1*R_x
        s_mat_list[:, 3, 1] = -1*R_y
        """
        This Operrator is constructed:
        [T_x  , 0    , R_x,    0],
        [ 0   , T_y  ,   0,  R_y],
        [-1*R_x, 0   , T_x,  0  ],
        [ 0    ,-1*R_y, 0  , T_y ]
        """
        return s_mat_list.reshape((1, self.wav_vec_len, 4, 4))

    def create_interface_rot(self, l_2, l_1):
        """
        Creates the interface S-Matrix for the transmission between
        two Layers in case of rotation, uses create_interface

        Parameters
        ----------
        l_1 :  NonMetaLayer or MetaLayer Objects
        l_2 :  NonMetaLayer or MetaLayer Objects

        Returns
        -------
        s_mat : Lx4x4 S-Matrix
        """
        vacuum_layer = NonMetaLayer(np.ones(self.wav_vec_len), height=None)
        s_mat1 = self.create_interface(vacuum_layer, l_2)
        s_mat2 = self.create_interface(l_1, vacuum_layer)
        s_mat = star_product_analyt(rot_smat(s_mat1, l_2.angle),
                                    rot_smat(s_mat2, l_1.angle))
        return s_mat

    def build(self):
        """
        Builds all the propagation and interface matrices and multiplies them.

        Returns
        -------
        s_mat : Lx4x4 or HxLx4x4 numpy array
                S-matrix describing the behavior of the whole stack. The
                dimension is HxLx4x4 when a height vector was given
        """
        # Create Layer-Objects for the cladding and substrate
        clad_layer = NonMetaLayer(self.cladding, height=None)
        subs_layer = NonMetaLayer(self.substrate, height=None)

        # add the substrate layer to the back
        self.layer_list.append(subs_layer)

        # create interface between the cladding and the first layer
        inter = self.create_interface(clad_layer, self.layer_list[0])
        s_mat_list = [inter]
        for i in range(len(self.layer_list) - 1):

            current_layer = self.layer_list[i]
            next_layer = self.layer_list[i+1]

            prop = self.create_propagator(current_layer)

            # This can be further optimized by a better differentiation between
            # the cases
            if (current_layer.angle != 0) or (next_layer.angle != 0):
                inter = self.create_interface_rot(current_layer, next_layer)
            else:
                inter = self.create_interface(current_layer, next_layer)

            s_mat_list.append(prop)
            s_mat_list.append(inter)
        # end building loop
        if self.__geo_bool:
            s_out = star_product_cascaded_geo(s_mat_list, self.geo_order).squeeze()
        else:
            s_out = star_product_cascaded(s_mat_list).squeeze()

        # remove subs_layer from the layer list
        del self.layer_list[-1]
        return s_out

    def build_geo(self, order):
        """
        A version of build using star_product_cascaded_geo(), change this doc_str

        Returns
        -------
        s_mat : Lx4x4 or HxLx4x4 numpy array
                S-matrix describing the behavior of the whole stack. The
                dimension is HxLx4x4 when a height vector was given
        """
        self.geo_order = order
        self.geo_bool = True
        s_mat = self.build()
        self.geo_bool = False
        return s_mat

    def order(self, order):
        """
        Returns the nth order S-Matrix of the starproduct developt via the
        geometric series.

        Parameters
        ----------
        order : int

        Returns
        -------
        s_out : H x L x 4 x 4 numpy Array
                S-Matrix of the order'th series developt
        """
        self.geo_bool = True
        previous_smat = 0
        if order > 1:
            # calculate previous S-matrix
            self.geo_order = order - 1
            previous_smat = self.build()
        # calculate current S-matrix
        self.geo_order = order
        current_smat = self.build()
        s_out = current_smat - previous_smat
        self.geo_bool = False

        return s_out

    def order_up_to(self, order):
        """
        Builds a list of S-matrices up to the target order.

        Parameters
        ----------
        order : int

        Returns
        -------
        s_list : list of  HxLx4x4 numpy Arrays
        """
        """
        currently cant get this working will use the stupid way,
        maybe the optimisation is unnecassry
        s_list = []
        self.geo_bool = True
        self.order = order
        previous_order = self.build()
        for i in range(order-1, 0, -1):
            self.order = i
            print(i,end=":")
            current_order = self.build()
            print(previous_order[0,0,0,0])
            print(current_order[0,0,0,0])
            s_list.insert(0, previous_order - current_order)
            previous_order = current_order

        self.geo_bool = False
        return s_list
        """
        s_list = []
        self.geo_bool = True
        for i in range(1, order+1):
            s_list.append(self.order(i))
        return s_list
