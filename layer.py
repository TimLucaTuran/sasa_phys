import numpy as np

class Layer:
    def __init__(self, wav_vec ):
        self.wav_vec = wav_vec
        self.wav_vec_len = self.wav_vec.size



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
        self.n_material = n_material



class Stack:

    def __init__(self, layer_list):
        """
        layer_list: list of layer-objects
        """
        self.layer_list = layer_list

    def create_propagator(self, non_meta_layer):
        #isotropic material
        if len(non_meta_layer.n_material) == 1:
            n_x = non_meta_layer.n_material[0]
            n_y = non_meta_layer.n_material[0]
        #anisotropic material
        elif len(non_meta_layer.n_material) == 2:
            n_x = non_meta_layer.n_material[0]
            n_y = non_meta_layer.n_material[1]

        #Height is a scalar
        if non_meta_layer.height_len == 1:
            non_meta_layer.height = np.array([non_meta_layer.height])
            """
            prop_x = np.exp(1j * n_x * non_meta_layer.height * 2*np.pi /non_meta_layer.wav_vec)
            prop_y = np.exp(1j * n_y * non_meta_layer.height * 2*np.pi /non_meta_layer.wav_vec)
            s_mat = np.zeros((non_meta_layer.wav_vec_len,4,4)).astype(complex)
            s_mat[:,0,0] = prop_x
            s_mat[:,1,1] = prop_y
            s_mat[:,2,2] = prop_x
            s_mat[:,3,3] = prop_y
            return s_mat
            """

        s_mat_list = np.zeros((non_meta_layer.height_len, non_meta_layer.wav_vec_len,4,4)).astype(complex)
        for i in range(non_meta_layer.height_len):
            prop_x = np.exp(1j * n_x * non_meta_layer.height[i] * 2*np.pi /non_meta_layer.wav_vec)
            prop_y = np.exp(1j * n_y * non_meta_layer.height[i] * 2*np.pi /non_meta_layer.wav_vec)
            s_mat_list[i,:,0,0] = prop_x
            s_mat_list[i,:,1,1] = prop_y
            s_mat_list[i,:,2,2] = prop_x
            s_mat_list[i,:,3,3] = prop_y

        return s_mat_list

    def create_interface(self)



layer = NonMetaLayer(1/np.arange(1,5), [5.3, 4.3], np.arange(4), 2*np.arange(4))
stack = Stack([layer])
s_mat = stack.create_propagator(stack.layer_list[0])
for i in range(4):
    print(s_mat[1,2,i,i])
