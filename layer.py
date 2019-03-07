import numpy as np

class Layer:
    def __init__(self, wav_vec ):
        self.wav_vec = wav_vec



class MetaLayer(Layer):
    def __init__(self, wav_vec, s_mat, n_embed):
        Layer.__init__(self, wav_vec)
        self.s_mat = s_mat
        self.n_embed = n_embed

class NonMetaLayer(Layer):
    def __init__(self, wav_vec, height, n_material):
        Layer.__init__(self, wav_vec)
        self.height = height
        self.n_material = n_material
