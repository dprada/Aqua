import numpy as np

def get_norm_vector(vector):

    return np.sqrt(np.dot(vector,vector))

def normalize_vector(vector):

    norm_vector = get_norm_vector(vector)
    if norm_vector > 0.0:
        return vector/norm_vector
    else:
        raise Exception('Norm of vector equal to Zero')
