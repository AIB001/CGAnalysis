import numpy as np

def apply_pbc(distance_vector, box_size):
    """
    Apply periodic boundary conditions to a distance vector.

    Parameters:
    - distance_vector: Numpy array of shape (3,) representing the distance.
    - box_size: Box size along each dimension.

    Returns:
    - Distance vector with PBC applied.
    """
    return distance_vector - box_size * np.round(distance_vector / box_size)
