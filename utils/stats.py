import numpy as np
import scipy

def make_mesh(n,data_range):
    """Returns a mesh of arbitrary dimensions based on data_range

    Args:
        n (integer): number of points per dimensions
        data_range (list): list of tuples. Length of list is the number of dimensions.
        Each tuple in the list defines the range of each dimension.
    """
    
    vectors = []
    for range in data_range:
        v = np.linspace(range[0],range[1],num=n)
        vectors.append(v)

    vectors = tuple(vectors)

    grid = np.meshgrid(*vectors)

    return grid

def interpolate(x,y,grid=None,n=None):
    """Using scipy interpolate griddata, interpolates x,y over grid values.
    D dimensional data

    Args:
        x (length D tuple of 1-D ndarrays): data points
        y (ndarray): values associated with x data points
        grid (length D tuple of ndarrays): mesh grid to interpolate over
    """

    # if type(grid) is list:
    #     grid = tuple(grid)

    if grid is None:
        grid = []
        if n is None:
            n = 100
        for x_t in x:
            vect = np.linspace(min(x_t),max(x_t),num=n)
    
    values = scipy.interpolate.griddata(x,y,grid)

    return values