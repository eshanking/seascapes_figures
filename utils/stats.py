import numpy as np
import scipy.interpolate

def interpolate(x,y,grid=None,n=None,**kwargs):
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
            grid.append(vect)
    
    grid = tuple(grid)
    grid = np.meshgrid(*grid)
    grid = tuple(grid)

    values = scipy.interpolate.griddata(x,y,grid,**kwargs)

    return values