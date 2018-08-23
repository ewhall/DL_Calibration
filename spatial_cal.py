import numpy as np

def yz_data(x, y, z, binsize=0.01, retbin=True, retloc=True):

    # get extrema values.
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    
    # make coordinate arrays.
    xi      = np.arange(xmin, xmax+binsize, binsize)
    yi      = np.arange(ymin, ymax+binsize, binsize)
    xi, yi = np.meshgrid(xi,yi)
    
    # make the grid.
    grid           = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape
    if retbin: bins = np.copy(grid)
    
    # create list in same shape as grid to store indices
    if retloc:
        wherebin = np.copy(grid)
        wherebin = wherebin.tolist()

    # fill in the grid.
    for row in range(nrow):
        for col in range(ncol):
            xc = xi[row, col]    # x coordinate.
            yc = yi[row, col]    # y coordinate.
            
            # find the position that xc and yc correspond to.
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)
            ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)
            ind  = np.where(ibin == True)[0]
            
            # fill the bin.
            bin = z[ibin]
            if retloc: wherebin[row][col] = ind
            if retbin: bins[row, col] = bin.size
            if bin.size > 5:
                binval         = np.average(bin)
                grid[row, col] = binval
            else:
                grid[row, col] = np.nan   # fill empty bins with nans.

    # return the grid
    if retbin:
        if retloc:
            return grid, bins, wherebin, xi, yi
        else:
            return grid, bins
    else:
        if retloc:
            return grid, wherebin
        else:
            return grid

def x_data(x, z, binsize=0.01):

    xmin, xmax = x.min(), x.max()

    xi = np.arange(xmin, xmax+binsize, binsize)
    grid = np.zeros(xi.shape, dtype=x.dtype)
    nrow = len(xi)

    for row in range(len(xi)):
        xc = xi[row]

        posx = np.abs(x - xc)
        ibin = np.logical_and(posx < binsize/2., posx < binsize/2.)
        ind = np.where(ibin == True)[0]

        bin = z[ibin]
        if bin.size > 5:
            binval = np.average(bin)
            grid[row] = binval
        else:
            grid[row] = np.NaN

    return grid, xi
