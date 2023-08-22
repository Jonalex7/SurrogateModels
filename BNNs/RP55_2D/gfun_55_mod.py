def gfun_55_mod(x):
    """Performance function for reliability problem 55.

    Parameters
    ----------
        x : numpy.array of float(s)
            Values of independent variables: columns are the different parameters/random variables (x1, x2,...xn) and rows are different parameter/random variables sets for different calls.

    Returns
    -------
        g_val_sys : numpy.array of float(s)
            Performance function value for the system.
        g_val_comp : numpy.array of float(s)
            Performance function value for each component.
        msg : str
            Accompanying diagnostic message, e.g. warning.
    """
    import numpy as np
    # expected number of random variables/columns
    nrv_e = 2

    g, g1, g2, g3, g4 = float('nan'), float('nan'), float('nan'), float('nan'), float('nan')
    msg = 'Ok'
    x = np.array(x, dtype='f')

    n_dim = len(x.shape)
    if n_dim == 1:
        x = np.array(x)[np.newaxis]
    elif n_dim > 2:
        msg = 'Only available for 1D and 2D arrays.'
        return float('nan'), float('nan'), msg

    nrv_p = x.shape[1]
    if nrv_p != nrv_e:
        msg = f'The number of random variables (x, columns) is expected to be {nrv_e} but {nrv_p} is provided!'
    else:
        g1 = 0.5 + 0.6 * (x[:, 0] - x[:, 1]) ** 4 - (x[:, 0] - x[:, 1]) / np.sqrt(2)
        g2 = 0.5 + 0.6 * (x[:, 0] - x[:, 1]) ** 4 + (x[:, 0] - x[:, 1]) / np.sqrt(2)
        g3 = (x[:, 0] - x[:, 1]) + 5 / np.sqrt(2) - 1.6
        g4 = (x[:, 1] - x[:, 0]) + 5 / np.sqrt(2) - 1.6
        g = np.amin(np.stack((g1, g2, g3, g4)), 0)

    g_val_sys = g
    g_val_comp = np.stack((g1, g2, g3, g4))
    return g_val_sys
