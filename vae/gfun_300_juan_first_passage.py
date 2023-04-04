def gfun_300(x):
    """Performance function for reliability problem 300.

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
    from code_rp.gfun.juan_first_passage.juan_first_passage import juan_first_passage
    from timeit import default_timer as timer

    # expected number of random variables/columns
    nrv_e = 14

    g = float('nan')
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
        start = timer()
        g = juan_first_passage(x)
        end = timer()
        # print(f'Elapsed time: {end - start} sec')

    g_val_sys = g
    g_val_comp = g
    return g_val_sys, g_val_comp, msg
