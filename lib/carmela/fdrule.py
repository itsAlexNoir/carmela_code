import numpy as np

def calculate_fdweights(xi,x,order):
    """
    Calculate finite difference weights for numerical differentiation.

    This function computes the weights for finite difference approximations
    of derivatives of a given order at a specific point `xi`, based on a set
    of grid points `x`.

    Parameters:
        xi (float): The point at which the derivative is to be approximated.
        x (array-like): A sequence of grid points used for the finite difference
            approximation. The length of `x` determines the number of points in
            the rule.
        order (int): The order of the derivative to approximate.

    Returns:
        numpy.ndarray: A 2D array of shape (len(x), order + 1) containing the
        finite difference weights. The entry `c[i, k]` corresponds to the weight
        for the k-th derivative at the i-th grid point.

    Notes:
        - The function assumes that the grid points in `x` are distinct.
        - The weights are computed using a recursive algorithm that ensures
            numerical stability.

    Example:
        >>> import numpy as np
        >>> xi = 0.0
        >>> x = np.array([-1.0, 0.0, 1.0])
        >>> order = 2
        >>> weights = calculate_fdweights(xi, x, order)
        >>> print(weights)
    """

    rulepts = len(x)
    c = np.zeros((rulepts,order+1))
    c1 = 1.0
    c4 = x[0] - xi

    c[0,0] = 1.0
    for i in range(1,rulepts):
        mn = min(i,order)
        c2 = 1.0
        c5 = c4
        c4 = x[i] - xi
        for j in range(0,i):
            c3 = x[i] - x[j]
            c2 = c2 * c3
            for k in range(mn,0,-1):
                c[i,k] = c1 * (k * c[i-1,k-1] - c5 * c[i-1,k]) / c2
            c[i,0] = -c1 * c5 * c[i-1,0] / c2
            for k in range(mn,0,-1):
                c[j,k] = (c4 * c[j,k] - k * c[j,k-1]) / c3
            c[j,0] = c4 * c[j,0] / c3
        c1 = c2
    return c
