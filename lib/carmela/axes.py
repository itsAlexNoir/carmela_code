########################################
##           axes.py
########################################
## This file contains the definition
## of the class axes.
########################################

import numpy as np

class axes():
    class axes:
        """
        A class to represent a 1D axis with evenly spaced points.

        Attributes:
            maxptsx (int): The total number of points along the axis.
            dx (float): The spacing between consecutive points.
            maxx (float): The maximum value of the axis in one direction.
            ax (numpy.ndarray): A 1D array representing the axis values.

        Methods:
            __init__(maxptsx, dx):
                Initializes the axis with the specified number of points and spacing.
        """
    def __init__(self,maxptsx,dx):

        self.maxptsx = maxptsx
        self.dx = dx
        self.maxx = float(self.maxptsx - 1)/2. * self.dx
        self.ax = np.linspace(-self.maxx,self.maxx,self.maxptsx)

class kaxes():
    """
    kaxes class is used to define and initialize axes for a given system, including spatial and energy axes.

    Attributes:
        maxptsx (int): The maximum number of points along the x-axis, inherited from the provided `axes` object.
        dx (float): The spacing between points along the x-axis, calculated based on the `axes` object.
        maxx (float): The maximum value along the x-axis, derived from `maxptsx` and `dx`.
        maxE (float): The maximum energy value for the energy axis. Default is 1.0.
        dE (float): The spacing between points along the energy axis. Default is 1.0.
        maxptsE (int): The number of points along the energy axis, calculated based on `maxE` and `dE`.
        ax (numpy.ndarray): A 1D array representing the spatial axis, ranging from `-maxx` to `maxx` with `maxptsx` points.
        E_ax (numpy.ndarray): A 1D array representing the energy axis, ranging from 0.0 to `maxE` with `maxptsE` points.

    Args:
        axes (object): An object containing the attributes `maxptsx` and `dx`, which are used to initialize the spatial axis.
        maxE (float, optional): The maximum energy value for the energy axis. Default is 1.0.
        dE (float, optional): The spacing between points along the energy axis. Default is 1.0.
    """
    def __init__(self,axes,maxE=1.0,dE=1.0):

        self.maxptsx = axes.maxptsx
        self.dx = 2.0 * np.pi / axes.maxptsx / axes.dx
        self.maxx = float(self.maxptsx - 1) * 0.5 * self.dx
        self.maxE = maxE
        self.dE = dE
        self.maxptsE = int(self.maxE / self.dE) + 1

        self.ax = np.linspace(-self.maxx,self.maxx,self.maxptsx)
        self.E_ax = np.linspace(0.0,self.maxE,self.maxptsE)
        ########################################################

###########################################################
def gaussleg(x1: float, x2: float, degree: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the Gauss-Legendre quadrature points and weights for numerical integration.

    This function calculates the nodes (points) and weights for performing 
    numerical integration using the Gauss-Legendre quadrature method over 
    the interval [x1, x2].

    Parameters:
    -----------
    x1 : float
        The lower bound of the integration interval.
    x2 : float
        The upper bound of the integration interval.
    degree : int
        The number of quadrature points (degree of the polynomial).

    Returns:
    --------
    xx : numpy.ndarray
        The quadrature points mapped to the interval [x1, x2].
    ww : numpy.ndarray
        The corresponding weights for the quadrature points.

    Notes:
    ------
    - The Gauss-Legendre quadrature is a numerical integration method that 
        approximates the integral of a function as a weighted sum of function 
        values at specific points within the integration interval.
    - The `degree` parameter determines the accuracy of the approximation. 
        A higher degree results in more points and greater accuracy.

    Example:
    --------
    >>> x1, x2, degree = 0, 1, 3
    >>> xx, ww = gaussleg(x1, x2, degree)
    >>> print(xx)  # Quadrature points
    >>> print(ww)  # Quadrature weights
    """

    x,w = np.polynomial.legendre.leggauss(degree)

    xl = (x2-x1) * 0.5

    xx = 0.5 * (x + 1) * (x2 - x1) + x1
    ww = xl * w

    return xx, ww
