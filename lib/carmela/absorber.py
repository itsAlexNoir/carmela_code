############################
##       absorber.py      ##
############################
import numpy as np


class absorber():
    """
    A class to model an absorber function that modifies a wavefunction 
    based on a specified split point and absorption edge.

    Attributes:
        xsplit (float): The x-coordinate at which the absorption begins.
        Medge (float): The minimum edge value for the absorption.
        absorbx (numpy.ndarray): The absorption factor array applied to the wavefunction.
        sigma (float): The scaling factor for the absorption function.

    Methods:
        __init__(axes, xsplit, Medge):
            Initializes the absorber with the given axes, split point, and edge value.
        make_split(wavefunc):
            Applies the absorption factor to the given wavefunction and returns the modified wavefunction.
    """
    def __init__(self, axes, xsplit, Medge):
        """
        Initializes the absorber object.

        Parameters:
        -----------
        axes : object
            An object containing the attribute `ax`, which is an array-like structure 
            representing the axis values, and `maxx`, which is the maximum value of the axis.
        xsplit : float
            The split point along the axis beyond which the absorption effect is applied.
        Medge : float
            The minimum edge value for the absorption effect, used to calculate the absorption 
            profile.

        Attributes:
        -----------
        xsplit : float
            Stores the split point along the axis.
        Medge : float
            Stores the minimum edge value for the absorption effect.
        absorbx : numpy.ndarray
            An array representing the absorption profile along the axis. Values are 1 within 
            the `xsplit` range and decrease exponentially beyond it.
        sigma : float
            A scaling factor for the absorption profile, calculated based on `xsplit`, `maxx`, 
            and `Medge`.
        """

        self.xsplit = xsplit
        self.Medge = Medge
        self.absorbx = np.ones(np.shape(axes.ax))

        indx = np.where(abs(axes.ax)>=xsplit)

        self.sigma = (axes.maxx - self.xsplit) / \
                     np.sqrt(-np.log(self.Medge))
        self.absorbx[indx] = \
                            np.exp( - ((abs(axes.ax[indx]) - self.xsplit) / \
                                       self.sigma )**2)

    def make_split(self, wavefunc):
        """
        Splits the given wave function by applying the absorber function.

        This method multiplies the provided wave function by the absorber
        coefficient (`absorbx`) to produce a modified wave function.

        Args:
            wavefunc (array-like): The input wave function to be modified.

        Returns:
            array-like: The modified wave function after applying the absorber.
        """
        wavefunc_split = self.absorbx * wavefunc
        return wavefunc_split

