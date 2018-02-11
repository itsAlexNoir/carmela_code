#!/usr/bin/env python

############################
##       absorber.py      ##
############################

import numpy as np


class absorber():
    def __init__(self,axes,xsplit,Medge):

        self.xsplit = xsplit
        self.Medge = Medge
        self.absorbx = np.ones(np.shape(axes.ax))

        indx = np.where(abs(axes.ax)>=xsplit)

        self.sigma = (axes.maxx - self.xsplit) / \
                     np.sqrt(-np.log(self.Medge))
        self.absorbx[indx] = \
                            np.exp( - ((abs(axes.ax[indx]) - self.xsplit) / \
                                       self.sigma )**2)

    def make_split(self,wavefunc):
        wavefunc_split = self.absorbx * wavefunc
        return wavefunc_split

