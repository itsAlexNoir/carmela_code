#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
##           axes.py
########################################
## This file contains the definition
## of the class axes.
########################################

import numpy as np

class axes():
    def __init__(self,maxptsx,dx):

        self.maxptsx = maxptsx
        self.dx = dx
        self.maxx = float(self.maxptsx - 1)/2. * self.dx
        ####################################################
        self.ax = np.linspace(-self.maxx,self.maxx,self.maxptsx)
        ######################################################
        ######################################################

class kaxes():
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
def gaussleg(x1,x2,degree):

    x,w = np.polynomial.legendre.leggauss(degree)

    xl = (x2-x1) * 0.5
    xm = (x2+x1) * 0.5

    xx = 0.5 * (x + 1) * (x2 - x1) + x1
    ww = xl * w

    return xx, ww
