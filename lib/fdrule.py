#!/usr/bin/env python

import numpy as np

def calculate_fdweights(xi,x,order):

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
