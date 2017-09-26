#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename: prevot1993.py
# Author: "Thomas Weiß"
# Date:   2017-07-25 16:25:48
# Last Modified by:   "Thomas Weiß"
# Last Modified time: 2017-07-26 13:55:20

"""
This module implements the Prevot et al. (1993)
empirical surface backscattering model component

References
----------
Prevot et al. (1992): Estimation the characteristics of vegetation canopy with airborne radar measurements
"""

import numpy as np
import matplotlib.pyplot as plt

from . scatter import SurfaceScatterWaterCloudModel
import pdb

class Prevot93_surface(SurfaceScatterWaterCloudModel):
    def __init__(self, mv, theta, C1, C2, D):
        """
        Parameters
        ----------
        mv : float
            volumetric soil moisture [m**3/m**3]
        C1 : float
            calibration constant [dB]
        C2 : float
            angular sensitivity of the soil signal [dB/°] (related to roughness)
        theta : float
            incidence angle [rad]
        D : float
            sensitivity of the signal to soil moisture
        """
        super(Prevot93_surface, self).__init__(mv, theta, C1=C1, C2=C2, D=D)
        pdb.set_trace()

        # calculate backscatter
        self.C = self.C1 - self.C2 * self.theta
        self.surface = self.C + self.D * self.mv

    def plot(self):
        f = plt.figure()
        ax = f.add_subplot(111)
        t = np.rad2deg(self.theta)
        ax.plot(t, 10.*np.log10(self.surface), color='blue', label='surface component')
        ax.grid()
        # ax.set_ylim(-25.,0.)
        # ax.set_xlim(0.,70.)
        ax.legend()
        ax.set_xlabel('incidence angle [deg]')
        ax.set_ylabel('backscatter [dB]')
        plt.show()
        pdb.set_trace()







