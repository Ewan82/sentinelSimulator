#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename: prevot1993.py
# Author: "Thomas Weiß"
# Date:   2017-07-26 13:25:10
# Last Modified by:   "Thomas Weiß"
# Last Modified time: 2017-07-26 13:55:33

"""
This module implements the Prevot et al. (1993) empirical vegetation backscattering model component

References
----------
Prevot et al. (1992): Estimation the characteristics of vegetation canopy with airborne radar measurements
"""

import numpy as np
import matplotlib.pyplot as plt

from . scatter import VegetationScatterWaterCloudModel

class Prevot93_vegetation(VegetationScatterWaterCloudModel):
    def __init__(self, V, A, B, theta):
        """
        Parameters
        ----------
        V : float
            water content of the canopy [kg/m**2]
        A : float
            correspond to the albedo of the vegetation
        B : float
            attenuation factor
        theta : float
            incidence angle [rad]
        """
        super(Prevot93_vegetation, self).__init__(theta, V=V, A=A, B=B)

        # calculate backscatter
        self.tau = np.exp(-2 * self.B * self.V / np.cos(np.radians(self.theta)))
        self.vegetation = self.A * np.cos(np.radians(self.theta)) * (1 - self.tau)

    def plot(self):
        f = plt.figure()
        ax = f.add_subplot(111)
        t = np.rad2deg(self.theta)
        ax.plot(t, 10.*np.log10(self.vegetation), color='blue', label='vegetation component')
        ax.grid()
        # ax.set_ylim(-25.,0.)
        # ax.set_xlim(0.,70.)
        ax.legend()
        ax.set_xlabel('incidence angle [deg]')
        ax.set_ylabel('backscatter [dB]')
        plt.show()
        pdb.set_trace()




