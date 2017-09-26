import numpy as np
import matplotlib.pyplot as plt

class Fresnel0(object):
    def __init__(self, e):
        """
        calculate the Nadir Fresnel reflectivity
        e.g. Ulaby (2014), eq. 10.36

        Parameters
        ----------
        e : complex
            complex relative dielectric permitivity
        """
        self.x = self._calc(e)

    def _calc(self, e):
        return np.abs( (1.-np.sqrt(e))/(1.+np.sqrt(e))   )**2.



class Reflectivity(object):
    """
    calculate the reflectivity for H and V polarization
    """
    def __init__(self, eps, theta):
        """
        table 2.5 Ulaby (2014)
        assumes specular surface
        Parameters
        ----------
        eps : complex
            relative dielectric permitivity
        theta : float, ndarray
            incidence angle [rad]
            can be specified
        """
        self.eps = eps
        self.theta = theta

        self._calc_reflection_coefficients()

        self.v = np.abs(self.rho_v)**2.
        self.h = np.abs(self.rho_h)**2.


    def _calc_reflection_coefficients(self):
        """
        calculate reflection coefficients
        Woodhouse, 2006; Eq. 5.54, 5.55
        """
        # OLD
        co = np.cos(self.theta)
        si2 = np.sin(self.theta)**2.
        self.rho_v = (self.eps*co-np.sqrt(self.eps-si2))/(self.eps*co+np.sqrt(self.eps-si2))
        self.rho_h = (co-np.sqrt(self.eps-si2))/(co+np.sqrt(self.eps-si2))

#        srv = self.rho_v
#        srh = self.rho_h
#
#        # FROM PRISM1_FORWARDMODEL-1.m
#        n1 = np.sqrt(1.)
#        n2 = np.sqrt(self.eps)
#        costh2 = np.sqrt(1-(n1*np.sin(self.theta)/2.)**2)
#
#        self.rho_v = -(n2*np.cos(self.theta) - n1*costh2)/(n2*np.cos(self.theta) + n1*costh2)
#        self.rho_h = (n1*np.cos(self.theta) - n2*costh2)/(n1*np.cos(self.theta) + n2*costh2)

#        plt.plot(np.rad2deg(self.theta), self.rho_v-srv, label = 'v_diff')
#        plt.plot(np.rad2deg(self.theta), self.rho_h-srh, label = 'h_diff')
#        plt.legend()
        #doesn't make much difference in results!

    def plot(self):
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.plot(np.rad2deg(self.theta), self.v, color='red', linestyle='-', label='V')
        ax.plot(np.rad2deg(self.theta), self.h, color='blue', linestyle='--', label='H')
        ax.grid()
        ax.legend()













