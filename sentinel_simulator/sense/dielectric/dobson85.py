"""
Dielectric mixing model for soils
after Dobson et al. (1985)
coding after Ulaby (2014), Chapter 4
"""

import numpy as np
from . epsmodel import EpsModel

class Dobson85(EpsModel):
    def __init__(self, **kwargs):
        super(Dobson85, self).__init__(**kwargs)

        self._init_model_parameters()
        self.ew = self._calc_ew()
        self.eps = self._calc_eps()

    def _calc_ew(self, debye=False):
        """
        calculate dielectric permittivity of free water
        using either the Debye model or a more simplistic approach
        """
        if debye:
            return self._debye()
        else:
            # simplistic approach using Eq. 4.69
            return self._simple_ew()

    def _simple_ew(self):
        """ eq. 4.69 """
        f0 = 18.64   # relaxation frequency [GHz]
        hlp = self.f/f0
        e1 = 4.9 + (74.1)/(1.+hlp**2.)
        e2 =(74.1*hlp)/(1.+hlp**2.) + 6.46 * self.sigma/self.f 
        return e1 + 1.j * e2

    def _debye(self):
        assert False


    def _init_model_parameters(self):
        """
        model parameters, eq. 4.68, Ulaby (2014)
        """
        self.alpha = 0.65
        self.beta1 = 1.27-0.519*self.sand - 0.152*self.clay
        self.beta2 = 2.06 - 0.928*self.sand -0.255*self.clay
        self.sigma = -1.645 + 1.939*self.bulk - 2.256*self.sand + 1.594*self.clay

    def _calc_eps(self):
        """
        calculate dielectric permittivity
        Eq. 4.66 (Ulaby et al., 2014)
        """

        e1 = (1.+0.66*self.bulk+self.mv**self.beta1*np.real(self.ew)**self.alpha - self.mv)**(1./self.alpha)
        e2 = np.imag(self.ew)*self.mv**self.beta2
        return e1 + 1.j*e2
        
