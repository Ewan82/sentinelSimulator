"""
Definition of scatter types
"""
import numpy as np

class Scatterer(object):
    def __init__(self, **kwargs):
        # NOTE THAT THE arguments are not necessarily the
        # particle scattering cross sections!!!!
        # need to be made in a clearer way !!!!
        self.sigma_s_hh = kwargs.get('sigma_s_hh', None)  # particle scattering cross area
        assert self.sigma_s_hh is not None, 'Particle HH scattering cross section needs to be specified [m**-2]'

        self.sigma_s_vv = kwargs.get('sigma_s_vv', None)  # particle scattering cross area
        assert self.sigma_s_vv is not None, 'Particle VV scattering cross section needs to be specified [m**-2]'

        self.sigma_s_hv = kwargs.get('sigma_s_hv', None)  # particle scattering cross area
        assert self.sigma_s_hv is not None, 'Particle HV scattering cross section needs to be specified [m**-2]'


class ScatIso(Scatterer):
    """
    Isotropic scatterer definition
    see 11.2 in Ulaby (2014)
    """
    def __init__(self, **kwargs):
        super(ScatIso, self).__init__(**kwargs)
        # scattering coefficients need to be the same in isotropic case (eq. 11.19)
        assert self.sigma_s_hh == self.sigma_s_vv
        assert self.sigma_s_hh == self.sigma_s_hv

    def sigma_v_back(self):
        """
        volume backscattering coefficient
        for the isotropic case this corresponds to the
        volume scattering coefficient ks

        not that this is NOT the scattering cross section of a single particle!
        """
        return {'hh' : self.sigma_s_hh, 'vv' : self.sigma_s_vv, 'hv' : self.sigma_s_hv}

    def sigma_v_bist(self):
        # same as volume backscattering coefficient (Eq. 11.19)
        return self.sigma_v_back()



class ScatRayleigh(Scatterer):
    """
    Isotropic scatterer definition
    see 11.2 in Ulaby (2014)
    """
    def __init__(self, **kwargs):
        super(ScatRayleigh, self).__init__(**kwargs)

    def sigma_v_back(self):
        # sigma_s_pp is assumed to correspond to volume extinction coefficient
        return {'hh' : 1.5*self.sigma_s_hh, 'vv' : 1.5*self.sigma_s_vv, 'hv' : np.nan}

    def sigma_v_bist(self):
        # same as sigma_v_back (Eq. 11.22)
        return self.sigma_v_back()
