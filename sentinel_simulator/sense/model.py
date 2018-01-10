"""
Basic class for scattering modelling
"""

import numpy as np
from . surface import Oh92, Dubois95
from . util import f2lam
from . scatterer import ScatIso, ScatRayleigh
from . core import Reflectivity

class Model(object):
    def __init__(self, **kwargs):
        self.theta = kwargs.get('theta', None)

        self._check1()

    def _check1(self):
        assert self.theta is not None, 'ERROR: no incidence angle was specified!'

    def sigma0(self, **kwargs):
        """
        calculate sigma

        Parameters
        ----------
        dB : bool
            return results in decibel
        pol : list
            list with polarizations pq
            whereas p=receive, q=transmit
            p,g can be either H or V
        """

        self.dB = kwargs.get('dB', False)
        self.pol = kwargs.get('pol', [])
        self._check_pol()

        if self.dB:
            assert False, 'Not supported for dictionaries yet!'
            return 10.*np.log10(self._sigma0())
        else:
            return self._sigma0()

    def _sigma0(self, **kwargs):
        assert False, 'routine should be implemented in child class!'

    def _check_pol(self):
        if len(self.pol) == 0:
            assert 'ERROR: polarization needs to be specified'
        for k in self.pol:
            if k not in ['HH','VV','HV','VH']:
                assert False, 'Invalid polarization: ' + k



class SingleScatRT(Model):
    def __init__(self, **kwargs):
        """
        Single scattering model according to Ulaby and Long (2014)
        Eq. 11.17

        Parameters
        ----------
        surface : Surface description
            object describing the surface
        canopy : Canopy description
            object describing the canopy
        models : dict
            dictionary with configuration of scattering models
        """
        super(SingleScatRT, self).__init__(**kwargs)
        self.surface = kwargs.get('surface', None)
        self.canopy = kwargs.get('canopy', None)
        self.models = kwargs.get('models', None)
        self.freq = kwargs.get('freq', None)
        self.coherent = kwargs.get('coherent', True)  # use coherent simulations as default

        self._check()

    def _check(self):
        assert self.surface is not None
        assert self.canopy is not None
        assert self.models is not None
        assert self.freq is not None

        for k in ['surface', 'canopy']:
            assert k in self.models.keys()  # check that all models have been specified

        assert self.freq == self.surface.f, "Different frequencies in model and soil definition"
            # check that frequencies are the same!

    def _sigma0(self):
        """
        basic calculation of Sigma0
        based on Eq. 11.17 in Ulaby and Long (2014)
        """

        # ground backscatter = attenuated surface
        self.G = Ground(self.surface, self.canopy, self.models['surface'], self.models['canopy'], theta=self.theta, freq=self.freq)
        self.s0g = self.G.sigma()  # returns dictionary with different components

        # canopy contribution
        self.s0c = self.G.rt_c.sigma_c()   # returns a dictionary

        # total canopy ground contribution
        self.s0cgt = self.G.sigma_c_g(self.coherent)

        # ground-canopy-ground interaction
        self.s0gcg = self.G.sigma_g_c_g()

        # combine backscatter values
        self.stot = {}
        for k in ['hh', 'vv', 'hv']:
            self.stot.update({k : self._combine(k)})

    def _combine(self, k):
        """
        combine previous calculated backscatter values
        """
        if self.s0g[k] is None:
            return None
        if self.s0c[k] is None:
            return None
        return np.nansum(np.array([self.s0g[k], self.s0c[k], self.s0gcg[k], self.s0cgt[k]]))

class Ground(object):
    """
    calculate the (attenuated) ground contribution
    sigma_pq
    where p is receive and q is transmit polarization
    """
    def __init__(self, S, C, RT_s, RT_c, theta=None, freq=None):
        """
        calculate the attenuated ground contribution
        to the scattering

        Parameters
        ----------
        S : object
            descibing the surface properties
        C : object
            describing the canopy properties
        RT_s : str
            key describing the surface scattering model
        RT_c : str
            key specifying the canopy scattering model
        freq : float
            frequency[GHz]
        """
        self.S = S
        self.C = C
        self.theta = theta
        self._check(RT_s, RT_c)
        self.freq = freq
        assert self.freq is not None, 'Frequency needs to be provided'
        self._set_models(RT_s, RT_c)
        self._calc_rho()

    def _calc_rho(self):
        """
        calculate coherent p-polarized
        reflectivity
        Ref: Eq. 11.11 (Ulaby, 2014)

        Note that the specular reflectivity is corrected by a roughness term
        if ks>0.2

        however, a sensitivity analysis showed that even for ks==0.2
        deviations can be up to 15% for typical incidence angles
        Only in case that ks << 0.1, the correction can be neglected.
        We therefore always use the roughness correction factor!

        TODO: unclear so far how this relates to surface (soil) scattering models
        """

        R = Reflectivity(self.S.eps, self.theta)
        self.rho_v = R.v * np.exp(-4.*np.cos(self.theta)**2.*(self.S.ks**2.))
        self.rho_h = R.h * np.exp(-4.*np.cos(self.theta)**2.*(self.S.ks**2.))

    def _set_models(self, RT_s, RT_c):
        # set surface model
        if RT_s == 'Oh92':
            self.rt_s = Oh92(self.S.eps, self.S.ks, self.theta)
        elif RT_s == 'Dubois95':
            self.rt_s = Dubois95(self.S.eps, self.S.ks, self.theta, lam=f2lam(self.freq))
        else:
            assert False, 'Unknown surface scattering model'


        # set canopy models
        if RT_c == 'turbid_isotropic':  # turbid media (homogenous vegetation)
            self.rt_c = CanopyHomoRT(ke_h=self.C.ke_h, ke_v=self.C.ke_v, ks_h=self.C.ks_h, ks_v=self.C.ks_v, d=self.C.d, theta=self.theta, stype='iso')
        elif RT_c == 'turbid_rayleigh':
            self.rt_c = CanopyHomoRT(ke_h=self.C.ke_h, ke_v=self.C.ke_v, ks_h=self.C.ks_h, ks_v=self.C.ks_v, d=self.C.d, theta=self.theta, stype='rayleigh')
        else:
            assert False, 'Invalid canopy scattering model: ' + RT_c

    def _check(self, RT_s, RT_c):
        valid_surface = ['Oh92', 'Dubois95']
        valid_canopy = ['turbid_rayleigh', 'turbid_isotropic']
        assert RT_s in valid_surface, 'ERROR: invalid surface scattering model was chosen!'
        assert RT_c in valid_canopy, 'ERROR: invalid canopy model: ' + RT_c
        assert self.theta is not None

    def sigma_g_c_g(self):


        s_vv = self.rt_c.sigma_vol_back['vv']*np.cos(self.theta)*self.rho_v*self.rho_v*(self.rt_c.t_v*self.rt_c.t_v-self.rt_c.t_v**4.) / (self.C.ke_v + self.C.ke_v)
        s_hh = self.rt_c.sigma_vol_back['hh']*np.cos(self.theta)*self.rho_h*self.rho_h*(self.rt_c.t_h*self.rt_c.t_h-self.rt_c.t_h**4.) / (self.C.ke_h + self.C.ke_h)
        s_hv = self.rt_c.sigma_vol_back['hv']*np.cos(self.theta)*self.rho_h*self.rho_v*(self.rt_c.t_h*self.rt_c.t_v-self.rt_c.t_h**2.*self.rt_c.t_v**2.) / (self.C.ke_h + self.C.ke_v)

        return {'vv' : s_vv, 'hh' : s_hh, 'hv' : s_hv}


    def sigma_c_g(self, coherent=None):
        """
        calculate canopy ground scattering coefficient
        This is based on Eq. 11.17 (last term) in Ulaby (2014)
        and 11.14 in Ulaby (2014)

        for co-pol, coherent addition can be made as an option

        Parameters
        ----------
        coherent : bool
            do coherent calculation for co-pol calculations
        """
        assert coherent is not None, 'ERROR: please specify explicity if coherent calculations should be made.'
        if coherent:
            n = 2.
        else:
            n = 1.

        s_vv = n  * self.rt_c.sigma_vol_bistatic['vv'] * self.C.d *(self.rho_v + self.rho_v)*self.rt_c.t_v*self.rt_c.t_v
        s_hh = n  * self.rt_c.sigma_vol_bistatic['hh'] * self.C.d *(self.rho_h + self.rho_h)*self.rt_c.t_h*self.rt_c.t_h
        s_hv = 1. * self.rt_c.sigma_vol_bistatic['hv'] * self.C.d *(self.rho_v + self.rho_h)*self.rt_c.t_h*self.rt_c.t_v


        return {'vv' : s_vv, 'hh' : s_hh, 'hv' : s_hv}



    def sigma(self):
        """
        calculate the backscattering coefficient
        Eq. 11.4, p.463 Ulaby (2014)
        """

        # canopy transmisivities
        t_h = self.rt_c.t_h
        t_v = self.rt_c.t_v

        # backscatter
        s_hh = self.rt_s.hh*t_h*t_h
        s_vv = self.rt_s.vv*t_v*t_v

        if self.rt_s.hv is None:
            s_hv = None
        else:
            s_hv = self.rt_s.hv*t_v*t_h


        return {'vv' : s_vv, 'hh' : s_hh, 'hv' : s_hv}


class CanopyHomoRT(object):
    """
    homogeneous canopy RT model
    assumes homogeneous vertical distribution of scatterers

    in that case the Lambert Beer law applies

    NOTE that this model is only for BACKSCATTERING GEOMETRY!
    """
    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        ke_h, ke_v : float
            volume extinction coefficient [Np/m]
        d : float
            height of canopy layer
        theta : float, ndarray
            incidence angle [rad]
        """
        self.ke_h = kwargs.get('ke_h', None)
        self.ke_v = kwargs.get('ke_v', None)
        self.ks_h = kwargs.get('ks_h', None)
        self.ks_v = kwargs.get('ks_v', None)
        self.theta = kwargs.get('theta', None)
        self.d = kwargs.get('d', None)
        #self.Nv = kwargs.get('Nv', 1.)
        self.stype = kwargs.get('stype', None)  # scatterer type

        self._check()

        self.tau_h = self._tau(self.ke_h)
        self.tau_v = self._tau(self.ke_v)
        self.t_h = np.exp(-self.tau_h)
        self.t_v = np.exp(-self.tau_v)

        self._set_scat_type()
        self.sigma_vol_back = self._calc_back_volume()
        self.sigma_vol_bistatic = self._calc_sigma_bistatic()

    def _check(self):
        assert self.stype is not None

        assert self.ke_h is not None
        assert self.ke_v is not None
        assert self.ks_h is not None
        assert self.ks_v is not None

        assert self.ke_h >=0.
        assert self.ke_v >=0.
        assert self.ks_h >=0.
        assert self.ks_v >=0.


        assert self.ks_h <= self.ke_h
        assert self.ks_v <= self.ke_v

    def _set_scat_type(self):
        """ set scatterer type """
        if self.stype == 'iso':
            self.SC = ScatIso(sigma_s_hh=self.ks_h, sigma_s_vv=self.ks_v, sigma_s_hv=self.ks_v)   # note that the cross pol scatt. coeff. is the same as the copol due to isotropic behavior
        elif self.stype == 'rayleigh':
            self.SC = ScatRayleigh(sigma_s_hh = self.ks_h, sigma_s_vv=self.ks_v, sigma_s_hv=self.ks_v)  # eq. 11.22
        elif self.stype == 'cloud':
            assert False  # here implemenatation of 11.5 then
        else:
            assert False, 'Invalid scatterer type specified: ' + self.stype

    def _calc_back_volume(self):
        """
        calculate the volume backscattering coefficient sigma_v
        This is a function of the scatterer type chosen (e.g. isotropic,
        rayleigh, cloud model, ...)
        """
        return self.SC.sigma_v_back()

    def _calc_sigma_bistatic(self):
        """
        calculate volume bistatic scattering coefficient
        of scatterer
        """
        return self.SC.sigma_v_bist()




    def _tau(self, k):
        # assumption: extinction is isotropic
        return k*self.d/np.cos(self.theta)

    def sigma_gcg(self, G_v, G_h):
        """
        calculate ground-canopy-ground interactions
        Eq. 11.16, Ulaby(2014)

        Parameters
        ----------
        G_v : float
            v-polarized coherent Fresnel reflectivity under rough conditions
            see eq. 11.11 for explanations. As this depends on the
            surface model used, these should be provided here explicitely
        G_h : float
            same as above, but for h-polarization.
        """
        return G_v*G_h*(self.t_h*self.t_v-self.t_h**2.*self.t_v**2.)*(self.sigma_vol*np.cos(self.theta))/(self.ke_h+self.ke_v)

    def sigma_c(self):
        """
        calculate canopy volume contribution only
        Eq. 11.10 + 11.16 as seen in 11.17, Ulaby (2014)
        """

        s_hh = (1.-self.t_h*self.t_h)*(self.sigma_vol_back['hh']*np.cos(self.theta))/(self.ke_h+self.ke_h)
        s_vv = (1.-self.t_v*self.t_v)*(self.sigma_vol_back['vv']*np.cos(self.theta))/(self.ke_v+self.ke_v)
        s_hv = (1.-self.t_h*self.t_v)*(self.sigma_vol_back['hv']*np.cos(self.theta))/(self.ke_h+self.ke_v)

        # this seems o.k. here
#        a=self.sigma_vol_back['hh']
#        b=1.5*self.ks_h
#        print a,b,a-b, a/b, self.ks_h

        return {'hh' : s_hh, 'vv' : s_vv, 'hv' : s_hv}


# 502-503



