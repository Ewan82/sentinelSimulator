"""
Class specifying a soil
"""
import numpy as np
from . util import f2lam
from . dielectric import Dobson85

class Soil(object):
    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        eps : complex
            relative permitivity, if this is not given, then mv needs to be given
        s : float
            surface rms height [m]
        mv : float
            volumetric soil moisture [m**3/m**3]; either eps or mv needs to be given
        f : float
            frequency [GHz]
        l : float
            optional: autocorrelation length
        acl : str
            identifier for sphape of autocorrelation fucntion 
            G = Gaussian
            E = Exponential
        mv : float
            volumetric moisture content [m**3/m**3]
        clay : float
            optional fractional clay content
        sand : float
            optional fraction sand content
        """

        self.eps = kwargs.get('eps', None)
        self.mv = kwargs.get('mv', None)
        self.f = kwargs.get('f', None)
        self.s = kwargs.get('s', None)
        self.l = kwargs.get('l', None)
        self.acl = kwargs.get('acl', None)
        self.clay = kwargs.get('clay', None)
        self.sand = kwargs.get('sand', None)
        self.dc_model = kwargs.get('dc_model', 'Dobson85')
        self._check()
        if self.eps is not None:
            self._convert_eps2mv()
        if self.mv is not None:
            self._convert_mv2eps()

        # wavenumber
        self.k = 2.*np.pi / f2lam(self.f)  # note that wavenumber is in meter and NOT in cm!
        
        # roughness parameters
        self.ks = self.s*self.k
        if self.l is not None:
            self.kl = self.k*self.l
        else:
            self.kl = None


    def _convert_mv2eps(self):
        """
        convert mv to eps 
        using dielectric model
        """
        if (self.clay is None) or (self.sand is None):
            self.eps = None
            print('WARNING: Permittivity can not be calculated due to missing soil texture!')
        if self.dc_model == 'Dobson85':
            DC = Dobson85(clay=self.clay, sand=self.sand, mv=self.mv, freq=self.f) 
        else:
            assert False, 'Invalid DC model! ' + self.dc_model
        
        self.eps = DC.eps



    def _convert_eps2mv(self):
        """
        This routine converts soil moisture into
        dielectric properties and vice versa

        future implementations will comprise e.g. the Dobson model
        and others ...


        """
        assert self.eps is not None, 'Currently cnversion not implemented yet; you need to provide the DC directly!'

    def _check(self):
        if self.acl is not None:
            assert self.acl in ['G','E'], 'Invalid form of autocorrelation function specified'
        assert self.s is not None

        if self.eps is None:
            assert self.mv is not None, 'Either EPS or MV need to be given!'
        if self.mv is None:
            assert self.eps is not None, 'Either EPS or MV need to be given!'




