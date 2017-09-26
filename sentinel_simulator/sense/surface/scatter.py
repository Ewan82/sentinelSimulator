"""
Major surface scatter class
"""
class SurfaceScatter(object):
    def __init__(self, eps, ks, theta, kl=None, **kwargs):
        self.eps = eps
        self.ks = ks
        self.theta = theta
        self.kl = kl
        self._check()

    def _check(self):
        assert isinstance(self.eps, complex)

class SurfaceScatterWaterCloudModel(object):
    def __init__(self, mv, theta, C=None, C1=None, C2=None, D=None):
        self.mv = mv
        self.theta = theta
        self.C = C
        self.C1 = C1
        self.C2 = C2
        self.D = D



