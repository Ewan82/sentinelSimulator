"""
Module for some utilty functions
"""

c0=299792458.  # speed of light [m/s]

def f2lam(f):
    """
    given the frequency in GHz,
    return the wavelength [m]
    """
    return c0/(f*1.E9)
