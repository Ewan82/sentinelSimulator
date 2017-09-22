#!/usr/bin/env python

import re
from copy import copy
import numpy as np


class UnknownFileType(Exception):
    """Exception class for unknown filetypes
    """
    pass


class spectra(object):
    """Spectra class for sentinel simulator.
    """
    def __init__(self, fname=None, ftype="SVC", wavlCol=0, reflCol=1, hdrLines=1):

        self.refl = np.array([])
        self.wavl = np.array([])
        self.ftype = ftype

        if fname != None:
            self.loadSpectra(fname, wavlCol, reflCol, hdrLines)

    def loadSpectra(self, fname, wavlCol=0, reflCol=1, hdrLines=1):
        """Load in the spectra from a given file using a
        method appropriate to the type of file.

        .. note:: Current supported formats are:

         SVC - SCV .sig ascii file

         CSV - standard ascii comma seperated values

        :param fname: Valid filename containing spectra.
        :type fname: str
        :param wavlCol: Column containing wavelengths.
        :type wavlCol: int
        :param reflCol: Column containing reflectance data.
        :type reflCol: int
        :param hdrLines: Number of lines to skip at start of file.
        :type hdrLines: int
        """

        with open(fname) as f:
            if self.ftype == "SVC":
                self.loadSVCSig(f)
            elif self.ftype == "CSV":
                self.loadCSV(f, wavlCol, reflCol, hdrLines)
            else:
                raise UnknownFileType(self.ftype)
            f.close()

    def loadCSV(self, f, wavlCol=0, reflCol=1, hdrLines=1):
        """Read in data from a standard CSV file object.

        :param f: File object.
        :type f: file
        :param wavlCol: Column containing wavelengths.
        :type wavlCol: int
        :param reflCol: Column containing reflectance data.
        :type reflCol: int
        :param hdrLines: Number of lines to skip at start of file.
        :type hdrLines: int
        """

        tmp = np.loadtxt(f, delimiter=",", skiprows=hdrLines, usecols=(wavlCol, reflCol))
        self.wavl = tmp[:, 0]
        self.refl = tmp[:, 1]

    def loadSVCSig(self, f):
        """Read in data from an SVC .sig ascii file.

        :param f: File object.
        :type f: file
        """

        getData = False

        for line in f:
            if getData:
                reflTmp = np.append(self.refl, float(line.split()[3]) / 100.)
                wavlTmp = np.append(self.wavl, float(line.split()[0]))
                self.refl = copy(reflTmp)
                self.wavl = copy(wavlTmp)
            if re.match('data=', line):
                getData = True

    def interpolate(self, resltn=0.1):
        """Interpolate spectra to the given resolution.
        Overwites exisiting data.

        :param resltn: resolution of the interpolation.
        :type resltn: float
        """

        # find the starting and ending wavelengths
        begWavl = np.ceil(self.wavl[0] / resltn) * resltn
        endWavl = np.floor(self.wavl[-1] / resltn) * resltn

        # print self.wavl[0], begWavl
        # print self.wavl[-1], endWavl

        # generate new wavelength and relfectance arrays
        wavlTmp = np.arange(begWavl, endWavl + resltn, resltn)
        reflTmp = np.zeros(np.shape(wavlTmp))

        # perfrom a linear interpolation:
        m = 0
        for (n, wavl) in enumerate(wavlTmp):

            while self.wavl[m] < wavl and wavl != wavlTmp[-1]:
                m += 1

            if self.wavl[m] == wavl:
                reflTmp[n] = self.refl[m]
            else:
                w1 = self.wavl[m - 1]
                w2 = self.wavl[m]
                r1 = self.refl[m - 1]
                r2 = self.refl[m]
                f = (w2 - wavl) / (w2 - w1)
                reflTmp[n] = r1 * f + r2 * (1 - f)

        # copy in interploated data
        self.wavl = copy(wavlTmp)
        self.refl = copy(reflTmp)

    def trim(self, wlmin, wlmax):
        """Trim the spectra so it is between two specified wavelengths. Destroys the original data.

        :param wlmin: The lowest wavelength of the new spectra.
        :type wlmin: float
        :param wlmax: The highest wavelength of the new spectra.
        :type wlmax: float
        """

        reflTmp = np.array([])
        wavlTmp = np.array([])

        for (n, wavl) in enumerate(self.wavl):
            if (wavl >= wlmin - 1e-09) and (wavl <= wlmax + 1e-09):
                wavlTmp2 = np.append(wavlTmp, self.wavl[n])
                reflTmp2 = np.append(reflTmp, self.refl[n])
                wavlTmp = copy(wavlTmp2)
                reflTmp = copy(reflTmp2)
        # copy over trimmed data
        self.wavl = copy(wavlTmp)
        self.refl = copy(reflTmp)


def convolve(s1orig, s2orig, resln=1.0, s2norm=True):
    """Convolve one spectra with another, for example
    to apply a band pass, or a spectral response function.

    :param s1orig: A spectra object.
    :type s1orig: object
    :param s2orig: A spectra object.
    :type s2orig: float
    :param resln: The spectral resolution to use.
    :type resln: float
    :param s2norm: If True normalise the second spectra (e.g. to apply a spectra response function).
    :type s2norm: bool
    :return: Convolved spectra.
    :rtype: object
    """

    # make copies so as not to alter
    # original data
    s1 = copy(s1orig)
    s2 = copy(s2orig)

    # interpolate to common resolution
    s1.interpolate(resln)
    s2.interpolate(resln)

    # trim spectra to encompass the exclusive
    # range of the two
    wlmin = np.max([s1.wavl[0], s2.wavl[0]])
    wlmax = np.min([s1.wavl[-1], s2.wavl[-1]])
    s1.trim(wlmin, wlmax)
    s2.trim(wlmin, wlmax)

    # convolve and normailse if required
    norm = 1.0
    if s2norm:
        norm = s2.refl.sum()
    return np.dot(s1.refl, s2.refl) / norm


def sentinel2(s, mission="a"):
    srf = []
    sen2 = spectra()

    if mission == "a":
        fname = "../srfData/S2a_SRF.csv"
    else:
        fname = "../srfData/S2b_SRF.csv"

    for n in xrange(1, 14):
        srf = spectra(fname=fname, ftype="CSV", reflCol=n)
        reflTmp = np.append(sen2.refl, convolve(s, srf, 0.1))
        sen2.refl = copy(reflTmp)
        m = np.argmax(srf.refl)
        wavlTmp = np.append(sen2.wavl, srf.wavl[m])
        sen2.wavl = copy(wavlTmp)

    return sen2


if __name__ == "__main__":

    from matplotlib import pyplot as plt

    doTest1 = True
    doTest2 = False

    if doTest1:
        # test simulation of S2 bands
        svc = spectra(fname="../testData/HRPDA.053017.0065_moc.sig")
        sen2a = sentinel2(svc)
        sen2b = sentinel2(svc, mission="b")

        plt.plot(svc.wavl, svc.refl)
        plt.plot(sen2a.wavl, sen2a.refl, 'o-')
        plt.plot(sen2b.wavl, sen2b.refl, 'o-')
        plt.xlabel('wavelength (nm)')
        plt.ylabel('relfectance (-)')
        plt.show()

    if doTest2:
        # test interpolation and trim routines
        s = spectra(fname="testData/HRPDA.053017.0065_moc.sig")
        plt.plot(s.wavl, s.refl)
        s.interpolate(50)
        plt.plot(s.wavl, s.refl, 'o')
        s.interpolate(100)
        s.trim(1200, 2000)
        plt.plot(s.wavl, s.refl, '--')
        plt.show()
