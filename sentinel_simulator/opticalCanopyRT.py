# python builtins:
import os
import sys
import subprocess
from copy import copy
from tempfile import mkstemp

# third party imports:
import netCDF4
import numpy as np

# sentienl synergy specific imports:
import spectra as sp
import satelliteGeometry as satGeo
import stateVector as sV


def canopyRTOptical(state, geom, resln=1.0):
    """A python wrapper to the SemiDiscrete optical
    canopy RT model of Nadine Gobron. Runs the
    model for the the whole of its valid spectra
    range at a resolution set by resln.

    :param state: Instance of the stateVector class.
    :type state: instance
    :param geom: Instance of the sensorGeomety class.
    :type geom: instance
    :param resln: the spectral resolution in nm [optional].
    :type resln: float
    :return: Instance of the spectra class.
    :rtype: instance
    """

    spect = sp.spectra()
    spect.wavl = np.arange(400, 2500 + resln, resln)

    # generate a tmp file and write
    # wavelengths and geometry into it
    fd, temp_path = mkstemp(prefix='/tmp/senSyntmp__', text=True)
    tmpFile = os.fdopen(fd, 'w')
    print >> tmpFile, "1 %d" % len(spect.wavl),
    for w in spect.wavl:
        print >> tmpFile, " %f" % w,
    print >> tmpFile, "\n%s %s %s %s" %(geom.vza,geom.vaa,geom.sza,geom.saa)
    tmpFile.close()

    # set up nadim command line
    cmd = "semiDiscrete"
    if state.lai != None:
        cmd = cmd + " -LAI %f -hc %f -rsl1 %f" %(state.lai, state.can_height, 0.2*(1.-0.5*(state.soil_moisture/100.)))
        # CHANGE SOIL MOISTURE IMPLEMENTATION HERE, SHOULD BE IN M3 M-3 NOT KG M-2!!!!!
    cmd = cmd + " < %s" % temp_path

    # run process
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    out = p.stdout.readlines()
    p.wait()

    # clean up
    os.remove(temp_path)

    # process RT model output
    rtOut = out[1].split()
    for i in xrange(4, len(rtOut)):
        reflTmp = np.append(spect.refl, float(rtOut[i]))
        spect.refl = copy(reflTmp)

    return spect


def canopyRTOptical_fast(state, geom, mode='fast'):
    """A python wrapper to the SemiDiscrete optical
    canopy RT model of Nadine Gobron. Runs the
    model for the the whole of its valid spectra
    range at a resolution set by resln.

    :param state: Instance of the stateVector class.
    :type state: instance
    :param geom: Instance of the sensorGeomety class.
    :type geom: instance
    :param mode: Run semiDiscrete in either fast ('fast') or slow ('slow') mode [optional].
    :type resln: str
    :return: Instance of the spectra class.
    :rtype: instance
    """

    # generate a tmp file and write geometry into it
    fd, temp_path = mkstemp(prefix='/tmp/senSyntmp__', text=True)
    tmpFile = os.fdopen(fd, 'w')
    print >> tmpFile, "%s %s %s %s" %(geom.vza,geom.vaa,geom.sza,geom.saa)
    tmpFile.close()

    # set up nadim command line
    if mode == 'fast':
        cmd = "semiD -srf ../srfData/s2a.srf -fast"
    else:
        cmd = "semiD -srf ../srfData/s2a.srf"
    if state.lai != None:
        cmd = cmd + " -LAI %f -hc %f -rsl1 %f" %(state.lai, state.can_height, 0.2*(1.-0.5*(state.soil_moisture/100.)))
        # Think about soil moisture implementation here
    cmd = cmd + " < %s" % temp_path

    # run process
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    out = p.stdout.readlines()
    p.wait()
    return np.array([float(val) for val in out[0].split()])


if __name__ == "__main__":
    """This example opens a test output file from the JULES
    model, reads in LAI data, runs these through the optical RT
    model at given geometry and convoles the resulting spectra
    with Sentinel 2 spectra response functions.
    """

    from matplotlib import pyplot as plt

    # read example LAI data:
    allLai = netCDF4.Dataset('../testData/crop_germ_gl4.day.nc').variables['lai'][-365:, 5, 0, 0]

    # main classes:
    state = sV.stateVector()
    geom = satGeo.sensorGeometry()
    geom.sza = 30.

    # container for output:
    allSpect = []
    allBRF = []

    # orbit revist:
    revist = 10

    # plotting variables:
    xpnts = []
    lgnd = []

    # loop over all states and call the
    # canopy RT model and convolve the
    # retruned spectra to S2 bands
    for (n, L) in enumerate(allLai[::revist]):
        state.lai = L
        spect = canopyRTOptical(state, geom)
        allSpect.append(sp.sentinel2(spect))
        allBRF.append(allSpect[n].refl)
        xpnts.append(n * revist)

    allBRF = np.array(allBRF)

    # sort out legend
    for i in xrange(len(allBRF[0, :])):
        lgnd.append('S2 band %d' % (i + 1))

    # do plots:
    lineObjs = plt.plot(xpnts, allBRF)
    plt.ylabel('reflectance (BRF)')
    plt.xlabel('Day of year')
    plt.xlim([0, 364])
    for i in xrange(7, len(lineObjs)):
        lineObjs[i].set_dashes([3, 1])
    plt.legend(iter(lineObjs), lgnd)
    plt.show()
    # plt.savefig('s2Sim_test1.png')
