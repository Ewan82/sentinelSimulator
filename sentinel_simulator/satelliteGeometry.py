from pyorbital.orbital import Orbital
from datetime import datetime, time
import pyorbital.astronomy as ast
import numpy as np


class sensorGeometry:
    """Class to hold sun-sensor geometry information.
    """
    def __init__(self):
        self.date_utc = None
        # view zenith angle
        self.vza = None
        # view azimuth angle
        self.vaa = None
        # solar zenith angle
        self.sza = None
        # solar azimuth angle
        self.saa = None

    def printGeom(self):
        """Prints currently specified class attributes.
        """
        print(self.date_utc, self.vza, self.vaa, self.sza, self.saa)


def getSentinel2Geometry(startDateUTC, lengthDays, lat, lon, alt=0.0, mission="Sentinel-2a",
                         tleFile="../TLE/norad_resource_tle.txt"):
    """Calculate approximate geometry for Sentinel overpasses.
    Approximate because it assumes maximum satellite elevation
    is the time at which target is imaged.

    :param startDateUTC: a datetime object specifying when to start prediction.
    :type startDateUTC: object
    :param lengthDays: number of days over which to perform calculations.
    :type lengthDays: int
    :param lat: latitude of target.
    :type lat: float
    :param lon: longitude of target.
    :type lon: float
    :param alt: altitude of target (in km).
    :type alt: float
    :param mission: mission name as in TLE file.
    :type mission: str
    :param tleFile: TLE file.
    :type tleFile: str
    :return: a python list containing instances of the sensorGeometry class arranged in date order.
    :rtype: list
    """
    orb = Orbital(mission, tleFile)
    passes = orb.get_next_passes(startDateUTC, 24 * lengthDays, lon, lat, alt)

    geomList = []

    for p in passes:
        look = orb.get_observer_look(p[2], lon, lat, alt)
        vza = 90 - look[1]
        vaa = look[0]
        sza = ast.sun_zenith_angle(p[2], lon, lat)
        saa = np.rad2deg(ast.get_alt_az(p[2], lon, lat)[1])

        if mission == 'Sentinel-1b':
            if 75 < vaa < 105 and 20. < vza < 45.:  # vaa usually [0, 180], testing observation times
                thisGeom = sensorGeometry()
                thisGeom.date_utc = p[2]
                thisGeom.vza = vza
                thisGeom.vaa = vaa
                thisGeom.sza = sza
                thisGeom.saa = saa
                geomList.append(thisGeom)
        elif mission == 'Sentinel-1a':
            if 75 < vaa < 105 and 20. < vza < 45.:  # vaa usually [0, 180], testing observation times
                thisGeom = sensorGeometry()
                thisGeom.date_utc = p[2]
                thisGeom.vza = vza
                thisGeom.vaa = vaa
                thisGeom.sza = sza
                thisGeom.saa = saa
                geomList.append(thisGeom)
        elif mission == 'Sentinel-2a':
            if sza < 90 and vza < 10.3:
                thisGeom = sensorGeometry()
                thisGeom.date_utc = p[2]
                thisGeom.vza = vza
                thisGeom.vaa = vaa
                thisGeom.sza = sza
                thisGeom.saa = saa
                geomList.append(thisGeom)

    return geomList


if __name__ == "__main__":
    """Simple test with S2a for Wallerfing
    """

    startDate = datetime(2017, 7, 1, 0, 0, 0, 0, None)
    lon, lat = 12.880, 48.684
    alt = 0.
    days = 30

    geomList = getSentinel2Geometry(startDate, days, lat, lon, mission="Sentinel-1a", alt=alt)
    for g in geomList:
        g.printGeom()
