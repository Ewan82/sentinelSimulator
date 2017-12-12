import netCDF4 as nc
import datetime as dt
import numpy as np
import pandas as pd


class stateVector:
    """Class to hold state vector data for optical
    and microwave canopy RT models.

    """
    def __init__(self):

        self.date_utc = None
        self.lai = None
        self.can_height = None
        self.leafChl = None
        self.leafWater = None
        self.soil_moisture = None
        self.soilAlbedo = None


def get_jules_state(date_utc, nc_file='jules/output/wallerfing_79_12.3_hourly.nc'):
    """Function that returns a stateVector instance for a given time.

    :param date_utc: datetime object of when to extract JULES output.
    :type date_utc: object
    :param nc_file: JULES output file from which to extract data.
    :type nc_file: str
    :return: Instance of stateVector class.
    :rtype: instance
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    t_idx = nc.date2index(date_utc, nc_dat.variables['time'], select='nearest')
    state_inst = stateVector()
    state_inst.date_utc = nc.num2date(nc_dat.variables['time'][t_idx], nc_dat.variables['time'].units)
    state_inst.lai = nc_dat.variables['croplai'][t_idx, 0, 0, 0]  # (m2 m-2)
    state_inst.can_height = nc_dat.variables['cropcanht'][t_idx, 0, 0, 0]  # (m)
    state_inst.soil_moisture = nc_dat.variables['smcl'][t_idx, 0, 0, 0]  # (kg m-2)
    nc_dat.close()
    return state_inst


def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))


def get_date_list(year, month=1, days=365):
    start_date = dt.datetime(year, month, 1, 12, 0)
    date_list = pd.date_range(start_date, periods=days).tolist()
    return date_list


def read(file_format='jules', file_str=None, year=None):
    """Reads output data to a dictionary of state vectors indexed by time.

    .. note:: This function requires sub-functions capable of reading specified file format.

    :param file_format: format of output to read.
    :type file_str: str
    :param file_str: location of file.
    :type file_str: str
    :param year: year of data to extract, if equal to None whole time series extracted
    :type year: int
    :return: state dictionary.
    :rtype: dict
    """
    if file_format == 'jules':
        state_dict = read_jules(file_str, year)
    else:
        state_dict = {}
    return state_dict


def read_jules(nc_file=None, year=None):
    """Reads jules output from netCDF file and writes it to a dictionary indexed by date.

    :param nc_file: location of nc_file.
    :type nc_file: str
    :param year: year of data to extract, if equal to None whole time series extracted.
    :type year: int
    :return: state dictionary.
    :rtype: dict
    """
    if nc_file is None:
        nc_file = 'jules/output/wallerfing_jules_1989_2012.nc'
        print("%s") % nc_file
    nc_dat = nc.Dataset(nc_file, 'r')
    state_dict = {}
    time = nc_dat.variables['time']
    if year is not None:
        strt_idx = nc.date2index(dt.datetime(year, 1, 1), time)
        end_idx = nc.date2index(dt.datetime(year, 12, 31), time)
        t_idx = np.arange(strt_idx, end_idx + 1)
        times = nc.num2date(time[strt_idx:end_idx], time.units)
    else:
        times = nc.num2date(time[:], time.units)
        t_idx = np.arange(len(times))
    for t in enumerate(times):
        state_dict[t[1]] = stateVector()
        state_dict[t[1]].lai = nc_dat.variables['croplai'][t_idx[t[0]], 0, 0, 0]  # (m2 m-2)
        state_dict[t[1]].can_height = nc_dat.variables['cropcanht'][t_idx[t[0]], 0, 0, 0]  # (m)
        state_dict[t[1]].soil_moisture = nc_dat.variables['smcl'][t_idx[t[0]], 0, 0, 0]  # (kg m-2)
        # state_dict[t[1]].soil_temp = nc_dat.variables['t_soil'][t_idx[t[0]], 0, 0, 0]  # (K)
        # figure out how to add others, which soil albedo to output?
    return state_dict
