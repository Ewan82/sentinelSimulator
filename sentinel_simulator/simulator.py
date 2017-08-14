
import numpy as np
import netCDF4 as nc
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

import opticalCanopyRT as op_can_rt
import satelliteGeometry as satgeo
import stateVector as sv
import spectra as sp


class Simulator:
    """Class to simulate Sentinel 2 observations over Wallerfing for a given year.
    """
    def __init__(self, year=2012, month=1, days=365):
        """Calculate class attributes for given year, month and number of days.

        :param year: Start year.
        :type year: int
        :param month: Start month.
        :type month: int
        :param days: number of days to run for.
        :type days: int
        """
        #  Setup satellite geometry list
        self.start_date = dt.datetime(year, month, 1)
        self.lon = 12.880
        self.lat = 48.684
        self.alt = 0.
        self.days = days
        self.geom_lst = satgeo.getSentinel2Geometry(self.start_date, self.days, self.lat, self.lon,
                                                    mission="Sentinel-2a", alt=self.alt)
        self.vza_lst = [geo.vza for geo in self.geom_lst]
        self.sza_lst = [geo.sza for geo in self.geom_lst]

        #  Setup JULES state vector list
        self.state_lst = [sv.get_jules_state(geo.date_utc, 'jules/output/sensitivity_runs/crp_g_77_6.8535_0.17727_'
                                                           '0.000573.3_hourly.nc') for geo in self.geom_lst]
        self.date_lst = [state.date_utc for state in self.state_lst]
        self.lai_lst = [state.lai for state in self.state_lst]
        self.canht_lst = [state.can_height for state in self.state_lst]
        self.soilm_lst = [state.soil_moisture / 100. for state in self.state_lst]

        #  Setup canopy optical RT spectra list
        self.spect_lst = [op_can_rt.canopyRTOptical(self.state_lst[x], self.geom_lst[x]) for x in
                          xrange(len(self.state_lst))]
        self.all_spect_lst = [sp.sentinel2(spect) for spect in self.spect_lst]
        self.all_BRF_arr = np.array([all_sp.refl for all_sp in self.all_spect_lst])
        self.band_labels = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']


def plot_class_var(date_lst, var, y_lab=None, line_type='-'):
    """Plot specified variable.
    :param var: Class attribute variable as list.
    :type var: list
    :param y_lab: Label for Y-axis.
    :type y_lab: str
    :return: Figure.
    :rtype: object
    """
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 1, 'lines.markersize': 10})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 5))
    #ax.xaxis_date()
    sns.set_style('whitegrid')
    ax.plot(date_lst, var, line_type)
    plt.ylabel(y_lab)
    plt.xlabel('Date')
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%B')
    ax.xaxis.set_major_formatter(myFmt)
    # plt.legend(loc=2)
    # plt.show()
    return fig


if __name__ == "__main__":
    """Run simulator and plot output.
    """
    sim_c = Simulator()

    fig = plot_class_var(sim_c.date_lst, sim_c.vza_lst, y_lab='View zenith angle (degrees)')
    fig.savefig('../docs/source/simulator/vza.png')
    fig = plot_class_var(sim_c.date_lst, sim_c.sza_lst, y_lab='Solar zenith angle (degrees)')
    fig.savefig('../docs/source/simulator/sza.png')
    fig = plot_class_var(sim_c.date_lst, sim_c.lai_lst, y_lab='Leaf area index')
    fig.savefig('../docs/source/simulator/lai.png')
    fig = plot_class_var(sim_c.date_lst, sim_c.canht_lst, y_lab='Canopy height (m)')
    fig.savefig('../docs/source/simulator/can_ht.png')
    fig = plot_class_var(sim_c.date_lst, sim_c.soilm_lst, y_lab=r'Soil moisture (m$^{3}~$m$^{-3}$)')
    fig.savefig('../docs/source/simulator/soil_m.png')
    for x in range(13):
        fig = plot_class_var(sim_c.date_lst, sim_c.all_BRF_arr[:,x], y_lab=sim_c.band_labels[x]+' reflectance')
        fig.savefig('../docs/source/simulator/'+sim_c.band_labels[x]+'.png')