import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import opticalCanopyRT as op_can_rt
import satelliteGeometry as satgeo
import stateVector as sv
import spectra as sp
# Import sense code
from sense import model as sense_mod
from sense import soil as sense_soil
from sense import canopy as sense_canopy

class Simulator(object):
    """Class to simulate Sentinel observations over Wallerfing for a given year.

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


class S1_simulator(Simulator):
    """Given Simulator class this subclass will simulate Sentinel 1 data.

    .. note:: This function requires the Community SAR ScattEring model (SenSE) to be installed on the system. This
     code is available from:

     - https://github.com/PMarzahn/sense

    """
    def __init__(self, **kwargs):
        """
        Initialize with same arguemnts as superClass 'Simulator'.

        """
        super(S1_simulator, self).__init__(**kwargs)
        # Setup SAR RT spectra list (Sentinel 1)
        self.freq = 5.
        self.theta = np.deg2rad(37.0)
        self.stype = 'turbid_rayleigh'
        # stype='turbid_isotropic'
        self.models = {'surface': 'Oh92', 'canopy': self.stype}
        self.s = 0.02
        d = 0.22
        ke = 1.
        # canopy = sense_canopy.OneLayer(ke_h=ke, ke_v=ke, d=0.1*self.state_lst[x].can_height, ks_h=omega * ke, ks_v=omega * ke)
        omega = 0.12

        self.SAR_list = [sense_mod.SingleScatRT(
            surface=sense_soil.Soil(mv=self.state_lst[x].soil_moisture, f=self.freq, s=self.s, clay=0.23, sand=0.27),
            canopy=sense_canopy.OneLayer(ke_h=0.1*self.state_lst[x].lai, ke_v=0.1*self.state_lst[x].lai,
                                         d=0.1*self.state_lst[x].can_height,
                                         ks_h=omega * 0.1*self.state_lst[x].lai,
                                         ks_v=omega * 0.1*self.state_lst[x].lai),
            models=self.models,
            theta=self.theta,
            freq=self.freq) for x in xrange(len(self.state_lst))]
        for s in self.SAR_list:
            s.sigma0()
        self.backscatter_keys = ['vv', 'hh', 'hv']


class S2_simulator(Simulator):
    """Given Simulator class this subclass will simulate Sentinel 2 data.

    """
    def __init__(self, **kwargs):
        """
        Initialize with same arguemnts as superClass 'Simulator'.

        """
        super(S2_simulator, self).__init__(**kwargs)
        # Setup canopy optical RT spectra list (Sentinel 2)
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
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 2, 'lines.markersize': 6})
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
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 2, 'lines.markersize': 6})
    sns.set_style('whitegrid')
    sim_c1 = S1_simulator()
    #sim_c2 = S2_simulator()

    fig = plot_class_var(sim_c1.date_lst, sim_c1.vza_lst, y_lab='View zenith angle (degrees)', line_type='o')
    fig.savefig('../docs/source/simulator/vza.png')
    fig = plot_class_var(sim_c1.date_lst, sim_c1.sza_lst, y_lab='Solar zenith angle (degrees)', line_type='o')
    fig.savefig('../docs/source/simulator/sza.png')
    fig = plot_class_var(sim_c1.date_lst, sim_c1.lai_lst, y_lab='Leaf area index')
    fig.savefig('../docs/source/simulator/lai.png')
    fig = plot_class_var(sim_c1.date_lst, sim_c1.canht_lst, y_lab='Canopy height (m)')
    fig.savefig('../docs/source/simulator/can_ht.png')
    fig = plot_class_var(sim_c1.date_lst, sim_c1.soilm_lst, y_lab=r'Soil moisture (m$^{3}~$m$^{-3}$)')
    fig.savefig('../docs/source/simulator/soil_m.png')
    for x in xrange(3):  # Currently using vv polarisation (vv, hh, hv) is this correct?
        fig = plot_class_var(sim_c1.date_lst,
                             [10*np.log10(sim_c1.SAR_list[s1c].__dict__['stot'][sim_c1.backscatter_keys[x]])
                              for s1c in xrange(len(sim_c1.SAR_list))],
                             y_lab='Backscatter ' + sim_c1.backscatter_keys[x] + ' polarisation (db)',
                             line_type='o')
        # Must also think about canopy height and extinction coefficient!!!
        fig.savefig('../docs/source/simulator/' + sim_c1.backscatter_keys[x] + '.png')
    #for x in range(13):
    #    fig = plot_class_var(sim_c2.date_lst, sim_c2.all_BRF_arr[:,x], y_lab=sim_c2.band_labels[x]+' reflectance',
    #                         line_type='o')
    #    fig.savefig('../docs/source/simulator/'+sim_c2.band_labels[x]+'.png')