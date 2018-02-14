import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import opticalCanopyRT as op_can_rt
import satelliteGeometry as satgeo
import stateVector as sv
import spectra as sp
# import sense code
from sense import model as sense_mod
from sense import soil as sense_soil
from sense import canopy as sense_canopy
# import smac code
from SMAC import smac
# import plot sub-module
import plot_sim as plt_s


class Simulator(object):
    """Class to simulate Sentinel observations over Wallerfing for a given year.

    """
    def __init__(self, year=2012, month=1, days=365,
                 jules_nc='jules/output/sensitivity_runs/crp_g_77_6.8535_0.17727_0.000573.3_hourly.nc',
                 mission="None"):
        """Calculate class attributes for given year, month, number of days and netCDF file with driving data.

        :param year: Start year.
        :type year: int
        :param month: Start month.
        :type month: int
        :param days: number of days to run for.
        :type days: int
        :param jules_nc: location of netCDF file to use
        :type jules_nc: str
        :param mission: Sentinel-1b or Sentinel-2a for satellite geometry
        :type mission: str
        :return: Instance of Simulator class.
        :rtype: object
        """
        #  Setup satellite geometry list
        self.jules_nc = jules_nc
        self.start_date = dt.datetime(year, month, 1)
        self.lon = 12.880
        self.lat = 48.684
        self.alt = 0.
        self.days = days
        if mission == "None":
            self.date_lst_in = sv.get_date_list(year, month, days)
        else:
            self.geom_lst = satgeo.getSentinel2Geometry(self.start_date, self.days, self.lat, self.lon,
                                                        mission=mission, alt=self.alt)
            self.vza_lst = [geo.vza for geo in self.geom_lst]
            self.vaa_lst = [geo.vaa for geo in self.geom_lst]
            self.sza_lst = [geo.sza for geo in self.geom_lst]
            self.saa_lst = [geo.saa for geo in self.geom_lst]
            self.date_lst_in = [geo.date_utc for geo in self.geom_lst]

        #  Setup JULES state vector list
        self.state_lst = [sv.get_jules_state(dat, self.jules_nc) for dat in self.date_lst_in]
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
    def __init__(self, lai_coef=0.1, s=0.015, omega=0.1, **kwargs):
        """
        Initialize with same arguemnts as superClass 'Simulator'.

        """
        super(S1_simulator, self).__init__(mission="Sentinel-1b", **kwargs)
        # Setup SAR RT spectra list (Sentinel 1)
        self.freq = 5.405
        self.theta = np.deg2rad(37.0)
        self.stype = 'turbid_rayleigh'
        #self.stype='turbid_isotropic'
        self.surf = 'Oh92'
        # self.surf = 'Dubois95'
        self.models = {'surface': self.surf, 'canopy': self.stype}
        self.s = s  # 0.02
        self.lai_coef = lai_coef  # 0.1
        self.eps = 15. - 0.j
        omega = omega  # 0.12  single scattering albedo

        self.SAR_list = [sense_mod.SingleScatRT(
            surface=sense_soil.Soil(mv=self.soilm_lst[x], f=self.freq, s=s, clay=0.23, sand=0.27,
                                    bulk=1.65),
            #surface=sense_soil.Soil(eps=self.eps, f=self.freq, s=self.s),
            canopy=sense_canopy.OneLayer(ke_h=self.lai_coef*self.state_lst[x].lai,  # extinction coefficient
                                         ke_v=self.lai_coef*self.state_lst[x].lai,
                                         d=self.state_lst[x].can_height,
                                         ks_h=omega * (self.lai_coef*self.state_lst[x].lai),  # volume scattering coeff
                                         ks_v=omega * (self.lai_coef*self.state_lst[x].lai)),
            models=self.models,
            # theta=np.deg2rad(self.vza_lst[x]),
            theta = self.theta,
            freq=self.freq) for x in xrange(len(self.state_lst))]
        for s in self.SAR_list:
            s.sigma0()
        self.backscatter_keys = ['vv', 'hh', 'hv']
        # Extract total backscatter (and convert to dB) in the three polarisations from SAR list
        self.stot_hv = np.array([10 * np.log10(self.SAR_list[s1c].__dict__['stot']['hv'])
                                for s1c in xrange(len(self.SAR_list))])
        self.stot_hh = np.array([10 * np.log10(self.SAR_list[s1c].__dict__['stot']['hh'])
                                for s1c in xrange(len(self.SAR_list))])
        self.stot_vv = np.array([10 * np.log10(self.SAR_list[s1c].__dict__['stot']['vv'])
                                for s1c in xrange(len(self.SAR_list))])


class S2_simulator(Simulator):
    """Given Simulator class this subclass will simulate Sentinel 2 data.

    """
    def __init__(self, **kwargs):
        """
        Initialize with same arguemnts as superClass 'Simulator'.

        """
        super(S2_simulator, self).__init__(mission="Sentinel-2a", **kwargs)
        # Setup canopy optical RT spectra list (Sentinel 2)
        self.all_BRF_arr = np.array([op_can_rt.canopyRTOptical_fast(self.state_lst[x], self.geom_lst[x]) for x in
                          xrange(len(self.state_lst))])
        self.band_labels = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8a', 'B9', 'B10', 'B11', 'B12']


class S2_simulator_old_semid(Simulator):
    """Given Simulator class this subclass will simulate Sentinel 2 data.

    """
    def __init__(self, **kwargs):
        """
        Initialize with same arguemnts as superClass 'Simulator'.

        """
        super(S2_simulator_old_semid, self).__init__(mission="Sentinel-2a", **kwargs)
        # Setup canopy optical RT spectra list (Sentinel 2)
        self.spect_lst = [op_can_rt.canopyRTOptical(self.state_lst[x], self.geom_lst[x]) for x in
                          xrange(len(self.state_lst))]
        self.all_spect_lst = [sp.sentinel2(spect) for spect in self.spect_lst]
        self.all_BRF_arr = np.array([all_sp.refl for all_sp in self.all_spect_lst])
        self.band_labels = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8a', 'B9', 'B10', 'B11', 'B12']


if __name__ == "__main__":
    """Run simulator and plot output.
    """
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 2, 'lines.markersize': 6})
    sns.set_style('whitegrid')
    sim_c1 = S1_simulator()  # lai_coef=0.01, s=0.006)
    sim_c2 = S2_simulator()
    # sim_c2_slow = S2_simulator_old_semid()
    output_dir = '../docs/source/simulator/test_fast/'

    fig = plt_s.plot_class_var(sim_c2.date_lst, sim_c2.vza_lst, y_lab='View zenith angle (degrees)', line_type='o')[0]
    fig.savefig(output_dir + 'vza_s2.png')
    fig = plt_s.plot_class_var(sim_c2.date_lst, sim_c2.sza_lst, y_lab='Solar zenith angle (degrees)', line_type='o')[0]
    fig.savefig(output_dir + 'sza_s2.png')
    fig = plt_s.plot_class_var(sim_c1.date_lst, sim_c1.vza_lst, y_lab='View zenith angle (degrees)', line_type='o')[0]
    fig.savefig(output_dir + 'vza_s1.png')
    fig = plt_s.plot_class_var(sim_c1.date_lst, sim_c1.sza_lst, y_lab='Solar zenith angle (degrees)', line_type='o')[0]
    fig.savefig(output_dir + 'sza_s1.png')
    fig = plt_s.plot_class_var(sim_c1.date_lst, sim_c1.lai_lst, y_lab='Leaf area index')[0]
    fig.savefig(output_dir + 'lai.png')
    fig = plt_s.plot_class_var(sim_c1.date_lst, sim_c1.canht_lst, y_lab='Canopy height (m)')[0]
    fig.savefig(output_dir + 'can_ht.png')
    fig = plt_s.plot_class_var(sim_c1.date_lst, sim_c1.soilm_lst, y_lab=r'Soil moisture (m$^{3}~$m$^{-3}$)')[0]
    fig.savefig(output_dir + 'soil_m.png')
    plt.close('all')
    sig_list = ['stot', 's0cgt', 's0c', 's0gcg', 's0g']
    sig = 'stot'  # for sig in sig_list:
    for x in xrange(3):
        fig = plt_s.plot_class_var(sim_c1.date_lst,
                             [10*np.log10(sim_c1.SAR_list[s1c].__dict__[sig][sim_c1.backscatter_keys[x]])
                              for s1c in xrange(len(sim_c1.SAR_list))],
                             y_lab='Backscatter ' + sim_c1.backscatter_keys[x] + ' polarisation (db)',
                             line_type='o')[0]
        # Must also think about canopy height and extinction coefficient!!!
        fig.savefig(output_dir + sim_c1.backscatter_keys[x] + '.png')
        plt.close()

    for x in range(13):
        fig = plt_s.plot_class_var(sim_c2.date_lst, sim_c2.all_BRF_arr[:,x], y_lab=sim_c2.band_labels[x] + ' reflectance',
                             line_type='o')[0]
        fig.savefig(output_dir + sim_c2.band_labels[x] + '.png')
        plt.close()
        coeffs_file = 'SMAC/COEFS/Coef_S2A_CONT_'+sim_c2.band_labels[x]+'.dat'
        coeffs = smac.coeff(coeffs_file)
        r_toa = smac.smac_dir(sim_c2.all_BRF_arr[:,x], np.mean(sim_c2.sza_lst), np.mean(sim_c2.saa_lst),
                              np.mean(sim_c2.vza_lst), np.mean(sim_c2.vaa_lst),
                              1013, 0.1, 0.3, 0.3, coeffs)
        fig = plt_s.plot_class_var(sim_c2.date_lst, r_toa, y_lab=sim_c2.band_labels[x]+' TOA reflectance',
                             line_type='o')[0]
        fig.savefig(output_dir + sim_c2.band_labels[x] + '_TOA.png')
        plt.close()