import numpy as np
import simulator as sim
import pickle
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import plot_sim as plt_s
from SMAC import smac


def S1_run(jules_nc='jules/output/sensitivity_runs/crp_g_77_6.8535_0.17727_0.000573.3_hourly.nc',
           fname='s1_sim_dict.p', year=2012, month=1, days=364, lai_coeff=0.1, s=0.015):
    simulator = sim.Simulator(jules_nc=jules_nc, year=year, month=month, days=days)
    S1_simulator = sim.S1_simulator(lai_coef=lai_coeff, s=s, omega=0.1, jules_nc=jules_nc, year=year, month=month,
                                 days=days)
    sim_dict = {}
    sim_dict['nc_file'] = jules_nc
    sim_dict['S1_dates'] = S1_simulator.date_lst
    sim_dict['dates'] = simulator.date_lst
    sim_dict['lai'] = simulator.lai_lst
    sim_dict['can_height'] = simulator.canht_lst
    sim_dict['soil_m'] = simulator.soilm_lst
    for x in xrange(3):
        sim_dict[S1_simulator.backscatter_keys[x]] = \
            [10*np.log10(S1_simulator.SAR_list[s1c].__dict__['stot'][S1_simulator.backscatter_keys[x]])
             for s1c in xrange(len(S1_simulator.SAR_list))]

    f = open(fname, 'wb')
    pickle.dump(sim_dict, f)
    f.close()
    return 'simulator pickled'


def S2_run(jules_nc='jules/output/sensitivity_runs/crp_g_77_6.8535_0.17727_0.000573.3_hourly.nc',
           fname='s2_sim_dict.p', year=2012, month=1, days=364):
    simulator = sim.Simulator(jules_nc=jules_nc, year=year, month=month, days=days)
    S2_simulator = sim.S2_simulator(jules_nc=jules_nc, year=year, month=month, days=days)
    sim_dict = {}
    sim_dict['nc_file'] = jules_nc
    sim_dict['S2_dates'] = S2_simulator.date_lst
    sim_dict['dates'] = simulator.date_lst
    sim_dict['lai'] = simulator.lai_lst
    sim_dict['can_height'] = simulator.canht_lst
    sim_dict['soil_m'] = simulator.soilm_lst
    for x in range(13):
        sim_dict[S2_simulator.band_labels[x]] = S2_simulator.all_BRF_arr[:,x]
    f = open(fname, 'wb')
    pickle.dump(sim_dict, f)
    f.close()
    return 'simulator pickled'


def clim_sen_S1(save_dir):
    for yr in xrange(1979, 2012+1):
        try:
            S1_run(jules_nc='jules/output/wallerfing_79_12.3_hourly.nc', fname=save_dir+'/'+str(int(yr))+'_s1_sim_dict.p',
                   year=yr)
        except ValueError:
            print 'caught a value error, doing next year %d' % (yr)
            continue
    return 'finished climatological sensitivity run'


def clim_sen_S2(save_dir):
    for yr in xrange(1979, 2012+1):
        try:
            S2_run(jules_nc='jules/output/wallerfing_79_12.3_hourly.nc', fname=save_dir+'/'+str(int(yr))+'_s2_sim_dict.p',
               year=yr)
        except ValueError:
            print 'caught a value error, doing next year %d' %(yr)
            continue
    return 'finished climatological sensitivity run'


def param_sen_S1(save_dir):
    glob_list = pickle.load(open('glob_list.p', 'rb'))
    for fname in glob_list:
        try:
            S1_run(jules_nc=fname, fname=save_dir+'/'+fname[30:46]+'_s1_sim_dict.p')
        except ValueError:
            print 'caught a value error, doing next year'
            continue
    return 'finished parameter sensitivity run'


def param_sen_S2(dir):
    glob_list = pickle.load(open('glob_list.p', 'rb'))
    for fname in glob_list:
        try:
            S2_run(jules_nc=fname, fname=dir+'/'+fname[30:46]+'_s2_sim_dict.p')
        except ValueError:
            print 'caught a value error, doing next year'
            continue
    return 'finished parameter sensitivity run'


def demo_sen_S1(direct):
    glob_list = glob.glob(direct + '/*.nc')
    for fname in glob_list:
        try:
            S1_run(jules_nc=fname, fname=direct + '/' + fname[-16:-13] + '_s1_sim_dict.p', lai_coeff=0.1, s=0.015)
        except ValueError:
            print 'caught a value error, doing next year'
            continue
    return 'finished parameter sensitivity run'


def demo_sen_S2(direct):
    glob_list = glob.glob(direct + '/*.nc')
    for fname in glob_list:
        try:
            S2_run(jules_nc=fname, fname=direct + '/' + fname[-16:-13] + '_s2_sim_dict.p')
        except ValueError:
            print 'caught a value error, doing next year'
            continue
    return 'finished parameter sensitivity run'


def plot_sens_s1(pik_dir, save_dir):
    glob_list = glob.glob(pik_dir+'/*s1*.p')
    glob_0 = pickle.load(open(glob_list[0], 'rb'))
    date_lst = glob_0['dates']
    s1_date_lst = glob_0['S1_dates']
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 2, 'lines.markersize': 6})
    sns.set_style('whitegrid')

    fig1, ax1 = plt_s.plot_class_var(date_lst, glob_0['lai'], y_lab='Leaf area index')
    fig2, ax2 = plt_s.plot_class_var(date_lst, glob_0['can_height'], y_lab='Canopy height (m)')
    fig3, ax3 = plt_s.plot_class_var(date_lst, glob_0['soil_m'], y_lab=r'Soil moisture (m$^{3}~$m$^{-3}$)')

    for x in xrange(1, len(glob_list)):
        print glob_list[x]
        glob_dict = pickle.load(open(glob_list[x], 'rb'))
        print len(glob_dict['lai'])
        len_plt = len(glob_dict['dates'])
        if len_plt > len(glob_0['dates']):
            len_plt = len(glob_0['dates'])
        plt_s.plot_class_var(date_lst[0:len_plt], glob_dict['lai'][0:len_plt], y_lab='Leaf area index', axes=ax1)
        plt_s.plot_class_var(date_lst[0:len_plt], glob_dict['can_height'][0:len_plt], y_lab='Canopy height (m)', axes=ax2)
        plt_s.plot_class_var(date_lst[0:len_plt], glob_dict['soil_m'][0:len_plt], y_lab=r'Soil moisture (m$^{3}~$m$^{-3}$)', axes=ax3)
    fig1.savefig(save_dir + '/lai.png')
    fig2.savefig(save_dir + '/can_ht.png')
    fig3.savefig(save_dir + '/soil_m.png')
    plt.close('all')
    sig_list = ['stot', 's0cgt', 's0c', 's0gcg', 's0g']
    sig = 'stot'  # for sig in sig_list:
    pol = ['vv', 'hh', 'hv']
    for x in xrange(3):
        fig, ax = plt_s.plot_class_var(s1_date_lst, glob_0[pol[x]], y_lab='Backscatter ' + pol[x] + ' polarisation (db)',
                                     line_type='o')
        # Must also think about canopy height and extinction coefficient!!!
        for i in xrange(1, len(glob_list)):
            glob_dict = pickle.load(open(glob_list[i], 'rb'))
            len_plt = len(glob_dict['S1_dates'])
            if len_plt > len(glob_0['S1_dates']):
                len_plt = len(glob_0['S1_dates'])
            plt_s.plot_class_var(s1_date_lst[0:len_plt], glob_dict[pol[x]][0:len_plt], y_lab='Backscatter ' + pol[x] + ' polarisation (db)', line_type='o',
                               axes=ax)
        fig.savefig(save_dir + '/' + pol[x] + '.png')
        plt.close()


def plot_sens_s2(pik_dir, save_dir):
    glob_list = glob.glob(pik_dir+'/*s2*.p')
    glob_0 = pickle.load(open(glob_list[0], 'rb'))
    date_lst = glob_0['dates']
    s2_date_lst = glob_0['S2_dates']
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 2, 'lines.markersize': 6})
    sns.set_style('whitegrid')

    fig1, ax1 = plt_s.plot_class_var(date_lst, glob_0['lai'], y_lab='Leaf area index')
    fig2, ax2 = plt_s.plot_class_var(date_lst, glob_0['can_height'], y_lab='Canopy height (m)')
    fig3, ax3 = plt_s.plot_class_var(date_lst, glob_0['soil_m'], y_lab=r'Soil moisture (m$^{3}~$m$^{-3}$)')
    print len(date_lst)
    for x in xrange(1, len(glob_list)):
        print glob_list[x]
        glob_dict = pickle.load(open(glob_list[x], 'rb'))
        print len(glob_dict['lai'])
        len_plt = len(glob_dict['dates'])
        if len_plt > len(glob_0['dates']):
            len_plt = len(glob_0['dates'])
        plt_s.plot_class_var(date_lst[0:len_plt], glob_dict['lai'][0:len_plt], y_lab='Leaf area index', axes=ax1)
        plt_s.plot_class_var(date_lst[0:len_plt], glob_dict['can_height'][0:len_plt], y_lab='Canopy height (m)', axes=ax2)
        plt_s.plot_class_var(date_lst[0:len_plt], glob_dict['soil_m'][0:len_plt], y_lab=r'Soil moisture (m$^{3}~$m$^{-3}$)', axes=ax3)
    fig1.savefig(save_dir + '/lai.png')
    fig2.savefig(save_dir + '/can_ht.png')
    fig3.savefig(save_dir + '/soil_m.png')
    plt.close('all')
    band_labels = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8a', 'B9', 'B10', 'B11', 'B12']
    band_labels2 = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8a', 'B9', 'B10', 'B11', 'B12']
    for x in xrange(13):
        fig, ax = plt_s.plot_class_var(s2_date_lst, glob_0[band_labels[x]], y_lab=band_labels[x] + ' reflectance',
                                     line_type='o')
        # Must also think about canopy height and extinction coefficient!!!
        for i in xrange(1, len(glob_list)):
            glob_dict = pickle.load(open(glob_list[i], 'rb'))
            len_plt = len(glob_dict['S2_dates'])
            if len_plt > len(glob_0['S2_dates']):
                len_plt = len(glob_0['S2_dates'])
            plt_s.plot_class_var(s2_date_lst[0:len_plt], glob_dict[band_labels[x]][0:len_plt],
                               y_lab=band_labels[x] + ' reflectance', line_type='o', axes=ax)
        fig.savefig(save_dir + '/' + band_labels[x] + '.png')
        plt.close()
    for x in xrange(13):
        coeffs_file = 'SMAC/COEFS/Coef_S2A_CONT_' + band_labels2[x] + '.dat'
        coeffs = smac.coeff(coeffs_file)
        fig, ax = plt_s.plot_class_var(s2_date_lst, smac.smac_dir(glob_0[band_labels[x]], 49.17, 163.18, 5.14, 195.5,
                                     1013, 0.1, 0.3, 0.3, coeffs), y_lab=band_labels[x] + ' TOA reflectance',
                                     line_type='o')
        # Must also think about canopy height and extinction coefficient!!!
        for i in xrange(1, len(glob_list)):
            glob_dict = pickle.load(open(glob_list[i], 'rb'))
            len_plt = len(glob_dict['S2_dates'])
            if len_plt > len(glob_0['S2_dates']):
                len_plt = len(glob_0['S2_dates'])
            plt_s.plot_class_var(s2_date_lst[0:len_plt], smac.smac_dir(glob_dict[band_labels[x]][0:len_plt], 49.17,
                                                                                  163.18, 5.14, 195.5, 1013, 0.1, 0.3,
                                                                                  0.3, coeffs),
                               y_lab=band_labels[x] + ' TOA reflectance', line_type='o', axes=ax)
        fig.savefig(save_dir + '/' + band_labels[x] + '_TOA.png')
        plt.close()

if __name__ == "__main__":
    """Run simulator and plot output.
    """