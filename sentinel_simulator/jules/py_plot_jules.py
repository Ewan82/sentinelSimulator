# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 11:26:19 2017

@author: Ewan Pinnington

Plotting JULES results over Germany
"""
import numpy as np
import datetime as dt
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import sys
import glob


def open_nc(f_name):
    """Opens a netCDF dataset.

    :param f_name: location and name of file.
    :type f_name: str
    :return: netCDF dataset object.
    :rtype: object
    """
    return nc.Dataset(f_name, 'r')


def extract_vars_nc(data_nc, var_name, strt_yr=1980, end_yr=2012, strt_day=1, strt_hr=0):
    """Extracts variables from a given netCDF dataset object.

    :param data_nc: netCDF dataset object.
    :type data_nc: object
    :param var_name: Name of variable to be extracted.
    :type var_name: str
    :param strt_yr: Start year of data.
    :type strt_yr: int
    :param end_yr: End year of data.
    :type end_yr: int
    :return: Latitude values, longitude values, variables values, time values.
    :rtype: tuple
    """
    lats = data_nc.variables['latitude'][:,0]
    lons = data_nc.variables['longitude'][0,:]
    time = data_nc.variables['time']
    var = data_nc.variables[var_name][:]
    return lats, lons, var, time


def plot_climatology(nc_file, var_name, level='none', ylab = 'None'):
    """Given Wallerfing climatological JULES model output plots given variable climatology.

    :param nc_file: File location for JULES output.
    :type nc_file: str
    :param var_name: Name of variable to plot.
    :type var_name: str
    :param level: If variable is 4D specify 4th dimension.
    :type level: int
    :param ylab: Label for Y-axis.
    :type ylab: str
    :return: Figure object to save.
    :rtype: object
    """
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 1, 'lines.markersize': 10})
    fig, ax = plt.subplots(nrows=1, ncols=1,) # figsize=(15, 5))
    sns.set_style('whitegrid')
    palette = sns.color_palette("colorblind", 11)
    dat = open_nc(nc_file)
    lats, lons, var, time = extract_vars_nc(dat, var_name)
    times = nc.num2date(time[:], time.units)
    idx = np.where([times[x].year == times[366*8].year for x in range(len(times))])[0]
    time_x = times[idx]
    plt_var = var[:]
    plt_var[plt_var > 1e18] = np.nan
    depths = [100, 250, 650, 2000]
    #depths = [150, 350, 650, 2000]
    labels = ['0 - 0.1m', '0.1 - 0.35m', '0.35 - 1m', '1 - 3m']
    #if level in [0,1,2,3]:
    #    ax.plot(times[0:365], plt_var[0:365, level]/depths[level], label='wfdei', color=palette[0])
    #else:
    for yr in xrange(times[0].year, times[-1].year):
        idx = np.where([times[x].year == yr for x in range(len(times))])[0]
        #  print len(idx)
        if var_name == 'smcl':
            ax.plot(time_x[0:364*8], plt_var[idx[0]:idx[364*8], level, 0, 0]/depths[level], )
        elif level != 'none':
            ax.plot(time_x[0:364*8], plt_var[idx[0]:idx[364*8], level, 0, 0],)
        else:
            ax.plot(time_x[0:364*8], plt_var[idx[0]:idx[364*8], 0, 0], )
    #plt.ylabel('Volumetric soil water content (m3 m-3)')
    plt.xlabel('Date')
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%B')
    ax.xaxis.set_major_formatter(myFmt)
    if ylab != 'None':
        plt.ylabel(ylab)
    #plt.legend(loc=2)
    #plt.show()
    return fig


def plot_sensitivity(nc_dir, var_name, level='none', ylab = 'None'):
    """Given Wallerfing sensitivity JULES model output plots given variable sensitivity.

    :param nc_dir: Directory location for JULES output.
    :type nc_dir: str
    :param var_name: Name of variable to plot.
    :type var_name: str
    :param level: If variable is 4D specify 4th dimension.
    :type level: int
    :param ylab: Label for Y-axis.
    :type ylab: str
    :return: Figure object to save.
    :rtype: object
    """
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 1, 'lines.markersize': 10})
    fig, ax = plt.subplots(nrows=1, ncols=1,) # figsize=(15, 5))
    sns.set_style('whitegrid')
    palette = sns.color_palette("colorblind", 11)
    for nc_file in glob.glob(nc_dir+'crp_g_*.3_hourly.nc'):
        dat = open_nc(nc_file)
        lats, lons, var, time = extract_vars_nc(dat, var_name, strt_yr=2012, end_yr=2012, strt_day=1, strt_hr=3)
        times = nc.num2date(time[:], time.units)
        plt_var = var[:]
        plt_var[plt_var > 1e18] = np.nan
        depths = [100, 250, 650, 2000]
        if var_name == 'smcl':
            ax.plot(times[:], plt_var[:, level, 0, 0]/depths[level],)
        elif level != 'none':
            ax.plot(times[:], plt_var[:, level, 0, 0],)
        else:
            ax.plot(times[:], plt_var[:, 0, 0],)
    #plt.ylabel('Volumetric soil water content (m3 m-3)')
    plt.xlabel('Date')
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%B')
    ax.xaxis.set_major_formatter(myFmt)
    if ylab != 'None':
        plt.ylabel(ylab)
    #plt.legend(loc=2)
    #plt.show()
    return fig


if __name__ == '__main__':
    dict_key = sys.argv[1]
    file_dir = sys.argv[2]
    save_dir = sys.argv[3]
    plot_dict = {'sensitivity': (plot_sensitivity, 'sens_sowd_b_smwilt_neff'),
                 'climatology': (plot_climatology, 'clim_79_12')}
    fig = plot_dict[dict_key][0](file_dir, 'croplai', level=0, ylab='Crop LAI')
    fig.savefig(save_dir + plot_dict[dict_key][1]+'_croplai.png', bbox_inches='tight')
    fig = plot_dict[dict_key][0](file_dir, 'cropcanht', level=0, ylab='Crop canopy height (m)')
    fig.savefig(save_dir + plot_dict[dict_key][1]+'_cropcanht.png', bbox_inches='tight')
    fig = plot_dict[dict_key][0](file_dir, 'smcl', level=0, ylab=r'Soil moisture top 10cm (m$^3~$m$^{-3}$)')
    fig.savefig(save_dir + plot_dict[dict_key][1]+'_smcl.png', bbox_inches='tight')
    fig = plot_dict[dict_key][0](file_dir, 't_soil', level=0, ylab='Soil temperature top 10cm (K)')
    fig.savefig(save_dir + plot_dict[dict_key][1]+'_tsoil.png', bbox_inches='tight')
    fig = plot_dict[dict_key][0](file_dir, 'tstar', level=5, ylab='Tile surface temperature')
    fig.savefig(save_dir + plot_dict[dict_key][1]+'_tstar.png', bbox_inches='tight')
    fig = plot_dict[dict_key][0](file_dir, 'fsmc', level=5, ylab='Beta')
    fig.savefig(save_dir + plot_dict[dict_key][1]+'_fsmc.png', bbox_inches='tight')
    if dict_key == 'climatology':
        fig = plot_dict[dict_key][0](file_dir, 'rainfall', ylab='Rainfall rate (kg m$^{-2}$ s$^{-1}$)')
        fig.savefig(save_dir + plot_dict[dict_key][1] + '_rain.png', bbox_inches='tight')