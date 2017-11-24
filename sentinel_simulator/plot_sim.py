import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

import simulator as sim


def plot_class_var(date_lst, var, y_lab=None, line_type='-', axes=None):
    """Plot specified variable.

    :param var: Class attribute variable as list.
    :type var: list
    :param y_lab: Label for Y-axis.
    :type y_lab: str
    :return: Figure.
    :rtype: object
    """
    sns.set_context('poster', font_scale=1.2, rc={'lines.linewidth': 2, 'lines.markersize': 6})
    if axes is not None:
        ax = axes
        ret_val = ax
    else:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 5))
        ret_val = fig, ax
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
    return ret_val


def plot_backscat(lai_coeff, s, omega):
    """Plot backscatter for given coefficients.

    :param lai_coeff: LAI constant used in calculation of extinction coefficient
    :type lai_coeff: float
    :param s: Surface roughness
    :type s: float
    :param omega: Single scattering albedo
    :type omega: float
    :return: Figure
    :rtype: object
    """
    """
    """
    sim_c1 = sim.S1_simulator(lai_coef=lai_coeff, s=s, omega=omega, jules_nc='jules/output/demo/test.3_hourly.nc')
    for x in xrange(3):
        fig = plot_class_var(sim_c1.date_lst,
                             [10*np.log10(sim_c1.SAR_list[s1c].__dict__['stot'][sim_c1.backscatter_keys[x]])
                              for s1c in xrange(len(sim_c1.SAR_list))],
                             y_lab='Backscatter ' + sim_c1.backscatter_keys[x] + ' polarisation (db)',
                             line_type='o')[0]
        # Must also think about canopy height and extinction coefficient!!!
        fig.savefig('jules/output/demo/epstest_' + sim_c1.backscatter_keys[x] + '_laicoeff' + str(np.round(lai_coeff,4))
                    + '.png')
        plt.close('all')
    plt.close('all')
    return 'done'