#!/usr/bin/python
"""
Description:
Produces a CDF plot comparing the DCA values of sets of events. This script was originally created
to evaluate dtheta as an event selection method compared to energy windowing or multiple energy windowing.

In this case we use
"""
__author__ = "Dennis Mackin <dsmackin@mdanderson.org>"
__date__ = "Feb. 15, 2016"
__version__ = "$Revision: 0.0.0$"
#------------------------------------------------------------------
# PYTHON IMPORT STATEMENTS
#------------------------------------------------------------------
import sys
import numpy
import pandas
import time
import cProfile
import matplotlib.pyplot as plt
from itertools import izip

# import pyximport; pyximport.install(reload_support=True)
import utilities
import energy_matcher


def read_csv_events(file_path):
    df = pandas.read_csv(file_path, names=['E1', 'x1', 'y1', 'z1', 'E2', 'x2', 'y2', 'z2'])

    return df


def make_DCA_CDF_plot(dataframes, moniker, title, output_folder, z=0):
    plt.clf()
    ax = plt.subplot(111)

    labels = [r'$\Delta$E', r'$\Delta$$\Theta$', 'window', 'none']
    for i, df in enumerate(dataframes):
        # df =df[df.dca < 100]
        x, cdf_vals = utilities.get_CDF(df.dca.as_matrix(), 100)

        plt.plot(x, cdf_vals, label="%s" % (labels[i]), lw=1)
        print "%s Min x:" % labels[i], min(x)

    plt.title(title, fontsize=11)
    legend = ax.legend(loc='lower right', title=r'Filter', fontsize=8)
    ax.grid(True,linestyle='-',which='both',color='0.750')

    plt.setp(ax.get_xticklabels(), fontsize=9)
    plt.setp(ax.get_yticklabels(), fontsize=9)
    plt.setp(legend.get_title(), fontsize=9)
    plt.xlabel("DCA [mm]", fontsize=10)
    # ax.set_xlim([75.0, 200.0])

    plt.ylabel("Selected Event CDF", fontsize=10)
    plt.gcf().set_size_inches(5,3.3)
    plt.gcf().subplots_adjust(bottom=0.15)

    plot_name = "%s/DCA_CDF_%s_%.1f.png" % (output_folder, moniker, z)
    plt.savefig(plot_name)

    return True



#------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------


def usage():
    print "USAGE: %s [event file] [comma seperated energies in MeV] [output folder]" % (sys.argv[0])
    print "\nExample:\n %s /home/dsmackin/public_html/Co60" % (sys.argv[0])
    sys.exit(-1)

def main():
    argv = sys.argv
    if not len(argv) == 2:
        usage()

    output_folder = sys.argv[1]

    # raw = read_csv_events("/home/dsmackin/projects/CORE/out/Co60c_raw/events.dat")
    # window = read_csv_events("/home/dsmackin/projects/CORE/out/Co60c_window/events.dat")
    # theta = read_csv_events("/home/dsmackin/projects/CORE/out/Co60c_dthetam/events.dat")
    # energy = read_csv_events("/home/dsmackin/projects/CORE/out/Co60c_dE/events.dat")
    #
    # make_DCA_CDF_plot([energy, theta, window, raw], "Co60centered", "Co-60 Source at (0 cm, 0 cm, 0 cm)", output_folder, max_dca=100)
    #
    # raw = read_csv_events("/home/dsmackin/projects/CORE/out/Co60s_raw/events.dat")
    # window = read_csv_events("/home/dsmackin/projects/CORE/out/Co60s_window/events.dat")
    # theta = read_csv_events("/home/dsmackin/projects/CORE/out/Co60s_dthetam/events.dat")
    # energy = read_csv_events("/home/dsmackin/projects/CORE/out/Co60s_dE/events.dat")
    #
    # make_DCA_CDF_plot([energy, theta, window, raw], "Co60shifted", "Co-60 Source at (0 cm, 0 cm,-60 cm)", output_folder, max_dca=100)
    #
    # raw = read_csv_events("/home/dsmackin/projects/CORE/out/Cs137_raw/events.dat")
    # window = read_csv_events("/home/dsmackin/projects/CORE/out/Cs137_window/events.dat")
    # theta = read_csv_events("/home/dsmackin/projects/CORE/out/Cs137_dthetam/events.dat")
    # energy = read_csv_events("/home/dsmackin/projects/CORE/out/Cs137_dE/events.dat")

    # window = read_csv_events("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-m__MAXIMUM_NUMBER_CONES-1000/events.dat")
    # raw = read_csv_events("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-c_dca__MAXIMUM_NUMBER_CONES-1000/events.dat")
    # energy = read_csv_events("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-c_dE__MAXIMUM_NUMBER_CONES-1000/events.dat")
    # theta = read_csv_events("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-c_dTheta__MAXIMUM_NUMBER_CONES-1000/events.dat")

    #Co60 Source
    fraction_of_events = 0.1 #keep the best 10% of events
    window_range = [0.7*1.17, 1.03*1.33]
    points = [numpy.array([0.0, 0.0, 0.0]), numpy.array([0.0, 0.0, 32.5]), numpy.array([0.0, 0.0, -32.5])]
    folders = ["/home/dsmackin/public_html/dtheta_centered", "/home/dsmackin/public_html/dtheta_mod2", "/home/dsmackin/public_html/dtheta_mod3"]
    for dca_point, folder in izip(points, folders):

        # dca_point = numpy.array([0.0, 0.0, 32.5])
        # folder = "/home/dsmackin/public_html/dtheta_mod2"
        basic = read_csv_events("%s/theta_m.csv" % folder)
        t0 = time.clock()
        for i in range(20):
            basic['dca'] = energy_matcher.get_dca_for_events(basic.as_matrix(), dca_point)
        print "DCA call took %s secs for 20 calls. . ." % (time.clock() - t0)
        mask = basic['dca'] > 0.0
        # basic.loc[mask, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]]
        # basic.loc[mask, ['dca']] = basic.loc[mask, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]]
        basic = basic[mask]

        window = basic[basic.E1 + basic.E2 <  window_range[1]]
        window = window[window.E1 + window.E2 >  window_range[0]]
        window['dca'] = energy_matcher.get_dca_for_events(window.as_matrix(), dca_point)

        num_events = int(fraction_of_events * len(basic))

        energy = read_csv_events("%s/theta_m_dE.csv" % folder)
        energy = energy[0:num_events]
        energy['dca'] = energy_matcher.get_dca_for_events(energy.as_matrix(), dca_point)

        theta = read_csv_events("%s/theta_m_dTheta.csv" % folder)
        theta = theta[0:num_events]
        theta['dca'] = energy_matcher.get_dca_for_events(theta.as_matrix(), dca_point)

        make_DCA_CDF_plot([energy, theta, window, basic], "Co60", "Co-60 Source at (%.0f, %.0f, %.0f) mm" % (dca_point[0], dca_point[1], dca_point[2]), output_folder, dca_point[2])

    #
    # #Cs137 Source
    # window_range = [0.7*0.662, 1.03*0.662]
    # basic = read_csv_events("/home/dsmackin/public_html/Cs/theta_m.csv")
    # basic['dca'] = energy_matcher.get_dca_for_events(basic.as_matrix(), dca_point)
    #
    # window = basic[basic.E1 + basic.E2 <  window_range[1]]
    # window = window[window.E1 + window.E2 >  window_range[0]]
    # window['dca'] = energy_matcher.get_dca_for_events(window.as_matrix(), dca_point)
    #
    # num_events = int(fraction_of_events * len(basic))
    #
    # energy = read_csv_events("/home/dsmackin/public_html/Cs/theta_m_dE.csv")
    # energy = energy[0:num_events]
    # energy['dca'] = energy_matcher.get_dca_for_events(energy.as_matrix(), dca_point)
    #
    # theta = read_csv_events("/home/dsmackin/public_html/Cs/theta_m_dTheta.csv")
    # theta = theta[0:num_events]
    # theta['dca'] = energy_matcher.get_dca_for_events(theta.as_matrix(), dca_point)
    #
    # make_DCA_CDF_plot([energy, theta, window, basic], "Cs137", "Cs-137 Source at (0, 0, 0) mm", output_folder, max_dca=100)
    #
    #

if __name__ == "__main__":

    cProfile.run("main()")
    # main()

