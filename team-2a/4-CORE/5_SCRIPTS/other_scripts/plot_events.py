#!/usr/bin/python
"""
Description:
Script converts output from the Polaris Compton camera imaging system to a simple csv output
suitable for the reconstruction software CORE.
"""
__author__ = "Dennis Mackin <dsmackin@mdanderson.org>"
__date__ = "Feb. 15, 2016"
__version__ = "$Revision: 0.0.0$"
#------------------------------------------------------------------
# PYTHON IMPORT STATEMENTS
#------------------------------------------------------------------
import sys, os

from math import sin, cos, pi
import pandas
import matplotlib
import matplotlib.pyplot as plt
import cProfile
import re
import csv
import StringIO
import math
import multiprocessing

import numpy

import pyximport; pyximport.install()
pyximport.install(setup_args={'include_dirs': numpy.get_include()})
import utilities


def make_1D_plot(params):
    utilities.plot_1D(*params)
    return

def make_2D_plot(params):
    utilities.plot_2D(*params)
    return

def make_plots(df, output_folder):
    '''
        Creates two lists of plot parameters for 1D and 2D plots. Multiprocessing is used to improve the performance.
    '''

    df['dS'] = numpy.sqrt((df.x1 - df.x2)**2 + (df.y1 - df.y2)**2 + (df.z1 - df.z2)**2)
    df['dT'] = numpy.sqrt((df.x1 - df.x2)**2 + (df.z1 - df.z2)**2)
    df['dY'] = df.y2 - df.y1

    params_list_1D = [(df.E1, 100, "E1", "Counts", 'E1', output_folder),
        (df.dca, 50, "DCA", "Counts", 'DCA(%.1f, %.1f, %.1f)' % (df.dca_x.iloc[1], df.dca_y.iloc[1], df.dca_z.iloc[1]), output_folder),
        (df.theta1, 100, "theta1", "Counts", '1st Scattering Angle', output_folder),
        (df.E2, 100, "E2", "Counts", 'E2', output_folder),
        (df.E, 100, "E", "Counts", 'Energy', output_folder),
        (df.E1 + df.E2, 100, "E", "Counts", 'E1_E2', output_folder),
        (df.x1, 100, "x1", "Counts", 'x1', output_folder),
        (df.x2, 100, "x2", "Counts", 'x2', output_folder),
        (df.y1, 100, "y1", "Counts", 'y1', output_folder),
        (df.dY[(df.dY < 50) & (df.dY > -50)], 100, "dY", "Counts", 'dY', output_folder),
        (df.dT[(df.dY < 50) & (df.dY > -50) & (df.dT < 50)], 100, "dT", "Counts", 'dT', output_folder),
        (df.dS[(df.dY < 50) & (df.dY > -50)], 100, "dS", "Counts", 'dS', output_folder),
        (df.y2, 100, "y2", "Counts", 'y2', output_folder),
        (df.z1, 100, "z1", "Counts", 'z1', output_folder),
        (df.z2, 100, "z2", "Counts", 'z2', output_folder),
        (df.px, 100, "px", "Counts", 'px', output_folder),
        (df.py, 100, "py", "Counts", 'py', output_folder),
        (df.pz, 100, "pz", "Counts", 'pz', output_folder),
        (df.alpha, 100, "alpha", "Counts", 'alpha', output_folder),
    ]

    print df[["theta1", "dca"]].tail()
    print numpy.sum(pandas.isnull(df.theta1)), numpy.sum(pandas.isnull(df.dca))
    params_list_2D = [
        (df.px, "x", 100, df.py, "y", 100, "x y projection", output_folder, False),
        (df.px, "x", 100, df.pz, "z", 100, "x z projection", output_folder, False),

        (df.py, "y", 100, df.pz, "z", 100, "y z projection", output_folder, False),
        (df.theta1, "theta 1", 50, df.dca, "DCA", 50, "DCA vs. Theta1", output_folder, True),
        (df.dca, "DCA", 50, numpy.sin(df.theta1), "sin_theta1", 50, "sin_theta1 vs. DCA", output_folder, False),

        (numpy.sin(df.theta1*math.pi/180.0), "sin_theta1", 100, df.E1, "E1", 100, "E1 vs. sin_theta1", output_folder, False),
        (df.theta1, "theta1", 100, df.E1, "E1", 100, "E1 vs. theta1", output_folder, False),
        (df.theta1, "theta1", 100, df.E1, "E1", 100, "E1 vs. theta1 Log", output_folder, True),
        (df.E1, "E1", 100, df.theta1, "theta1", 100, "theta1 vs. E1", output_folder, False),
        (df.E1, "E1", 100, df.theta1, "theta1", 100, "theta1 vs E1 Log", output_folder, True),
        (df.z1, "z1", 100, df.y1, "y1", 100, "z1_vs_y1", output_folder, False),
        (df.x1, "x1", 100, df.y1, "y1", 100, "x1_vs_y1", output_folder, False),
        (df.x1, "x1", 100, df.z1, "z1", 100, "x1_vs_z1", output_folder, False),
        (df.z2, "x2", 100, df.y2, "y2", 100, "x2_vs_y2", output_folder, False),
        (df.y2, "y2", 100, df.E2, "E2", 100, "E2_vs_y2", output_folder, False),
        (df.E1, "E1", 100, df.E2, "E2", 100, "E2_vs_E1", output_folder, False),
        (df.E1, "E1", 100, df.alpha, "alpha", 100, "alpha_vs_E1", output_folder, False),
        (df.theta1, "theta1", 100, df.alpha, "alpha", 100, "alpha_vs_theta1", output_folder, False),
        (1.33 - df.E1, "1.33minusE1", 100, df.E2, "E2", 100, "E2_vs_minusE1", output_folder, False),
        (df[df.E < 2.0].E, "E0", 50, df[df.E < 2.0].dca, "DCA", 50, "DCA vs. E", output_folder, False),
        (df.dS[(df.dY < 50) & (df.dY > -50)], "dS", 100, df.dca[(df.dY < 50) & (df.dY > -50)], "DCA", 100, "DCA vs. dS", output_folder, True),
        (df.dY[(df.dY < 50) & (df.dY > -50)], "dY", 100, df.dca[(df.dY < 50) & (df.dY > -50)], "DCA", 100, "DCA vs. dY", output_folder, True),
        (df.dT[(df.dY < 50) & (df.dY > -50)], "dT", 100, df.dca[(df.dY < 50) & (df.dY > -50)], "DCA", 100, "DCA vs. dT", output_folder, True), ]

    pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
    return_values = pool.map(make_1D_plot, params_list_1D)
    return_values = pool.map(make_2D_plot, params_list_2D)
    pool.terminate()


    df_slice = df[df.pz > -5]
    df_slice = df_slice[df.pz < 5]
    utilities.plot_2D(df_slice.px, "x", 100, df_slice.py, "y", 100, "x y slice", output_folder, True)


#------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------

def usage():
    print "USAGE: %s [event file] [output folder]" % (sys.argv[0])
    sys.exit(-1)



def main():

    argv = sys.argv
    if not len(argv) == 3:
        usage()

    output_folder = sys.argv[2]
    print "Reading in csv file %s . . ." % sys.argv[1]
    df  = pandas.read_csv(sys.argv[1])
    #df = add_calculated_fields(df)

    print "Plotting histograms . . ."

    make_plots(df, output_folder)



if __name__ == "__main__":

    #cProfile.run("main()")
    main()

