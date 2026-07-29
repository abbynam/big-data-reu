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
import numpy
from math import sin, cos, pi
import pandas
import matplotlib
import matplotlib.pyplot as plt
import cProfile
import re
import csv
import StringIO
import math

import utilities


def read_csv_events(file_path):
    df = pandas.read_csv(file_path, names=['E1', 'x1', 'y1', 'z1', 'E2', 'x2', 'y2', 'z2', 'E3', 'x3', 'y3', 'z3'])

    return df
#
# def add_calculated_fields(df):
#     colnames = df.columns.values
#     df["cos_theta1"] = numpy.cos(df.theta1 * math.pi/180.0)
#     d117 = numpy.abs(1.17 - df.E)
#     d133 = numpy.abs(1.33 - df.E)
#
#     df["E_known"] = 0.0*df.E2 \
#                 + (numpy.abs(1.33 - df.E) < numpy.abs(1.17 - df.E) )*1.33 \
#                 + (numpy.abs(1.33 - df.E) >= numpy.abs(1.17 - df.E) )*1.17
#     known_E2 = df["E_known"] - df.E1
#     print min(known_E2)
#     cos_theta_known = 1.0 - 0.511*(1.0/(df["E_known"] - df.E1) - 1.0/df["E_known"])
#     cos_theta_known[numpy.isinf(cos_theta_known)] = numpy.nan
#     df["cos_theta_known"] = cos_theta_known
#     bad_E1 = df[numpy.isnan(df["cos_theta_known"])].E1
#     print bad_E1, numpy.min(df["E_known"])
#     print sum(df["cos_theta_known"])
#     df["theta_known"] = numpy.arccos(df["cos_theta_known"])
#
#     df["dE"] = numpy.minimum(d117, d133)
#     df["cos_theta_117"] = 1 - 0.511*(1.0/(1.17 - df.E1) - 1.0/1.17)
#     df["cos_theta_133"] = 1 - 0.511*(1.0/(1.33 - df.E1) - 1.0/1.33)
#
#     df["theta_117"] = numpy.arccos(df["cos_theta_117"])
#     df["theta_133"] = numpy.arccos(df["cos_theta_133"])
#     df["delta_cos_theta_117"] = df["cos_theta1"] - df["cos_theta_117"]
#     df["delta_cos_theta_133"] = df["cos_theta1"] - df["cos_theta_133"]
#
#     df["delta_cos_theta"] = numpy.minimum(numpy.abs(df["delta_cos_theta_117"]), numpy.abs(df["delta_cos_theta_133"]))
#     df["delta_theta"] = numpy.minimum(numpy.abs(df["theta_117"] - df["theta1"]*(math.pi/180.0)), numpy.abs(df["theta_133"] - df["theta1"]*(math.pi/180.0)))
#
#     return df



def make_plots(df, output_folder):

    utilities.plot_1D(df.E1, 100, "E1", "Counts", 'E1', output_folder)
    utilities.plot_1D(df.E2, 100, "E2", "Counts", 'E2', output_folder)
    # utilities.plot_1D(df.E3, 100, "E3", "Counts", 'E3', output_folder)
    utilities.plot_1D(df.E1 + df.E2, 100, "E", "Counts", 'Energy', output_folder)
    utilities.plot_1D(df.x1, 100, "x1", "Counts", 'x1', output_folder)
    utilities.plot_1D(df.x2, 100, "x2", "Counts", 'x2', output_folder)
    utilities.plot_1D(df.y1, 100, "y1", "Counts", 'y1', output_folder)
    utilities.plot_1D(df.y2, 100, "y2", "Counts", 'y2', output_folder)
    utilities.plot_1D(df.z1, 100, "z1", "Counts", 'z1', output_folder)

    utilities.plot_2D(df.z1, "z1", 100, df.y1, "y1", 100, "z1_vs_y1", output_folder, False)
    utilities.plot_2D(df.x1, "x1", 100, df.y1, "y1", 100, "x1_vs_y1", output_folder, False)
    utilities.plot_2D(df.x1, "x1", 100, df.z1, "z1", 100, "x1_vs_z1", output_folder, False)
    utilities.plot_2D(df.z2, "x2", 100, df.y2, "y2", 100, "x2_vs_y2", output_folder, False)
    utilities.plot_2D(df.y2, "y2", 100, df.E2, "E2", 100, "E2_vs_y2", output_folder, False)

#------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------

def usage():
    print "USAGE: %s [event file] [output folder]" % (sys.argv[1])
    sys.exit(-1)

def main():

    argv = sys.argv
    if not len(argv) == 3:
        usage()

    output_folder = sys.argv[2]
    print "Reading in csv file %s . . ." % sys.argv[1]
    df  = read_csv_events(sys.argv[1])
    #df = add_calculated_fields(df)

    print "Plotting histograms . . ."

    make_plots(df, output_folder)

    # utilities.plot_1D(df_allevents.E1[df_allevents.E1 < 3.0], 100, 'E1', "Counts", "First Scatter Energy", output_folder)
    #
    # total_energy = (df_allevents.E1 + df_allevents.E2)
    # utilities.plot_1D(total_energy[total_energy < 3.0], 100, 'E1 + E2', "Counts", "Total Energy", output_folder)
    #
    # utilities.plot_2D(df_allevents.z2, "z", df_allevents.x2, "x", "Scatter2", output_folder, True)
    #
    # utilities.plot_1D(df_allevents.x1, 200, 'x1', "Counts", "All Events: x1", output_folder)


if __name__ == "__main__":

    #cProfile.run("main()")
    main()

