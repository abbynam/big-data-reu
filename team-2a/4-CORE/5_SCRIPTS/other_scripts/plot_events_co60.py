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
    df = pandas.read_csv(file_path)

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



def plot_dca(df, output_folder):

    # df = df[numpy.isnan(df.delta_theta) == False]

    # utilities.plot_1D(df.delta_theta, 100, "delta_theta", "Counts", 'delta_theta', output_folder)
    # utilities.plot_2D(df.theta1, "theta1", 100, df.delta_theta, "delta_theta", 100, "dTheta vs. Theta1", output_folder, True)
    # utilities.plot_2D(df.delta_theta, "dTheta", 100, df.dca, "DCA", 100, "DCA vs. dTheta", output_folder, True)
    # utilities.plot_2D(df.dE, "dE", 100, df.dca, "DCA", 100, "DCA vs. dE", output_folder, True)

    print len(df.dca)
    print len(df.dca_x)
    print len(df.dca_y)
    print len(df.dca_z)
    print df.dca_x.iloc[1]
    print df.dca_y.iloc[1]
    print df.dca_z.iloc[1]

    utilities.plot_2D(df.theta1, "theta1", 100, df.dca, "dca", 100, "dca vs. theta1 all", output_folder, False)

    utilities.plot_2D(df.theta1, "theta", 100, df.theta1 - df.theta1_c, "dtheta", 100, "dtheta", output_folder, False)
    utilities.plot_2D(df.theta1, "theta", 100, df.theta1_c, "thetac", 100, "thetavtheta", output_folder, False)
    utilities.plot_2D(numpy.sin((df.theta1 - df.theta1_c)*math.pi/180.0), "sin_dtheta", 100, df.dca, "dca", 100, "dca vs. sin(dtheta)", output_folder, False)
    utilities.plot_2D(numpy.sin((df.theta1 - df.theta1_c)*math.pi/180.0), "sin_dtheta", 100, df.E1, "E1", 100, "E1 vs. sin(dtheta)", output_folder, False)
    utilities.plot_2D(df.E1 - df.E1_c, "dE1", 100, df.dca, "dca", 100, "dca vs. dE1", output_folder, False)

    utilities.plot_1D(df.theta1_c, 100, "theta1c", "Counts", 'theta1_c', output_folder)
    utilities.plot_1D(df.theta1,   100, "theta1", "Counts", 'theta1', output_folder)
    utilities.plot_1D(df.theta1 - df.theta1_c, 100, "dtheta1", "Counts", 'dtheta1', output_folder)
    utilities.plot_1D(df.E1, 100, "E1", "Counts", 'E1', output_folder)

    utilities.plot_1D(df.E1, 100, "E1", "Counts", 'E1', output_folder)
    utilities.plot_1D(df.E1 - df.E1_c, 100, "dE1", "Counts", 'dE1', output_folder)

    utilities.plot_1D(df[df.E_k == 1.33].E1 - df[df.E_k == 1.33].E1_c, 100, "dE1", "Counts", 'dE1(1.33)', output_folder)
    utilities.plot_1D(df[df.E_k == 1.33].E - df[df.E_k == 1.33].E_k, 100, "dEkEm", "Counts", 'd(Ek,Em;1.33)', output_folder)
    utilities.plot_1D(df[df.E_k == 1.33].E, 100, "Em(Ek=1.33)", "Counts", 'E(1.33)', output_folder)

    utilities.plot_2D(df.E1, "E1", 100, df.E1_c, "E1_c", 100, "E1 vs E1c", output_folder, True)
    utilities.plot_2D(df.E1, "E1", 100, df.E1 - df.E1_c, "dE1", 100, "E1 vs dE1", output_folder, False)



    utilities.plot_1D(df.E1_c, 100, "E1c", "Counts", 'E1c', output_folder)
    utilities.plot_1D(df.dca, 100, "DCA", "Counts", 'DCA(%.1f, %.1f, %.1f)' % (df.dca_x.iloc[1], df.dca_y.iloc[1], df.dca_z.iloc[1]), output_folder)

    # utilities.plot_1D(df.dca_c, 100, "DCA_c", "Counts", 'DCA(%.1f, %.1f, %.1f)' % (df.dca_x.iloc[1], df.dca_y.iloc[1], df.dca_z.iloc[1]), output_folder)

    utilities.plot_1D(df.E2, 100, "E2", "Counts", 'E2', output_folder)
    df = df[df.E_k == 1.33]
    return


    # utilities.plot_1D(df.E3, 100, "E3", "Counts", 'E3', output_folder)
    utilities.plot_1D(df.E, 100, "E", "Counts", 'Energy', output_folder)
    utilities.plot_1D(df.x1, 100, "x1", "Counts", 'x1', output_folder)
    utilities.plot_1D(df.x2, 100, "x2", "Counts", 'x2', output_folder)
    utilities.plot_1D(df.y1, 100, "y1", "Counts", 'y1', output_folder)
    utilities.plot_1D(df.y2, 100, "y2", "Counts", 'y2', output_folder)
    utilities.plot_1D(df.z1, 100, "z1", "Counts", 'z1', output_folder)
    utilities.plot_1D(df.z2, 100, "z2", "Counts", 'z2', output_folder)
    utilities.plot_1D(df.px, 100, "px", "Counts", 'px', output_folder)
    utilities.plot_1D(df.py, 100, "py", "Counts", 'py', output_folder)
    utilities.plot_1D(df.pz, 100, "pz", "Counts", 'pz', output_folder)
    utilities.plot_1D(df.pz, 100, "alpha", "Counts", 'alpha', output_folder)
    utilities.plot_2D(df.theta1, "theta 1", 100, df.dca, "DCA", 100, "DCA vs. Theta1", output_folder, True)

    utilities.plot_1D(numpy.cos(df.theta1*math.pi/180.0), 100, "cos_theta1", "Counts", 'cos_theta1', output_folder)
    utilities.plot_2D(df.dca, "DCA", 100, numpy.cos(df.theta1), "cos_theta1", 100, "cos_theta1 vs. DCA", output_folder, False)
    utilities.plot_2D(numpy.cos(df.theta1*math.pi/180.0), "cos_theta1", 100, df.E1, "E1", 100, "E1 vs. cos_theta1", output_folder, False)
    utilities.plot_2D(df.theta1, "theta1", 100, df.E1, "E1", 100, "E1 vs. theta1", output_folder, False)
    utilities.plot_2D(df.theta1, "theta1", 100, df.E1, "E1", 100, "E1 vs. theta1 Log", output_folder, True)
    utilities.plot_2D(df.E1, "E1", 100, df.theta1, "theta1", 100, "theta1 vs. E1", output_folder, False)
    utilities.plot_2D(df.E1, "E1", 100, df.theta1, "theta1", 100, "theta1 vs E1 Log", output_folder, True)


    # utilities.plot_2D(df[df.delta_cos_theta < 0.3].theta1, "theta1", 100, \
    #                   df[df.delta_cos_theta < 0.3].E1, "E1", 100, "dCosTheta_LT_0.3 ", output_folder, False)
    # utilities.plot_2D(df[df.delta_cos_theta < 0.1].theta1, "theta1", 100, df[df.delta_cos_theta < 0.1].E1, "E1", 100, "dCosTheta_LT_0.3 Log", output_folder, True)
    # utilities.plot_2D(df[df.dE < 0.05].theta1, "theta1", 100, df[df.dE < 0.05].E1, "E1", 100, "dE_LT_0.05_Log", output_folder, True)

    # utilities.plot_1D(df[df.dE < 0.05].dca, 100, "DCA", "CountsDE", 'DCA(%.1f, %.1f, %.1f): dE_LT_0.05' % (df.dca_x[1], df.dca_y[1], df.dca_z[1]), output_folder)
    # utilities.plot_1D(df[df.delta_cos_theta < 0.3].dca, 100, "DCA", "CountsDcose", 'DCA(%.1f, %.1f, %.1f): dCosTheta_LT_0.3' % (df.dca_x[1], df.dca_y[1], df.dca_z[1]), output_folder)
    utilities.plot_2D(df.z1, "z1", 100, df.y1, "y1", 100, "z1_vs_y1", output_folder, False)
    utilities.plot_2D(df.x1, "x1", 100, df.y1, "y1", 100, "x1_vs_y1", output_folder, False)
    utilities.plot_2D(df.x1, "x1", 100, df.z1, "z1", 100, "x1_vs_z1", output_folder, False)
    utilities.plot_2D(df.z2, "x2", 100, df.y2, "y2", 100, "x2_vs_y2", output_folder, False)
    utilities.plot_2D(df.y2, "y2", 100, df.E2, "E2", 100, "E2_vs_y2", output_folder, False)

    utilities.plot_2D(df[df.E < 2.0].E, "E0", 100, df[df.E < 2.0].dca, "DCA", 100, "DCA vs. E", output_folder, False)
    # utilities.plot_2D(df[df.E < 2.0].dE, "dE", 100, df[df.E < 2.0].dca, "DCA", 100, "DCA vs. deltaE", output_folder, False)
    # utilities.plot_2D(df[df.E < 2.0].delta_cos_theta, "delta_cos_theta", 100, df[df.E < 2.0].dE, "dE", 100, "deltaE vs dCosTheta", output_folder, True)


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

    plot_dca(df, output_folder)

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

