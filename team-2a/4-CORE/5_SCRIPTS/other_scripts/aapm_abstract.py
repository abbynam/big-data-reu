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

def get_calculated_fields(df):

    df["E_known"] = 0.0*df.E2 \
                + (numpy.abs(1.33 - df.E) < numpy.abs(1.17 - df.E) )*1.33 \
                + (numpy.abs(1.33 - df.E) >= numpy.abs(1.17 - df.E) )*1.17
    df["E2_known"] = df["E_known"] - df.E1
    df["dE"] = df["E"] - df["E_known"]

    cos_theta_known = 1.0 - 0.511*(1.0/df["E2_known"] - 1.0/df["E_known"])
    cos_theta_known[numpy.isinf(cos_theta_known)] = numpy.nan
    df["cos_theta_known"] = cos_theta_known

    df["theta_known"] = numpy.arccos(df["cos_theta_known"])
    df["delta_theta"] = df["theta1"]*math.pi/180.0 - df["theta_known"]
    df["cos_delta_theta"] = numpy.cos(df.delta_theta)
    df["sin_delta_theta"] = numpy.sin(numpy.abs(df.delta_theta))

    df["cos_theta1"] = numpy.cos(df.theta1 * math.pi/180.0)
    df["delta_cos_theta"] = df["cos_theta1"] - df["cos_theta_known"]

    return df



def plot_data(df, output_folder):

    df = df[numpy.isnan(df.delta_theta) == False]
    df = df[df.dca < 50]
    utilities.plot_2D(df.delta_theta, "delta_theta", 100, df.dca, "dca", 100, "dTheta_dca", output_folder, False)
    utilities.plot_2D(df[numpy.abs(df.dE) < 0.2 ].delta_theta, "delta_theta", 100, df[numpy.abs(df.dE) < 0.2 ].dE, "dE", 100, "dTheta_dE", output_folder, True)
    utilities.plot_2D(df.cos_delta_theta, "cos_delta_theta", 100, df.dca, "dca", 100, "cosdTheta_dca", output_folder, False)
    utilities.plot_2D(df.sin_delta_theta, "sin_delta_theta", 100, df.dca, "dca", 100, "sindTheta_dca", output_folder, False)
    utilities.plot_2D(df.dE, "dE", 100, df.dca, "dca", 100, "dE_dca", output_folder, True)

    utilities.plot_1D(df.E_known, 100, "E_known", "Counts", 'E_Known', output_folder)
    utilities.plot_1D(df.E2_known, 100, "E2_known", "Counts", 'E2_Known', output_folder)
    utilities.plot_1D(df.dE, 100, "dE", "Counts", 'dE', output_folder)
    utilities.plot_1D(df[numpy.abs(df.cos_theta_known) <= 1.0].cos_theta_known, 100, "cos_theta_known", "Counts", 'cos_theta_known', output_folder)
    utilities.plot_1D(df.theta_known, 100, "theta_known", "Counts", 'theta_known', output_folder)
    utilities.plot_1D(df.delta_theta, 100, "delta_theta", "Counts", 'delta_theta', output_folder)
    utilities.plot_1D(df.cos_delta_theta, 100, "cos_delta_theta", "Counts", 'cos_delta_theta', output_folder)
    utilities.plot_1D(df.sin_delta_theta, 100, "sin_delta_theta", "Counts", 'sin_delta_theta', output_folder)
    utilities.plot_1D(df[numpy.abs(df.cos_theta_known) <= 1.0].delta_cos_theta, 100, "delta_cos_theta", "Counts", 'delta_cos_theta', output_folder)


    # utilities.plot_1D(df.cos_theta_known, 100, "cos_theta_known", "Counts", 'CosThetaKnown', output_folder)
    # utilities.plot_1D(df.theta_117, 100, "Theta117", "Counts", 'Theta117', output_folder)
    # utilities.plot_1D(df.theta_133, 100, "Theta133", "Counts", 'Theta133', output_folder)
    # utilities.plot_1D(df.delta_theta, 100, "DeltaTheta", "Counts", 'DeltaTheta', output_folder)
    # utilities.plot_1D(df.delta_cos_theta, 100, "dCosTheta", "Counts", 'DeltaCosTheta', output_folder)
    #
    # utilities.plot_1D(df.theta_known, 100, "theta_known", "Counts", 'ThetaKnown', output_folder)
    # utilities.plot_1D(df.dca, 100, "DCA", "Counts", 'DCA(%.1f, %.1f, %.1f)' % (df.dca_x[1], df.dca_y[1], df.dca_z[1]), output_folder)
    # utilities.plot_1D(df.theta1, 100, "theta1", "Counts", '1st Scattering Angle', output_folder)
    # utilities.plot_1D(df.E1, 100, "E1", "Counts", 'E1', output_folder)
    # utilities.plot_1D(df.E2, 100, "E2", "Counts", 'E2', output_folder)
    # utilities.plot_1D(df.E3, 100, "E3", "Counts", 'E3', output_folder)
    # utilities.plot_1D(df.E[df.E < 3.0], 100, "E", "Counts", 'Energy', output_folder)
    # utilities.plot_1D(df.y3, 100, "y3", "Counts", 'y3', output_folder)
    # utilities.plot_2D(df.theta1, "theta 1", 100, df.dca, "DCA", 100, "DCA vs. Theta1", output_folder, True)
    #
    # utilities.plot_1D(df.cos_theta1, 100, "cos_theta1", "Counts", 'cos_theta1', output_folder)
    # utilities.plot_2D(df.dca, "DCA", 100, df.cos_theta1, "cos_theta1", 100, "cos_theta1 vs. DCA", output_folder, False)
    # utilities.plot_2D(df.cos_theta1, "cos_theta1", 100, df.E1, "E1", 100, "E1 vs. cos_theta1", output_folder, False)
    # utilities.plot_2D(df[df.E1 < 2.0].theta1, "theta1", 100, df[df.E1 < 2.0].E1, "E1", 100, "E1 vs. theta1", output_folder, False)
    # utilities.plot_2D(df[df.E1 < 2.0].theta1, "theta1", 100, df[df.E1 < 2.0].E1, "E1", 100, "E1 vs. theta1 Log", output_folder, True)
    #
    # utilities.plot_2D(df[df.delta_cos_theta < 0.3].theta1, "theta1", 100, \
    #                   df[df.delta_cos_theta < 0.3].E1, "E1", 100, "dCosTheta_LT_0.3 ", output_folder, False)
    # utilities.plot_2D(df[df.delta_cos_theta < 0.1].theta1, "theta1", 100, df[df.delta_cos_theta < 0.1].E1, "E1", 100, "dCosTheta_LT_0.3 Log", output_folder, True)
    # utilities.plot_2D(df[df.dE < 0.05].theta1, "theta1", 100, df[df.dE < 0.05].E1, "E1", 100, "dE_LT_0.05_Log", output_folder, True)
    #
    # utilities.plot_1D(df[df.dE < 0.05].dca, 100, "DCA", "CountsDE", 'DCA(%.1f, %.1f, %.1f): dE_LT_0.05' % (df.dca_x[1], df.dca_y[1], df.dca_z[1]), output_folder)
    # utilities.plot_1D(df[df.delta_cos_theta < 0.3].dca, 100, "DCA", "CountsDcose", 'DCA(%.1f, %.1f, %.1f): dCosTheta_LT_0.3' % (df.dca_x[1], df.dca_y[1], df.dca_z[1]), output_folder)
    # utilities.plot_2D(df.z1, "z1", 100, df.y1, "y1", 100, "z1_vs_y1", output_folder, False)
    # utilities.plot_2D(df.x1, "x1", 100, df.y1, "y1", 100, "x1_vs_y1", output_folder, False)
    # utilities.plot_2D(df.x1, "x1", 100, df.z1, "z1", 100, "x1_vs_z1", output_folder, False)
    # utilities.plot_2D(df.z2, "x2", 100, df.y2, "y2", 100, "z2_vs_y2", output_folder, False)
    # utilities.plot_2D(df.y2, "y2", 100, df.E2, "E2", 100, "E2_vs_y2", output_folder, False)
    #
    # utilities.plot_2D(df[df.E < 2.0].E, "E0", 100, df[df.E < 2.0].dca, "DCA", 100, "DCA vs. E", output_folder, False)
    # utilities.plot_2D(df[df.E < 2.0].dE, "dE", 100, df[df.E < 2.0].dca, "DCA", 100, "DCA vs. deltaE", output_folder, False)
    # utilities.plot_2D(df[df.E < 2.0].delta_cos_theta, "delta_cos_theta", 100, df[df.E < 2.0].dE, "dE", 100, "deltaE vs dCosTheta", output_folder, True)

def save_filtered_events(df, output_folder):

    dE_cut = 0.2
    dCosTheta_cut = 0.3
    cosDTheta_cut = 0.2
    sinDTheta_cut = 0.2
    dTheta_cut = 0.2

    df = df[df.dca < 100.0]
    df = df[numpy.isnan(df.delta_theta) == False]

    colnames = df.columns.values
    theta_M_events = df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events.to_csv('%s/theta_m.csv' % output_folder, index=False, header=False)

    theta_C_events = df[['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_events.to_csv('%s/theta_c.csv' % output_folder, index=False, header=False)

    df["abs_dE"] = numpy.abs(df.dE)
    df = df.sort_values(by=['abs_dE'], ascending=[True])

    theta_M_events_cut = df[df.abs_dE  < dE_cut][['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_M_events_cut.to_csv('%s/theta_c_cutDE.csv' % output_folder, index=False, header=False)

    theta_M_events_cut = df[df.abs_dE < dE_cut][['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events_cut.to_csv('%s/theta_m_cutDE.csv' % output_folder, index=False, header=False)

    df = df.sort_values(by=['sin_delta_theta'], ascending=[True])
    theta_M_events_cut_sinDtheta = df[df.sin_delta_theta < sinDTheta_cut][['E1','x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events_cut_sinDtheta.to_csv('%s/theta_m_cutSinDelta.csv' % output_folder, index=False, header=False)

    theta_C_events_cut_sinDtheta = df[df.sin_delta_theta < sinDTheta_cut][['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_events_cut_sinDtheta.to_csv('%s/theta_c_cutSinDelta.csv' % output_folder, index=False, header=False)



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
    df = read_csv_events(sys.argv[1])
    df = get_calculated_fields(df)
    save_filtered_events(df, output_folder)

    print "Plotting histograms . . ."
    plot_data(df, output_folder)



if __name__ == "__main__":

    #cProfile.run("main()")
    main()

