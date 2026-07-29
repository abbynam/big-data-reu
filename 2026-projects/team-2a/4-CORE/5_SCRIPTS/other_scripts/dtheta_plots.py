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


def get_E_known(E, known_energies):

    #if there's only one energy, select it
    if len(known_energies) == 1:
        return E*0.0 + known_energies[0]

    known_energies = numpy.unique(known_energies)
    compton_maxes = known_energies - (0.511 * known_energies / (0.511 + 2.0 * known_energies))
    compton_maxes_array = numpy.repeat(compton_maxes, len(E))
    compton_maxes_array = compton_maxes_array.reshape((len(E), len(compton_maxes)), order='F')


    e1_known = e_known * e_known * (1.0 - cos_theta) / (0.511 + e_known * (1.0 - cos_theta))

    ke_array = numpy.repeat(known_energies, len(E))
    ke_array = ke_array.reshape((len(E), len(known_energies)), order='F')
    assert(len(ke_array) == len(E))

    E_array = numpy.tile(E, (1, len(known_energies)))[0]
    E_array = E_array.reshape((len(E), len(known_energies)), order='F')
    assert(len(E_array) == len(E))

    diff_array = numpy.abs((ke_array - E_array))

    min_vals = numpy.amin(diff_array, axis=1, keepdims=True)
    min_truths = (diff_array == min_vals)
    # E_known = E_array[min_truths]
    E_known = numpy.multiply(ke_array, min_truths)
    E_known = numpy.sum(E_known,axis=1)

    assert(len(E_known) == len(E))

    return E_known


def get_calculated_fields(df, known_energies):

    df["E"] = df.E1 + df.E2
    df["cos_theta1"] = 1.0 - 0.511*(1.0/df.E2 - 1.0/df.E)
    df["theta1"] = numpy.arccos(df["cos_theta1"])

    known_energies_arr = get_E_known(df.E, known_energies)
    df_len = len(df.E)
    ke_len = len(known_energies_arr)
    assert(len(known_energies_arr) == len(df))
    df["E_known"] = known_energies_arr

    df["E2_known"] = df["E_known"] - df.E1
    df["dE"] = df["E"] - df["E_known"]
    df["dE_scaled"] = df["dE"]/(df["E_known"] + 1E-10)

    cos_theta_known = 1.0 - 0.511*(1.0/df["E2_known"] - 1.0/df["E_known"])
    cos_theta_known[numpy.isinf(cos_theta_known)] = numpy.nan
    df["cos_theta_known"] = cos_theta_known

    df["theta_known"] = numpy.arccos(df["cos_theta_known"])
    #df["delta_theta"] = df["theta1"]*math.pi/180.0 - df["theta_known"]
    df["delta_theta"] = df["theta1"] - df["theta_known"]
    df["cos_delta_theta"] = numpy.cos(df.delta_theta)
    df["sin_delta_theta"] = numpy.sin(numpy.abs(df.delta_theta))

    df["delta_cos_theta"] = df["cos_theta1"] - df["cos_theta_known"]
    df["scatter_distance"] = numpy.sqrt(numpy.square(df.x1 - df.x2) + numpy.square(df.y1 - df.y2) + numpy.square(df.z1 - df.z2) )

    return df


def plot_data(df, output_folder):

    df = df[numpy.isnan(df.delta_theta) == False]
    # df = df[df.y1 > 180]
    # df = df[df.y2 > 180]
    # df = df[df.y1 < df.y2]
    df_bad = df[df.dca > 100]
    df_good = df[df.dca < 20]

    utilities.plot_1D(df.dca, 100, "DCA", "Counts", 'DCA(%.1f, %.1f, %.1f)' % (df.dca_x.iloc[1], df.dca_y.iloc[1], df.dca_z.iloc[1]), output_folder)
    for i,d in enumerate([df_good, df_bad]):
        dy = d.y1 - d.y2
        if len(dy) < 5: continue
        utilities.plot_1D(dy.as_matrix(), 100, "dy1y2", "Counts_%d" % i, "dy1y2_%d" % i, output_folder)
        utilities.plot_1D(d.theta1*180/pi, 100, "theta", "Counts_%d" % i, "theta1_%d" % i, output_folder)
        utilities.plot_2D(d.y1 - d.y2, "dY1Y2", 100, d.dca, "DCA_%d" % i, 100, "dY1Y2 vs. DCA", output_folder, False)
        utilities.plot_2D(d.y1 - d.y2, "dY1Y2", 100, d.theta1*180/pi, "theta_%d" % i, 100, "dY1Y2 vs. theta", output_folder, False)
        utilities.plot_2D(d.y1 - d.y2, "dY1Y2", 100, d.sin_delta_theta, "sindtheta_%d" % i, 100, "dY1Y2 vs. delta_theta", output_folder, False)
        utilities.plot_2D(d.y1 - d.y2, "dY1Y2", 100, d.E1 - d.E2, "dE1E2_%s" % i, 100, "dY1Y2 vs. dE1E2", output_folder, False)
        utilities.plot_2D(d.E1 - d.E2, "dE1E2", 100, d.dca, "DCA_%d" % i, 100, "dE1E2 vs. DCA", output_folder, False)
        utilities.plot_2D(d.E1 - d.E2, "dE1E2", 100, d.theta1*180/pi, "theta_%d" % i, 100, "dE1E2 vs. theta", output_folder, False)
        utilities.plot_2D(d.E1 - d.E2, "dE1E2", 100, d.sin_delta_theta, "sindtheta_%d" % i, 100, "dE1E2 vs. delta_theta", output_folder, False)
        utilities.plot_1D(d.E1 - d.E2, 100, "dE1E2", "Counts_%d" % i, "dE1E2_%d" % i, output_folder)

    utilities.plot_2D(df[numpy.abs(df.dE) < 0.2 ].delta_theta, "delta_theta", 100, df[numpy.abs(df.dE) < 0.2 ].dE, "dE", 100, "dTheta_dE", output_folder, True)

    utilities.plot_2D(df.sin_delta_theta, "sin_delta_theta", 100, df.dca, "dca", 100, "sindTheta_dca", output_folder, True)
    utilities.plot_2D(df.sin_delta_theta, "sin_delta_theta", 100, df.dca, "dca", 100, "sindTheta_dca", output_folder, False)
    utilities.plot_2D(df.dE, "dE", 100, df.dca, "dca", 100, "dE_dca", output_folder, False)

    utilities.plot_2D(df.z1, "z1", 100, df.y1, "y1", 100, "z1_vs_y1", output_folder, False)
    utilities.plot_2D(df.x1, "x1", 100, df.y1, "y1", 100, "x1_vs_y1", output_folder, False)
    utilities.plot_2D(df.x1, "x1", 100, df.z1, "z1", 100, "x1_vs_z1", output_folder, False)
    utilities.plot_2D(df.z1, "z1", 100, df.y1, "y1", 100, "z1_vs_y1", output_folder, True)
    utilities.plot_2D(df.x1, "x1", 100, df.y1, "y1", 100, "x1_vs_y1", output_folder, True)
    utilities.plot_2D(df.x1, "x1", 100, df.z1, "z1", 100, "x1_vs_z1", output_folder, True)


    return
    utilities.plot_2D(df.y1 - df.y2, "dY", 100, df.dca, "DCA", 100, "DCA vs. dY", output_folder, True)
    utilities.plot_2D(df.y1 - df.y2, "dY", 100, df.dca, "DCA", 100, "DCA vs. dY", output_folder, False)
    utilities.plot_2D(df.y1 - df.y2, "dY", 100, df.E1 - df.E2, "dE1E2", 100, "dE1E2 vs. dY", output_folder, True)
    utilities.plot_2D(df.E1 - df.E2, "dE1E2", 100, df.dca, "DCA", 100, "dE1E2 vs. DCA", output_folder, True)
    utilities.plot_2D(df.delta_theta, "delta_theta", 100, df.dca, "dca", 100, "dTheta_dca", output_folder, True)

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
    # utilities.plot_2D(df.z2, "x2", 100, df.y2, "y2", 100, "z2_vs_y2", output_folder, False)
    # utilities.plot_2D(df.y2, "y2", 100, df.E2, "E2", 100, "E2_vs_y2", output_folder, False)
    #
    # utilities.plot_2D(df[df.E < 2.0].E, "E0", 100, df[df.E < 2.0].dca, "DCA", 100, "DCA vs. E", output_folder, False)
    # utilities.plot_2D(df[df.E < 2.0].dE, "dE", 100, df[df.E < 2.0].dca, "DCA", 100, "DCA vs. deltaE", output_folder, False)
    # utilities.plot_2D(df[df.E < 2.0].delta_cos_theta, "delta_cos_theta", 100, df[df.E < 2.0].dE, "dE", 100, "deltaE vs dCosTheta", output_folder, True)

def save_filtered_events(df, output_folder):

    dE_cut = 0.002
    dCosTheta_cut = 0.3
    cosDTheta_cut = 0.2
    sinDTheta_cut = 1.0 #1.0 => no cut
    dTheta_cut = 0.2

    # df = df[df.dca < 100.0]
    df = df[numpy.isnan(df.delta_theta) == False]
    df["abs_dE"] = numpy.abs(df.dE)

    colnames = df.columns.values
    theta_M_events = df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events.to_csv('%s/theta_m.csv' % output_folder, index=False, header=False)

    theta_C_events = df[['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_events.to_csv('%s/theta_c.csv' % output_folder, index=False, header=False)

    df = df.sort_values(by=['sin_delta_theta'], ascending=[True])
    theta_M_events_cut_sinDtheta = df[df.sin_delta_theta < sinDTheta_cut][['E1','x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events_cut_sinDtheta.to_csv('%s/theta_m_cutSinDelta.csv' % output_folder, index=False, header=False)
    theta_M_events_cut_sinDtheta = df[df.sin_delta_theta < sinDTheta_cut][['E1','x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2', 'sin_delta_theta']]
    theta_M_events_cut_sinDtheta.to_csv('%s/theta_m_cutSinDelta_order.csv' % output_folder, index=False, header=False)

    theta_C_events_cut_sinDtheta = df[df.sin_delta_theta < sinDTheta_cut][['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_events_cut_sinDtheta.to_csv('%s/theta_c_cutSinDelta.csv' % output_folder, index=False, header=False)

    df = df.sort_values(by=['abs_dE'], ascending=[True])
    theta_M_events_cut = df[0:len(theta_M_events_cut_sinDtheta)]
    # theta_M_events_cut = df[df.abs_dE  < dE_cut][['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_M_events_cut[['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']].to_csv('%s/theta_c_cutDE.csv' % output_folder, index=False, header=False)

    # theta_M_events_cut = df[df.abs_dE < dE_cut][['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events_cut[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].to_csv('%s/theta_m_cutDE.csv' % output_folder, index=False, header=False)




#------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------


def usage():
    print "USAGE: %s [event file] [comma seperated energies in MeV] [output folder]" % (sys.argv[0])
    print "\nExample:\n %s CC_events.csv (1.17,1.33) ~/public_html/cc" % (sys.argv[0])
    sys.exit(-1)

def main():
    argv = sys.argv
    if not len(argv) == 4:
        usage()

    output_folder = sys.argv[3]
    energies = sys.argv[2].replace("(","").replace(")", "").split(",")
    energies = [float(e) for e in energies]

    print "Reading in csv file %s . . ." % sys.argv[1]
    df = read_csv_events(sys.argv[1])
    df = get_calculated_fields(df, energies)
    save_filtered_events(df, output_folder)
    df.to_csv("%s/dataframe.csv" % output_folder)

    print "Plotting histograms . . ."
    plot_data(df, output_folder)


if __name__ == "__main__":

    #cProfile.run("main()", sort=2)
    main()

