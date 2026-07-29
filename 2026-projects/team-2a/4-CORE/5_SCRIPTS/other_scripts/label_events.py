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

import pyximport; pyximport.install()
import utilities
import energy_matcher


def read_csv_events(file_path):
    df = pandas.read_csv(file_path, names = ('E1', 'x1', 'y1', 'z1', 'E2', 'x2', 'y2', 'z2'))

    return df


def get_compton_edge(E):
    return (2.0 * E * E) / (0.511 + 2.0 * E)


def get_calculated_fields(df, known_energies):

    if len(known_energies) > 2:
        raise Exception("ERROR: Filter events currently works with 2 energies (Co-60). You have %s: %s." % (len(known_energies), known_energies))

    df["E"] = df.E1 + df.E2
    df["cos_theta1"] = 1.0 - 0.511*(1.0/df.E2 - 1.0/df.E)
    df["theta1"] = numpy.arccos(df["cos_theta1"])
    df["E_known"] = known_energies[0]
    if len(known_energies) > 1:
        compton_edge = get_compton_edge(known_energies[0])
        # mask = df.loc[:, pandas.IndexSlice['E1', :]] >= compton_edge
        mask = df.ix[:, 'E1'] >= compton_edge
        df.ix[mask,'E_known'] = known_energies[1]
        mask = df.ix[:, 'E'] >= known_energies[0]*1.05
        df.ix[mask,'E_known'] = known_energies[1]

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

    df["DCA"] = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix(), point=numpy.array([0.0,0.0,32.7]))

    return df


def plot_data(df, output_folder):

    df = df[numpy.isnan(df.delta_theta) == False]

    utilities.plot_2D(df["sin_delta_theta"], "sin_dtheta", 100, df.dE, "dE", 100, "sin_dTheta_dE", output_folder, True)
    utilities.plot_2D(df["delta_theta"], "dtheta", 100, df.DCA, "DCA", 100, "dTheta_DCA", output_folder, True)
    utilities.plot_2D(df["dE"], "dE", 100, df.DCA, "DCA", 100, "dE_DCA", output_folder, True)
    utilities.plot_2D(df["delta_theta"], "dtheta", 100, df.E_known, "Ek", 100, "dtheta_Ek", output_folder, True)
    utilities.plot_2D(df["dE"], "dE", 100, df.E_known, "Ek", 100, "dE_Ek", output_folder, True)
    utilities.plot_2D(df.DCA, "dca", 100, df.E_known, "Ek", 100, "dca_Ek", output_folder, True)
    utilities.plot_2D(df[numpy.abs(df.dE) < 0.2 ].delta_theta, "delta_theta", 100, df[numpy.abs(df.dE) < 0.2 ].dE, "dE", 100, "dTheta_dE", output_folder, True)
    utilities.plot_2D(df.theta1, "theta", 100, df.theta_known, "theta known", 100, "theta_thetaTrue", output_folder, True)
    utilities.plot_2D(df.theta1, "theta", 100, df.theta_known, "theta known", 100, "theta_thetaTrue", output_folder, True)
    utilities.plot_2D(df.theta1, "theta", 100, df.E1, "E1", 100, "theta_E", output_folder, True)

    sin_dtheta = numpy.sin(df.delta_theta)
    utilities.plot_1D(sin_dtheta, 100, "sindtheta", "Counts", 'sindtheta', output_folder)
    utilities.plot_1D(df.DCA, 100, "DCA", "Counts", 'DCA', output_folder)
    dca = df.DCA[df.DCA < 30]
    dca = dca[dca >= 0]
    utilities.plot_1D(dca, 100, "DCA30", "Counts", 'DCA30', output_folder)
    utilities.plot_1D(df.E_known, 100, "E_known", "Counts", 'EKnown', output_folder)
    utilities.plot_1D(df.E, 100, "E_m", "Counts", 'Emeasured', output_folder)
    utilities.plot_1D(df.E2_known, 100, "E2_known", "Counts", 'E2_Known', output_folder)
    utilities.plot_1D(df.dE, 100, "dE", "Counts", 'dE', output_folder)
    utilities.plot_1D(df[numpy.abs(df.cos_theta_known) <= 1.0].cos_theta_known, 100, "cos_theta_known", "Counts", 'cos_theta_known', output_folder)
    utilities.plot_1D(df.theta_known, 100, "theta_known", "Counts", 'theta_known', output_folder)
    utilities.plot_1D(df.theta1, 100, "theta_measured", "Counts", 'theta_measured', output_folder)
    utilities.plot_1D(df.delta_theta, 100, "delta_theta", "Counts", 'delta_theta', output_folder)
    utilities.plot_1D(df.E1, 100, "E1", "Counts", 'E1', output_folder)
    utilities.plot_1D(df.E2, 100, "E2", "Counts", 'E2', output_folder)


def save_filtered_events(df, output_folder, point=numpy.array([0.0,0.0,0.0])):

    dE_cut = 2.0
    dCosTheta_cut = math.pi
    cosDTheta_cut = math.pi
    dTheta_cut = math.pi
    dTheta_cut = math.pi

    df = df[numpy.isnan(df.delta_theta) == False]

    colnames = df.columns.values
    theta_M_events = df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events.to_csv('%s/theta_m.csv' % output_folder, index=False, header=False)

    theta_C_events = df[['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_events.to_csv('%s/theta_c.csv' % output_folder, index=False, header=False)

    df.ix[:,'abs_dE'] = numpy.abs(df.dE)
    df = df.sort_values(by=['abs_dE'], ascending=[True])
    theta_M_dE = df[df.abs_dE < dE_cut][['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_dE.to_csv('%s/theta_m_dE.csv' % output_folder, index=False, header=False)
    theta_C_dE = df[df.abs_dE < dE_cut][['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_M_dE.to_csv('%s/theta_c_dE.csv' % output_folder, index=False, header=False)

    # df["abs_sin_delta_theta"] = numpy.abs(df.sin_delta_theta)
    # df = df.sort_values(by=['abs_sin_delta_theta'], ascending=[True])

    #DSM after deltaTheta gets bigger than pi/2, sin(delta theta) starts getting smaller.
    #  Therefore, we sort on delta_theta instead.
    df["abs_delta_theta"] = numpy.abs(df.delta_theta)
    df = df.sort_values(by=['abs_delta_theta'], ascending=[True])
    theta_M_dTheta = df[df.abs_delta_theta < dTheta_cut][['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_dTheta.to_csv('%s/theta_m_dTheta.csv' % output_folder, index=False, header=False)

    theta_C_dTheta = df[df.abs_delta_theta < dTheta_cut][['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_dTheta.to_csv('%s/theta_c_dTheta.csv' % output_folder, index=False, header=False)

    df["DCA"] = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix(), point)
    df = df.sort_values(by=['DCA'], ascending=[True])
    theta_m_dca = df[['E1','x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_m_dca.to_csv('%s/theta_m_dca.csv' % output_folder, index=False, header=False)

    df["DCA"] = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']].as_matrix(), point=numpy.array([0.0,0.0,0.0]))
    df = df.sort_values(by=['DCA'], ascending=[True])
    theta_c_dca = df[['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_c_dca.to_csv('%s/theta_c_dca.csv' % output_folder, index=False, header=False)


def save_labeled_events(df, output_folder):

    events = df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    events.to_csv('%s/theta_idealunsorted_m.csv' % output_folder, index=False, header=False)

    df["E2k"] = df.Ek - df.E1
    events = df[['E1', 'x1', 'y1', 'z1', "E2k", 'x2', 'y2', 'z2']]
    events.to_csv('%s/theta_idealunsorted_k.csv' % output_folder, index=False, header=False)

    df = df.sort_values(by=['DCA'], ascending=[True])
    events = df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    events.to_csv('%s/theta_ideal_m.csv' % output_folder, index=False, header=False)

    events = df[['E1', 'x1', 'y1', 'z1', "E2k", 'x2', 'y2', 'z2']]
    events.to_csv('%s/theta_ideal_k.csv' % output_folder, index=False, header=False)

    df["dE"] = df.Ek - df.E1 -df.E2
    mask = df.dE*df.dE > 0.09
    events = df.loc[mask, :]
    events_recov = events[['E1', 'x1', 'y1', 'z1', "E2k", 'x2', 'y2', 'z2']]
    events_recov.to_csv('%s/theta_recovered_k.csv' % output_folder, index=False, header=False)

    events = events[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    events.to_csv('%s/theta_recovered_m.csv' % output_folder, index=False, header=False)


def put_scatters_in_order(df, energies, point=numpy.array([0.0,0.0,0.0])):
    mask = df.E1 < df.E2
    df.ix[mask, ['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']] = df.ix[mask, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]]

    dca = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix(), point)
    mask = dca < 0.0
    df.ix[mask, ['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']] = df.ix[mask, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]]

    #DSM  this correction seems to add bias for PJ2 point source data (2016-12-05)
    # dca = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix(), point)
    # mask = dca > 200.0
    # df.ix[mask, ['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']] = df.ix[mask, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]]

    return df


def label_events(df, energies, dca_point=numpy.array([0.0,0.0,0.0])):

    df = df.copy()
    positions = [0,1]
    # dca_list = []
    dca_array = numpy.zeros([len(df), len(energies)*2])
    energy_array = numpy.repeat(energies, len(df), axis=0)
    energy_array = energy_array.reshape((len(df), len(energies)), order='F')

    for i, energy in enumerate(energies):
        for j, position in enumerate(positions):
            if position % 2 == 0:
                tmp_array = df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix()
            else:
                tmp_array = df[["E2", 'x2', 'y2', 'z2', 'E1', 'x1', 'y1', 'z1']].as_matrix()

            #Set the energy to the known value
            tmp_array[:,4] = energy - tmp_array[:,0]

            index = len(energies) * i + j
            dca_array[:, index] = energy_matcher.get_dca_for_events(tmp_array, point=dca_point)

    dca_array[dca_array < 0.0] = 1.0E6

    min_vals = numpy.amin(dca_array, axis=1, keepdims=True)
    min_index = numpy.argmin(dca_array, axis=1)
    df["DCA"] = min_vals
    df["order_flipped"] = min_index % 2
    df.ix[df["order_flipped"], ['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']] = df[["E2", 'x2', 'y2', 'z2', 'E1', 'x1', 'y1', 'z1']]

    df["Ek"] = energies[0]
    for i,e in enumerate(energies):
        df.ix[min_index/2 == i, "Ek"] = e

    return df



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

    dca_point = numpy.array([0.0,0.0,0.0])
    output_folder = sys.argv[3]
    energies = sys.argv[2].replace("(","").replace(")", "").split(",")
    energies = numpy.array([float(e) for e in energies])

    print "Reading in csv file %s . . ." % sys.argv[1]
    df = read_csv_events(sys.argv[1])
    df_labeled = label_events(df, energies)
    save_labeled_events(df_labeled, output_folder)

    df = put_scatters_in_order(df, energies, dca_point)
    df = get_calculated_fields(df, energies)
    save_filtered_events(df, output_folder)

    print "Plotting histograms . . ."
    plot_data(df, output_folder)


if __name__ == "__main__":
    # cProfile.run("main()", sort=2)
    main()

