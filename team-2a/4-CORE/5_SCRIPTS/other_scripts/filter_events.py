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


def get_E_known(E, known_energies):
    #DSM -- 2017-03-20 This method does not consider the compton edge.

    #if there's only one energy, select it
    if len(known_energies) == 1:
        return E*0.0 + known_energies[0]

    known_energies = numpy.unique(known_energies)

    ke_array = numpy.repeat(known_energies, len(E))
    ke_array = ke_array.reshape((len(E), len(known_energies)), order='F')

    # E_array = numpy.tile(E, (1, len(known_energies)))[0]
    E_array = numpy.tile(E, len(known_energies))
    E_array = E_array.reshape((len(E), len(known_energies)), order='F')

    diff_array = numpy.abs((ke_array - E_array))

    min_vals = numpy.amin(diff_array, axis=1, keepdims=True)
    min_truths = (diff_array == min_vals)
    print "truths = ", sum(sum(min_truths))

    return E_known


def get_compton_edge(E):
    return (2.0 * E * E) / (0.511 + 2.0 * E)


def get_calculated_fields(df, known_energies):

    if len(known_energies) > 2:
        raise Exception("ERROR: Filter events currently works with 2 energies (Co-60). You have %s: %s." % (len(known_energies), known_energies))

    df["E"] = df.E1 + df.E2
    df["cos_theta1"] = 1.0 - 0.511*(1.0/df.E2 - 1.0/df.E)
    df["theta1"] = numpy.arccos(df["cos_theta1"])
    df["E_known"] = known_energies[0]
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
    # utilities.plot_1D(df.cos_delta_theta, 100, "cos_delta_theta", "Counts", 'cos_delta_theta', output_folder)
    # utilities.plot_1D(df.delta_theta, 1000, "delta_theta", "Counts", 'delta_theta', output_folder)
    #utilities.plot_1D(df[numpy.abs(df.cos_theta_known) <= 1.0].delta_cos_theta, 100, "delta_cos_theta", "Counts", 'delta_cos_theta', output_folder)
    # utilities.plot_1D(df.delta_cos_theta, 100, "delta_cos_theta", "Counts", 'delta_cos_theta', output_folder)


    # utilities.plot_1D(df.cos_theta_known, 100, "cos_theta_known", "Counts", 'CosThetaKnown', output_folder)
    # utilities.plot_1D(df.theta_117, 100, "Theta117", "Counts", 'Theta117', output_folder)
    # utilities.plot_1D(df.theta_133, 100, "Theta133", "Counts", 'Theta133', output_folder)
    # utilities.plot_1D(df.delta_theta, 100, "DeltaTheta", "Counts", 'DeltaTheta', output_folder)
    # utilities.plot_1D(df.delta_cos_theta, 100, "dCosTheta", "Counts", 'DeltaCosTheta', output_folder)
    #
    # utilities.plot_1D(df.theta_known, 100, "theta_known", "Counts", 'ThetaKnown', output_folder)
    # utilities.plot_1D(df.dca, 100, "DCA", "Counts", 'DCA(%.1f, %.1f, %.1f)' % (df.dca_x[1], df.dca_y[1], df.dca_z[1]), output_folder)
    # utilities.plot_1D(df.theta1, 100, "theta1", "Counts", '1st Scattering Angle', output_folder)
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

    print "%d events to save . . ." % (len(df))
    dE_cut = 2.0
    dCosTheta_cut = math.pi
    cosDTheta_cut = math.pi
    dTheta_cut = math.pi
    dTheta_cut = math.pi

    df = df[numpy.isnan(df.delta_theta) == False]

    print "%d events are not nan . . ." % (len(df))

    colnames = df.columns.values
    theta_M_events = df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_events.to_csv('%s/theta_m.csv' % output_folder, index=False, header=False)

    theta_C_events = df[['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_events.to_csv('%s/theta_c.csv' % output_folder, index=False, header=False)

    df['abs_dE'] = numpy.abs(df.dE)
    df = df.sort_values(by=['abs_dE'], ascending=[True])
    theta_M_dE = df[df.abs_dE < dE_cut][['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_dE.to_csv('%s/theta_m_dE.csv' % output_folder, index=False, header=False)
    theta_C_dE = df[df.abs_dE < dE_cut][['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_dE.to_csv('%s/theta_c_dE.csv' % output_folder, index=False, header=False)

    # df["abs_sin_delta_theta"] = numpy.abs(df.sin_delta_theta)
    # df = df.sort_values(by=['abs_sin_delta_theta'], ascending=[True])

    #DSM after deltaTheta gets bigger than pi/2, sin(delta theta) starts getting smaller.
    #  Therefore, we sort on delta_theta instead.
    df["abs_delta_theta"] = numpy.abs(df.delta_theta)
    df = df.sort_values(by=['abs_delta_theta'], ascending=[True])
    theta_M_dTheta = df[df.abs_delta_theta < dTheta_cut][['E1','x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_M_dTheta.to_csv('%s/theta_m_dTheta.csv' % output_folder, index=False, header=False)

    theta_C_dTheta = df[df.abs_delta_theta < dTheta_cut][['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_C_dTheta.to_csv('%s/theta_c_dTheta.csv' % output_folder, index=False, header=False)

    df["DCA"] = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix(), point=numpy.array([0.0,0.0,0.0]))
    df = df.sort_values(by=['DCA'], ascending=[True])
    theta_m_dca = df[['E1','x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']]
    theta_m_dca.to_csv('%s/theta_m_dca.csv' % output_folder, index=False, header=False)

    df["DCA"] = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']].as_matrix(), point=numpy.array([0.0,0.0,0.0]))
    df = df.sort_values(by=['DCA'], ascending=[True])
    theta_c_dca = df[['E1','x1', 'y1', 'z1', "E2_known", 'x2', 'y2', 'z2']]
    theta_c_dca.to_csv('%s/theta_c_dca.csv' % output_folder, index=False, header=False)


def put_scatters_in_order(df):

    # make first scatter the higher energy scatter
    mask = df.E1 < df.E2
    print "E1 is smaller than E2 %.3f of events . . ." % (float(sum(mask))/len(mask))
    df.ix[mask] = df.ix[mask, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]]
    mask = df.E1 < df.E2
    print "E1 is smaller than E2 %.3f of events . . ." % (float(sum(mask))/len(mask))

    #if the DCA < 0 then swap the events to try the other order.
    dca = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix(), point=numpy.array([0.0,0.0,0.0]))
    mask = numpy.asarray(dca) < 0
    print "DCA < 0 for %.3f of events . . ." % (float(sum(mask))/len(mask))
    df.ix[mask] = df.ix[mask, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]]
    dca = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']].as_matrix(), point=numpy.array([0.0,0.0,0.0]))
    mask = numpy.asarray(dca) < 0
    print "DCA < 0 for %.3f of events . . ." % (float(sum(mask))/len(mask))

    return df


def filter_events(df):
    events = df.as_matrix()
    n0 = len(events)
    events = events[events[:, 4] > 0.01]
    events = events[events[:, 0] > 0.01]
    events = events[events[:, 4] < 1.5]
    events = events[events[:, 0] < 1.5]
    #
    # events = events[events[:, 0] + events[:, 4] < 1.5]
    # events = events[events[:, 0] + events[:, 4] > 0.1]

    n1 = len(events)
    print "%d of %d events passed energy cut" % (n1, n0)

    events = events[numpy.abs(events[:, 2] - events[:, 6]) > 0.5]
    print "%d of %d events passed dY cut" % (len(events), n1)

    df_new = pandas.DataFrame(events, columns=df.columns.values)

    return df_new


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
    df = filter_events(df)

    print "%d events after rough filtering  . . ." % (len(df))
    df = put_scatters_in_order(df)
    df = get_calculated_fields(df, energies)
    save_filtered_events(df, output_folder)

    print "Plotting histograms . . ."
    plot_data(df, output_folder)



if __name__ == "__main__":

    # cProfile.run("main()", sort=2)
    main()

