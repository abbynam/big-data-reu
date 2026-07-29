#!/usr/bin/python
"""
Description:
When reconstructing Co60 events from a Compton camera, there is missing information. First, we are
missing the labels for the gamma energies. Historically, we have attempted to fill in this information
based on the scatter energies. We are also missing information about the uncertainties of the detected
energies because we don't know the initial energy of the gamma in most cases. Further, we do not know
the proper order of the scatters.

In this script, we attempt to recover this missing information by framing the problem as an
Expectation Maximization problem. For the expectation step, we provide guesses of the energy labels.
For the maximization step, we maximize the parameters mu, and sigma that describe the distribution of
the uncertainties on E1. We use the known position of the gamma source to determine the order of the
Compton scatter events.
"""
__author__ = "Dennis Mackin <dsmackin@mdanderson.org>"
__date__ = "Aug. 16, 2016"
__version__ = "$Revision: 0.0.0$"

#-------- PYTHON IMPORT STATEMENTS -----------------------
import sys, os
import numpy as np
from math import sin, cos, pi
import pandas
# import matplotlib
# import matplotlib.pyplot as plt
import cProfile
# import re
# import csv
# import StringIO
import math
import itertools

from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

import pyximport; pyximport.install()
import utilities
import energy_matcher

def get_calculated_fields(df):

    df["E"] = df.E1 + df.E2
    df["cos_theta1"] = 1.0 - 0.511*(1.0/df.E2 - 1.0/df.E)
    df["theta1"] = np.arccos(df["cos_theta1"])

    return df


def plot_data(df, output_folder):

    # utilities.plot_1D(df.dca, 100, "DCA", "Counts", 'DCA(%.1f, %.1f, %.1f)' % (df.dca_x.iloc[1], df.dca_y.iloc[1], df.dca_z.iloc[1]), output_folder)
    # utilities.plot_2D(df.x1, "x1", 100, df.z1, "z1", 100, "x1_vs_z1", output_folder, True)

    return


def save_events(df, moniker, output_folder):
    df.to_csv('%s/%s_events.csv' % (output_folder, moniker), index=False, header=False)
    return df


def read_csv_events(file_path):
    return pandas.read_csv(file_path, names=('E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2'))


def filter_events(events):
    n0 = len(events)
    events = events[events[:, 4] > 0.04]
    events = events[events[:, 0] > 0.04]
    events = events[events[:, 4] < 1.5]
    events = events[events[:, 0] < 1.5]
    #
    # events = events[events[:, 0] + events[:, 4] < 1.5]
    # events = events[events[:, 0] + events[:, 4] > 0.1]

    n1 = len(events)
    print "%d of %d events passed energy cut" % (n1, n0)

    events = events[np.abs(events[:, 2] - events[:, 6]) > 0.5]
    print "%d of %d events passed dY cut" % (len(events), n1)

    return events


def make_plots(df, colnames, moniker, output_folder):

    em_events = df[df[:, colnames.DCA] >= 0.0]
    weights = em_events[:, colnames.WEIGHT]

    kleinnishina = np.array([energy_matcher.get_KN_xs(energy, theta) for energy, theta in
                    itertools.izip(em_events[:, colnames.EK], em_events[:, colnames.theta_k])])

    kleinnishina[np.isnan(kleinnishina)] = 0.0

    utilities.plot_1D_weighted(kleinnishina, em_events[:, colnames.WEIGHT], 100, "wkn_xs", "arbitrary",
                               'wknxs_%s' % moniker, output_folder)
    utilities.plot_1D(kleinnishina, 100, "kn_xs", "arbitrary", 'knxs_%s' % moniker, output_folder)
    utilities.plot_2D(em_events[:, colnames.E1], "E1", 100, kleinnishina, "knxs", 100, 'E1_v_knxs_%s' % moniker, output_folder, False)
    utilities.plot_2D_weighted(em_events[:, colnames.E1], "wE1", 100, kleinnishina, "knxs", 100, em_events[:, colnames.WEIGHT], 'wE1_v_knxs_%s' % moniker, output_folder, False)
    #

    utilities.plot_1D_weighted(em_events[:, colnames.DCA], em_events[:, colnames.WEIGHT], 100, "dca_weighted", "Counts",'dca_weighted_%s' % moniker, output_folder)
    utilities.plot_1D(em_events[:, colnames.DCA], 100, "DCA", "Counts", 'DCA_%s' % moniker, output_folder)

    utilities.plot_1D_weighted(em_events[:, colnames.E1], weights, 100, "wE1", "Counts", 'wE1_%s' % moniker, output_folder)
    utilities.plot_1D_weighted(em_events[:, colnames.E2], weights, 100, "wE2", "Counts", 'wE2_%s' % moniker, output_folder)
    utilities.plot_1D_weighted(em_events[:, colnames.Y2] - em_events[:, colnames.Y1], weights, 100, "wDy", "Counts", 'wDy_%s' % moniker, output_folder)

    dE1 = em_events[:, colnames.E1K] - em_events[:, colnames.E1]
    utilities.plot_1D_weighted(dE1[dE1*dE1 < 0.09], weights[dE1 * dE1 < 0.09], 100, "wdE1", "Counts",'wdE1_%s' % moniker, output_folder)

    dE2 = em_events[:, colnames.E2K] - em_events[:, colnames.E2]
    utilities.plot_1D_weighted(dE2, weights, 100, "wdE2", "Counts",'wdE2_%s' % moniker, output_folder)
    utilities.plot_2D(dE1, "dE1", 100, dE2, "dE2", 100, 'dE1_v_dE2_%s' % moniker, output_folder, True)
    utilities.plot_2D(em_events[:, colnames.E1], "E1", 100, dE2, "dE2", 100, 'E1_v_dE2_%s' % moniker, output_folder, False)
    utilities.plot_2D_weighted(em_events[:, colnames.E1], "E1", 50, dE1, "dE1", 50, weights, 'E1_v_dE1w_%s' % moniker, output_folder, False)
    utilities.plot_2D_weighted(em_events[:, colnames.E1], "E1", 50, dE2, "dE2", 50, weights, 'E1_v_dE2w_%s' % moniker, output_folder, False)

    utilities.plot_2D(em_events[:, colnames.E1], "E1", 100, em_events[:, colnames.E1K], "E1K", 100, 'E1_v_E1K_%s' % moniker, output_folder, False)
    utilities.plot_2D_weighted(em_events[:, colnames.E1], "E1", 100, em_events[:, colnames.E1K], "E1K", 100, weights, 'E1_v_E1Kw_%s' % moniker, output_folder, False)

    dE1_selected = dE1*dE1 < 0.2*0.2
    E1_selected = em_events[:, colnames.E1]
    E1_selected = E1_selected[dE1_selected]
    w_selected = weights[dE1_selected]
    utilities.plot_2D(E1_selected, "E1", 50, dE1[dE1_selected], "dE1", 50, 'E1_v_dE1_%s' % moniker, output_folder, False)
    utilities.plot_2D_weighted(E1_selected, "E1", 50, dE1[dE1_selected], "dE1", 50, w_selected, 'E1_v_dE1wtrimed_%s' % moniker, output_folder, False)
    utilities.plot_2D(em_events[:, colnames.E1], "E1", 50, dE2, "dE2", 50, 'E1_v_dE2_%s' % moniker, output_folder, False)

    try:
        utilities.plot_1D_weighted((180.0 / math.pi) * (em_events[:, colnames.theta_c] - em_events[:, colnames.theta_k]),
                               weights, 100, "wdThetack", "Counts", 'wdThetack_%s' % moniker, output_folder)
        utilities.plot_1D_weighted((180.0 / math.pi) * em_events[:, colnames.theta_k], weights, 100, "wtheta_k", "Counts",
                                   'wtheta_k_%s' % moniker, output_folder)
        utilities.plot_1D_weighted((180.0 / math.pi) * em_events[:, colnames.theta_c], weights, 100, "wtheta_c", "Counts",
                               'wtheta_c_%s' % moniker, output_folder)
        utilities.plot_1D((180.0 / math.pi) * (em_events[:, colnames.theta_c] - em_events[:, colnames.theta_k]), 100,
                          "dTheta", "Counts", 'dTheta_%s' % moniker, output_folder)
        utilities.plot_1D((180.0 / math.pi) * em_events[:, colnames.theta_k], 100, "theta_k", "Counts", 'theta_k_%s' % moniker,
                          output_folder)
        utilities.plot_1D((180.0 / math.pi) * em_events[:, colnames.theta_c], 100, "theta_c", "Counts", 'theta_c_%s' % moniker,
                          output_folder)
    except IndexError:
        print "nan in theta_c or theta_k"



    utilities.plot_1D(em_events[:, colnames.E1], 100, "E1", "Counts", 'E1_%s' % moniker, output_folder)
    utilities.plot_1D(em_events[:, colnames.E2], 100, "E2", "Counts", 'E2_%s' % moniker, output_folder)
    utilities.plot_1D(em_events[:, colnames.E1K] - em_events[:, colnames.E1], 100, "dE1", "Counts", 'dE1_%s' % moniker,
                      output_folder)
    utilities.plot_1D(em_events[:, colnames.E2K] - em_events[:, colnames.E2], 100, "dE2", "Counts", 'dE2_%s' % moniker,
                      output_folder)


    # # # Measured e1 and e2 often produce nan for the scattering angle -- Get rid of those nan events
    theta_m = em_events[:, colnames.theta_m]
    notnan_events = em_events[np.isnan(theta_m) == False]
    theta_k = notnan_events[:, colnames.theta_k]
    theta_m = notnan_events[:, colnames.theta_m]
    theta_c = notnan_events[:, colnames.theta_c]
    weights_notnan = notnan_events[:, colnames.WEIGHT]
    E1_notnan = notnan_events[:, colnames.E1]
    utilities.plot_1D((180.0 / math.pi) * theta_m, 100, "theta_m", "Counts", 'theta_m_%s' % moniker, output_folder)

    dthetakm = (180.0 / math.pi) * (theta_k - theta_m)
    events_for_dtheta = dthetakm * dthetakm < 50 * 50
    weight_dtheta = weights_notnan[events_for_dtheta]
    dtheta_selected = dthetakm[events_for_dtheta]
    utilities.plot_1D_weighted((180.0 / math.pi) *theta_m[events_for_dtheta], weight_dtheta, 100, "wtheta_m", "Counts", 'wtheta_m_%s' % moniker, output_folder)
    utilities.plot_1D(dthetakm, 100, "dthetakm", "Counts", 'dthetakm_%s' % moniker, output_folder)
    utilities.plot_1D_weighted(dtheta_selected, weight_dtheta, 100, "wdtheta", "Counts", 'wdthetakm_%s' % moniker, output_folder)

    utilities.plot_2D_weighted(E1_notnan, "E1", 25, dthetakm, "dthetakm", 25, weights_notnan, 'wE1_v_dthetakm_%s' % moniker, output_folder, False)

    dthetack = (180.0 / math.pi) * (theta_c - theta_k)
    events_for_dtheta = dthetack * dthetack < 50 * 50
    weight_dtheta = weights_notnan[events_for_dtheta]
    dtheta_selected = dthetack[events_for_dtheta]
    utilities.plot_1D(dthetack, 100, "dthetack", "arbitrary",'thetack_%s' % moniker, output_folder)
    utilities.plot_1D_weighted(dtheta_selected, weight_dtheta, 100, "wdtheta", "Counts", 'wdtheta_%s' % moniker, output_folder)
    utilities.plot_1D(dthetack, 100, "dthetack", "Counts", 'dthetack_%s' % moniker, output_folder)
    utilities.plot_2D_weighted(E1_notnan, "E1", 25, dthetack, "dthetack", 25, weights_notnan, 'wE1_v_dthetack_%s' % moniker, output_folder, False)


    dthetacm = (180.0 / math.pi) * (theta_c - theta_m)
    events_for_dtheta = dthetacm * dthetacm < 50 * 50
    weight_dtheta = weights_notnan[events_for_dtheta]
    dtheta_selected = dthetacm[events_for_dtheta]
    kn_selected = np.array([energy_matcher.get_KN_xs(energy, theta) for energy, theta in itertools.izip(notnan_events[:, colnames.EK], notnan_events[:, colnames.theta_k])])

    utilities.plot_2D_weighted(dthetacm, "dthetacm", 50, dthetack, "dthetack", 50, weights_notnan, 'dthetacm_v_dthetackw_%s' % moniker, output_folder, False)
    utilities.plot_2D(dthetacm, "dthetacm", 50, dthetack, "dthetack", 50, 'dthetacm_v_dthetack_%s' % moniker, output_folder, False)

    utilities.plot_1D(dthetacm, 100, "dthetacm", "Counts", 'dthetacm_%s' % moniker, output_folder)
    utilities.plot_1D_weighted(dtheta_selected, weight_dtheta, 100, "wdthetacm", "Counts", 'wdthetacm_%s' % moniker, output_folder)
    utilities.plot_1D(dthetacm, 100, "dthetacm", "Counts", 'dthetacm_%s' % moniker, output_folder)

    utilities.plot_2D_weighted(E1_notnan, "E1", 25, dthetacm, "dthetacm", 25, weights_notnan, 'wE1_v_dthetacm_%s' % moniker, output_folder, False)
    utilities.plot_2D_weighted((180.0 / math.pi) *theta_m, "wtheta_m", 50, kn_selected, "knxs", 50, weights_notnan, 'wthetam_v_knxs_%s' % moniker, output_folder, False)

    utilities.plot_2D(np.log(weights_notnan), "weight", 50, kn_selected, "knxs", 50, 'weight_v_knxs_%s' % moniker, output_folder, False)


def build_rf_dataframe(df):

    df = df[df['DCA'] > 0]
    X = df[['EK', 'E1', 'E2', 'X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2']]
    X['dY'] = df['Y2'] - df['Y1']
    X['dX'] = df['X2'] - df['X1']
    X['dZ'] = df['Z2'] - df['Z1']
    X['dS'] = np.sqrt(X['dX']*X['dX'] + X['dY']*X['dY'] + X['dZ']*X['dZ'])

    X['dE2'] = (df['EK'] - df['E1']) - df['E2']
    X['dTheta'] = df['theta_k'] - df['theta_m']

    Y = df['DCA']
    # Y = Y[Y > 0]
    # X = X[Y > 0, :]
    # Y.loc[Y < 0] = max(Y)

    nans = np.isnan(X['dTheta'])
    X.loc[nans, 'dTheta'] = math.pi

    return X, Y, df['WEIGHT']


def train_random_forest(df, output_folder):

    X_train, Y_train, weights = build_rf_dataframe(df)

    rfr = RandomForestRegressor(n_estimators=100, random_state=43, oob_score=True, n_jobs=-1, )
    print "Training random forest . . ."
    score = rfr.fit(X_train, Y_train)
    # score = rfr.fit(X_train[:1000], Y_train[:1000], weights.as_matrix()[:1000])
    print "\nTrained random forest . . ."
    print "OOB accuracy:", rfr.oob_score_
    print "OOB prediction:", rfr.oob_prediction_
    print "OOB score:", rfr.oob_score_
    print "Feature importances:\n\t%s\n\t%s\n" % (list(X_train.columns.values), rfr.feature_importances_)

    utilities.plot_2D(Y_train, "DCA", 100, rfr.oob_prediction_, "oobP(DCA)", 100, "oobDCAvsDCA", output_folder, False)
    utilities.plot_2D(Y_train, "DCA", 100, rfr.oob_prediction_, "oobP(DCA)", 100, "oobDCAvsDCA", output_folder, True)

    return rfr


def apply_random_forest(rf, df, output_folder, distance_param=10):

    df = df[df['DCA'] > 0]
    X, Y, weights = build_rf_dataframe(df)
    Y_predicted = rf.predict(X)

    utilities.plot_2D(Y, "DCA", 100, Y_predicted, "P(DCA)", 100, "DCAvsPDCA", output_folder, False)
    utilities.plot_2D(Y, "DCA", 100, Y_predicted, "P(DCA)", 100, "DCAvsPDCA", output_folder, True)

    df.loc[:, 'WEIGHT'] = (1.0/distance_param) * np.exp(-Y_predicted/distance_param)
    df.loc[:, 'PDCA'] = Y_predicted

    return df

#------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------


def usage():
    print "USAGE: %s [training event file] [training origin] [test event file] [test origin] [output folder]"  % (sys.argv[0])
    print "\nExample:\n %s CC_events1.csv 0.0,0.0,0.0  CC_events2.csv 0.0,0.0,6.0  ~/public_html/em" % (sys.argv[0])
    sys.exit(-1)


def main():
    if len(sys.argv) == 6:
        ftrain = sys.argv[1]
        otrain = [float(x) for x in sys.argv[2].split(",")]
        ftest = sys.argv[3]
        otest = [float(x) for x in sys.argv[4].split(",")]
    else:
        print "ERROR: Expected 6 parameters. Got %d . . ." % len(sys.argv)
        print "arguments: ", sys.argv
        usage()
    output_folder = sys.argv[5]

    print "Reading in csv file %s . . ." % sys.argv[1]

    em_events = read_csv_events(ftrain).as_matrix()
    em_events = filter_events(em_events)

    df_train = energy_matcher.run_known_origin_algorithm(em_events, np.array(otrain), np.array([1.17, 1.33]), np.array([0.5, 0.5]), 0.1)
    matrix_train = df_train.as_matrix()
    colnames = energy_matcher.column_names_enum()
    make_plots(matrix_train[matrix_train[:, colnames.EK] == 1.17], colnames, '117', output_folder)
    make_plots(matrix_train[matrix_train[:, colnames.EK] == 1.33], colnames, '133', output_folder)

    randomforest = train_random_forest(df_train, output_folder)

    df_test = read_csv_events(ftest)
    test_events = filter_events(df_test.as_matrix())
    df_test = energy_matcher.run_known_origin_algorithm(test_events, np.array(otest), np.array([1.17, 1.33]), np.array([0.5, 0.5]), 0.1)

    print "Applying random forest . . ."
    df_trained = apply_random_forest(randomforest, df_test, output_folder)
    df_trained.to_csv("%s/trained_all.csv" % output_folder)

    print "Sorting by PDCA . . ."
    df_trained.sort_values(by=['PDCA'], ascending=True, inplace=True)

    print "Saving em . . ."
    df_trained.to_csv("%s/trained_em.csv" % output_folder, columns=['E1', 'X1', 'Y1', 'Z1', 'E2', 'X2', 'Y2', 'Z2'], index=False)

    print "Saving EK . . ."
    df_trained['E2Corr'] = df_trained['EK'] - df_trained['E1']
    df_trained.to_csv("%s/trained_ek.csv" % output_folder, columns=['E1', 'X1', 'Y1', 'Z1', 'E2Corr', 'X2', 'Y2', 'Z2'], index=False)

    return

    em_events = df_train.as_matrix()
    colnames = energy_matcher.column_names_enum()
    index = colnames.EK
    tmp = em_events[1, colnames.EK]
    make_plots(em_events[em_events[:, colnames.EK] == 1.17], colnames, '117', output_folder)
    make_plots(em_events[em_events[:, colnames.EK] == 1.33], colnames, '133', output_folder)



    save_events(df_train, "orig", output_folder)
    save_events(pandas.DataFrame(em_events), "corr", output_folder)



if __name__ == "__main__":

    # cProfile.run("main()", sort=2)
    main()

