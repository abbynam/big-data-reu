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

#-------- GLOBAL VARIABLES -----------------------
BIG_DCA = 700


#-------- PYTHON IMPORT STATEMENTS -----------------------
import sys, os
import numpy as np
from math import sin, cos, pi
import matplotlib
import matplotlib.pyplot as plt
import pandas
import math
import itertools
import collections
import numpy
import cPickle
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

import pyximport; pyximport.install(setup_args={'include_dirs': numpy.get_include()})
import utilities
import energy_matcher

def save_events(df, moniker, output_folder):
    df.to_csv('%s/%s_events.csv' % (output_folder, moniker), index=False, header=False)
    return df


def read_csv_events(file_path):
    return pandas.read_csv(file_path, names=('E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2'))


def build_training_rf_dataframe(df):

    predictors = build_rf_dataframe(df)
    response = df.DCAk

    return predictors, response


def build_rf_dataframe(df):

    # df.theta_M[pandas.isnull(df.theta_M)] = -100
    # df.theta_K[pandas.isnull(df.theta_K)] = -100
    df.loc[pandas.isnull(df.theta_M), ['theta_M']] = -100
    df.loc[pandas.isnull(df.theta_K), ['theta_K']] = -100
    # predictors = df[['Ek', 'E1', 'E2', 'x1', 'y1', 'z1', 'x2', 'y2', 'z2', 'theta_M', 'theta_K']]
    predictors = df[['Ek', 'E1', 'E2', 'x1', 'x2','y1', 'y2', 'z1', 'z2', 'theta_M', 'theta_K']]
    predictors.loc[:, 'dY'] = df['y2'] - df['y1']
    # predictors.loc[:, 'dX'] = df['x2'] - df['x1']
    # predictors.loc[:, 'dZ'] = df['z2'] - df['z1']
    predictors.loc[:, 'dS'] = np.sqrt((df.x1 - df.x2)*(df.x1 - df.x2) + (df.y1 - df.y2)*(df.y1 - df.y2) + (df.z1 - df.z2)*(df.z1 - df.z2))

    predictors.loc[:, 'dTheta'] = df['theta_K'] - df['theta_M']
    predictors.loc[:, 'dE'] = df.Ek - df.E1 - df.E2
    predictors.loc[:, 'dE1E2'] = df.E1 - df.E2

    # predictors.loc[:, 'Ewindow'] = (df.Ek - df.E1 - df.E2 < predictors.Ek * 0.10)


    #DSM variables are not useful so I removed them 2017-04-05
    # predictors.loc[:, 'x1_abs'] = numpy.abs(df.x1)
    # predictors.loc[:, 'x2_abs'] = numpy.abs(df.x2)
    # predictors.loc[:, 'z1_abs'] = numpy.abs(df.z1)
    # predictors.loc[:, 'z2_abs'] = numpy.abs(df.z2)

    for p in list(predictors):
        # print "Checking for nan: ", p, " . . ."
        assert(pandas.notnull(predictors[p]).all())

    return predictors


def get_compton_edge(E):
    return (2.0 * E * E) / (0.511 + 2.0 * E)


def label_events_ideal(event_permutations, output_folder=".",
                       distance_param=10, verbose=False):
    num_events = event_permutations[0].shape[0]
    num_permutations = len(event_permutations)

    dca_arr = numpy.zeros([num_events, num_permutations])

    for i, df in enumerate(event_permutations):
        # print "Setting column %d to %.3f, Ek = %.3f" % (i, numpy.mean(df.DCAk), numpy.mean(df.Ek))
        dca_arr[:, i] = df.DCAk

    nans_test = np.isnan(dca_arr)
    best_dcas = numpy.argmin(dca_arr, axis=1)
    #
    # print "avg_best_dcas = %.3f" % numpy.mean(best_dcas)
    # for i in range(30):
    #     print "BESTS", best_dcas[i], dca_arr[i,:]

    #Make initial data permutation 0. Copy in the events from other permutations if they have the
    # minimum predicted DCA values
    df = event_permutations[0]
    for i in range(1, num_permutations):
        df.ix[best_dcas == i, :] = event_permutations[i].ix[best_dcas == i, :].values

    # print df.head(30).DCAk

    n0 = df.shape[0]
    any_na1 = df.isna().any()
    df.theta_M[df.theta_M.isna()] = 0.0
    df.dTheta[df.dTheta.isna()] = math.pi
    any_na2 = df.isna().any()
    df = df.dropna()
    if verbose: print "Removed non-physical double scatters (%d/%d = %.2f)  . . ." \
                      % (n0 - df.shape[0], n0, float(n0 - df.shape[0])/float(n0))
    df_small = df[df.DCAk < 100]
    # utilities.plot_1D(df.DCAk, 100, "DCA_test", "Counts", "DCA_test_ideal", output_folder)
    # utilities.plot_1D(df_small.DCAk, 100, "DCA_test_small", "Counts", "DCA_test_small_ideal", output_folder)
    # utilities.plot_1D(df_small.theta_M - df_small.theta_K, 100, "dtheta_test", "Counts", "dtheta_test_ideal", output_folder)
    # utilities.plot_1D(df_small.theta_M, 100, "thetaM_test", "Counts", "thetaM_test_ideal", output_folder)
    # utilities.plot_1D(df_small.theta_K, 100, "thetaK_test", "Counts", "thetaK_test_ideal", output_folder)
    # utilities.plot_1D(df_small.Ek, 100, "Ek_test", "Counts", "Ek_test_ideal", output_folder)
    # utilities.plot_1D(df_small.y2 - df_small.y1, 100, "dY_test", "Counts", "dY_test_ideal", output_folder)
    # utilities.plot_1D(df_small.E2 - df_small.E1, 100, "dE_test", "Counts", "dE_test_ideal", output_folder)
    # utilities.plot_1D(df_small.E2 + df_small.E1, 100, "EM_test", "Counts", "EM_test_ideal", output_folder)

    return df


def get_known_energies(E1, E2, energies, cross_sections = [0.5, 0.5]):
    '''
    Note: This version of the algorithm only works for 2 energies of the same
    cross section. For more energies implement an algorithm in cython where
        1) Detemine the physically allowed energies.
        2) Re-calculated the probabilities based on the allowed energies and cross sections
        3) Randomly choose one of the energies.
    Since we may never use this method, I implemented a simpler solution for Co60 which has only two energies with
    the same cross section.
    :param E1: Scatter 1 measured energy
    :param E2: Scatter 2 measured energy
    :param energies: Scatter 3
    :param cross_sections:
    :return:
    '''

    energies.sort()

    Ek = E1*0.0 - 100.0 #if E1 is vector Ek is vector
    Em = E1 + E2
    Ek_max = numpy.max(energies)

    #Set energy equal to lowest energy that is physically allowed
    for e in reversed(energies):
        ce = get_compton_edge(e)
        Ek[(E1 < ce) & (Em <= e * 1.05)] = e

    random_mask = (numpy.random.random_sample(len(E1)) < 0.5)

    #If total energy is close to a known energy, then choose that known energy
    for e in reversed(energies):
        ce = get_compton_edge(e)
        Ek[(E1 < ce) & ((Em - e)*(Em - e) <= (0.05 * e)*(0.05 * e))] = e

    return Ek


def label_events_random_forest(rf, event_permutations, output_folder, distance_param=100):
    '''
    Uses the random forest
    :param rf: random forest that predicts DCA values
    :param event_permutations: A list of data frames, each with one of the possible combinations of event energies and scatter sequences
    :param output_folder:
    :param distance_param: Arbitrary value that is used for event weighting. 10 (1 cm) is probably OK.
    :return: dataframe with the best predicted event energies and scatter sequences.
    '''
    num_events = event_permutations[0].shape[0]
    num_permutations = len(event_permutations)
    predictions = numpy.zeros([num_events, num_permutations])
    dca_arr = numpy.zeros([num_events, num_permutations])
    dca_arr = numpy.zeros([num_events, num_permutations])

    for i, df in enumerate(event_permutations):
        X = build_rf_dataframe(df)
        predictions[:,i] = rf.predict(X)
        dca_arr[:, i] = df.DCAk

    min_dcas = numpy.argmin(dca_arr, axis=1)
    min_dcas_pred = numpy.argmin(predictions, axis=1)
    best_dcas = numpy.amin(dca_arr, axis=1)
    best_dcas_pred = numpy.amin(predictions, axis=1)

    #Make initial data permutation 0. Copy in the events from other permutations if they have the
    # minimum predicted DCA values
    df = event_permutations[0]
    for i in range(1, num_permutations):
        df.ix[min_dcas_pred == i, :] = event_permutations[i].ix[min_dcas_pred == i, :].values

    #------ HOW GOOD DID WE DO AT PICKING BEST PERMUTATION? ------------
    correct_selection = (min_dcas == min_dcas_pred)
    dca_diff = df.DCAk.as_matrix() - best_dcas

    bad_miss = 0
    good_get = 1
    dca_cutoff = 50.0
    for i in range(len(correct_selection)):
        if (not correct_selection[i]) and (best_dcas[i] < dca_cutoff):
            # print "DCA comp,", i, df.E.as_matrix()[i], df.Ek.as_matrix()[i], df.DCAk.as_matrix()[i], best_dcas[i], dca_diff[i], predictions[i], dca_arr[i]
            bad_miss += 1
        elif correct_selection[i] and (best_dcas[i] < dca_cutoff):
            good_get += 1
    print ("Mislabelled %.2f of (DCA < %.0f) events (%d/%d)." % (float(bad_miss)/float(bad_miss + good_get), dca_cutoff, bad_miss, bad_miss + good_get))


    df.loc[:, 'WEIGHT'] = (1.0/distance_param) * np.exp(-best_dcas_pred/distance_param)
    df.loc[:, 'PDCA'] = best_dcas_pred

    # df_small = df[df.DCAk < 100]
    # utilities.plot_2D(df_small.DCAk, "DCA", 100, df_small.PDCA, "P_DCA", 100, "point_DCAvsPDCA", output_folder, False)
    # utilities.plot_2D(df_small.DCAk, "DCA", 100, df_small.PDCA, "P_DCA", 100, "point_DCAvsPDCA", output_folder, True)
    # utilities.plot_1D(df.DCAk, 100, "DCA_test", "Counts", "DCA_test", output_folder)
    # utilities.plot_1D(df_small.DCAk, 100, "DCA_test_small", "Counts", "DCA_test_small", output_folder)
    # utilities.plot_1D(df_small.PDCA, 100, "DCA_pred", "Counts", "DCA_test_pred", output_folder)
    # utilities.plot_1D(df_small.theta_M - df_small.theta_K, 100, "dtheta_test", "Counts", "dtheta_test", output_folder)
    # utilities.plot_1D(df_small.theta_M, 100, "thetaM_test", "Counts", "thetaM_test", output_folder)
    # utilities.plot_1D(df_small.theta_K, 100, "thetaK_test", "Counts", "thetaK_test", output_folder)
    # utilities.plot_1D(df_small.Ek, 100, "Ek_test", "Counts", "Ek_test", output_folder)
    # utilities.plot_1D(df_small.y2 - df_small.y1, 100, "dY_test", "Counts", "dY_test", output_folder)
    # utilities.plot_1D(df_small.E2 - df_small.E1, 100, "dE_test", "Counts", "dE_test", output_folder)
    # utilities.plot_1D(df_small.E2 + df_small.E1, 100, "EM_test", "Counts", "EM_test", output_folder)

    return df


def train_random_forest(df, output_folder, rf_file_name, rebuild_forest=False, distance_param=20):


    if not rebuild_forest:
        try:
            print "Unpickling (%s) . . ." % rf_file_name
            rfr = cPickle.load(open(rf_file_name, 'rb'))
            print "Using previously trained random forest (%s) . . ." % rf_file_name
            if rfr != None: return rfr
        except:
            print "File ", rf_file_name, " not found. Rebuilding forest . . ."

    X_train, Y_train = build_training_rf_dataframe(df)
    #print X_train.corr(method='spearman')
    weights = 1.0 * numpy.exp(-0.5 * numpy.array(Y_train)/distance_param)
    # print "range X0:", min(df.X0), max(df.X0)
    # print "range Y0:", min(df.Y0), max(df.Y0)
    # print "range Z0:", min(df.Z0), max(df.Z0)
    # utilities.plot_1D(df.X0, 100, "X0", "Counts", "X0", output_folder)
    # utilities.plot_1D(df.Y0, 100, "Y0", "Counts", "Y0", output_folder)
    # utilities.plot_1D(df.Z0, 100, "Z0", "Counts", "Z0", output_folder)
    # utilities.plot_1D(weights, 100, "trainingweights", "Counts", "weights", output_folder)
    # utilities.plot_1D(numpy.log2(1 + weights), 100, "logweights", "Counts", "weightslog", output_folder)
    # utilities.plot_1D(Y_train, 100, "Y_train", "Counts", "Y_train", output_folder)
    # utilities.plot_1D_weighted(Y_train, weights, 100, "wY_train", "Counts",'ytrain', output_folder)

    # m = X_train.shape[1]/3
    m = X_train.shape[1]/3
    print "Selecting ", m, " features at a time . . ."
    rfr = RandomForestRegressor(n_estimators=100, max_features=m, random_state=43, oob_score=True, n_jobs=19, min_samples_leaf=5)
    # rfr = RandomForestRegressor(n_estimators=100, max_features=m, random_state=43, oob_score=False,
    #                             n_jobs=15, min_samples_leaf=5)
    print "Training random forest . . ."
    score = rfr.fit(X_train, Y_train, sample_weight=weights)
    # score = rfr.fit(X_train, Y_train)

    # print "\nTrained random forest . . ."
    # print "OOB accuracy:", rfr.oob_score_
    # print "OOB prediction:", rfr.oob_prediction_
    # print "OOB score:", rfr.oob_score_
    # print "Feature importances:\n\t%s\n\t%s\n" % (list(X_train.columns.values), rfr.feature_importances_)

    # utilities.plot_2D(Y_train, "DCA", 100, rfr.oob_prediction_, "oobP_DCA", 100, "oobDCAvsDCA", output_folder, False)
    # utilities.plot_2D(Y_train, "DCA", 100, rfr.oob_prediction_, "oobP_DCA", 100, "oobDCAvsDCA", output_folder, True)
    # utilities.plot_1D(df.DCA[df.DCA > 0], 100, "DCA_train", "Counts", "DCA_train", output_folder)
    # utilities.plot_1D(Y_train, 100, "DCAk_train", "Counts", "DCAk_train", output_folder)
    # utilities.plot_1D(rfr.oob_prediction_, 100, "DCA_train", "Counts", "DCA_train_oob", output_folder)
    # utilities.plot_1D(X_train.dTheta, 100, "dTheta_train", "Counts", "dtheta_train", output_folder)
    # utilities.plot_1D(X_train.theta_M, 100, "thetaM_train", "Counts", "thetaM_train", output_folder)
    # utilities.plot_1D(X_train.theta_K, 100, "thetaK_train", "Counts", "thetaK_train", output_folder)
    # utilities.plot_1D(X_train.Ek, 100, "Ek_train", "Counts", "Ek_train", output_folder)
    # utilities.plot_1D(X_train.dY, 100, "dY_train", "Counts", "dY_train", output_folder)

    # df_2d = pandas.DataFrame({'DCA':Y_train.as_matrix(), 'prediction': rfr.oob_prediction_}, columns = ['DCA', 'prediction'])

    # df_2d = df_2d[df_2d.DCA < 100]
    # utilities.plot_2D(df_2d.DCA, "DCA", 100, df_2d.prediction, "oobP_DCA", 100, "predictionvsDCA_small", output_folder, True)
    # utilities.plot_2D(df_2d.DCA, "DCA", 100, df_2d.prediction, "oobP_DCA", 100, "predictionvsDCA_small", output_folder, False)

    cPickle.dump(rfr, open(rf_file_name, 'wb'), protocol=-1)

    return rfr


def read_training_events(file_path):
    return energy_matcher.read_training_events(file_path)


def make_efficiency_plots(df, output_folder):

    df = df[numpy.abs(df.X0) < 40]
    df = df[numpy.abs(df.Y0) < 0.1]
    df = df[numpy.abs(df.Z0) < 130]

    utilities.plot_1D(df.DCA,  100, "DCA", "arbitrary", 'DCA', output_folder)
    utilities.plot_1D(df.DCAk, 100, "DCA Known", "arbitrary", 'DCAk', output_folder)
    utilities.plot_1D(df.X0, 100, "x0 emission", "arbitrary", 'x0 emission', output_folder)
    utilities.plot_1D(df.Y0, 100, "y0 emission", "arbitrary", 'y0 emission', output_folder)
    utilities.plot_1D(df.Z0, 100, "z0 emission", "arbitrary", 'z0 emission', output_folder)

    utilities.plot_2D(df.Z0, "Z0", 50, df.X0, "X0", 50, 'Emission Plane', output_folder, False)

    df_good = df[df.DCAk < 10]
    utilities.plot_2D(df_good.Z0, "Z0", 30, df_good.X0, "X0", 30, 'DCA < 10 Emission Plane', output_folder, False)
    utilities.plot_1D(df_good.X0, 100, "x0", "arbitrary", 'x0_DCA10', output_folder)

    utilities.plot_1D(df_good.DCA, 100, "DCA", "arbitrary", 'DCA10', output_folder)
    utilities.plot_1D(df_good.DCAk, 100, "DCA Known", "arbitrary", 'DCAk10', output_folder)


def make_training_sample(df, energies = [1.17, 1.33]):
    '''
    Need to create a sample that is not too heavily biased toward bad events
    '''
    permutations = get_event_permutations(df, energies)
    permutations = map(get_calculated_fields_train, permutations)
    permutations = [p.dropna() for p in permutations]

    df_train = pandas.concat([permutations[0],permutations[2],permutations[1],permutations[3]])

    return df_train


def get_event_permutations(df, energies = [1.17, 1.33]):
    '''

    :param df:
    :return: df with events permuted for each energy and scatter sequence
    '''
    p1 = df
    p2 = df.copy()
    p2.ix[:, ['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']] = p1.ix[:, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]].values
    p3 = p1.copy()
    p4 = p2.copy()
    p1['Ek'] = energies[0]
    p2['Ek'] = energies[0]
    p3['Ek'] = energies[1]
    p4['Ek'] = energies[1]

    return [p1,p2, p3, p4]


def get_dca(df, point, theta_cone):
    '''
    Smarter get_dca method that uses dot products instead of point rotations
    to determine the DCA.

    This should be faster and more straightforward.
    '''
    a = df.loc[:,['x1','y1','z1','x2','y2','z2']].as_matrix()
    v2_arr = a[:, 3:] - a[:, :3]

    v1_arr = np.tile(point, a.shape[0]).reshape([a.shape[0], len(point)])
    v1_arr -= a[:,:3]

    v1_mag = np.linalg.norm(v1_arr, axis=1, keepdims=False)

    v1_arr /= np.linalg.norm(v1_arr, axis=1, keepdims=True)
    v2_arr /= np.linalg.norm(v2_arr, axis=1, keepdims=True)

    dot_arr = (v1_arr * v2_arr).sum(axis=1)
    theta = math.pi - np.arccos(dot_arr)

    y = v1_mag * np.cos(theta)
    radius = y * np.tan(theta_cone)
    distance_from_cone_axis = y * np.tan(theta)

    dca = np.abs((distance_from_cone_axis - radius) * np.cos(theta_cone))
    dca[np.isnan(dca)] = -10.0

    return dca


def get_calculated_fields(df):
    # type: (pandas.core.frame.DataFrame) -> Union[pandas.core.frame.DataFrame, None]

    df["E"] = df.E1 + df.E2
    df["cos_theta_M"] = 1.0 - 0.511*(1.0/df.E2 - 1.0/df.E)
    df.cos_theta_M[numpy.isinf(df.cos_theta_M)] = numpy.nan
    df["theta_M"] = np.arccos(df["cos_theta_M"])

    df["E2k"] = df.Ek - df.E1
    df["dE"] = df["E"] - df["Ek"]
    df["cos_theta_k"] = 1.0 - 0.511*(1.0/df.E2k - 1.0/df.Ek)
    df.cos_theta_k[numpy.isinf(df.cos_theta_k)] = numpy.nan
    df["theta_K"] = np.arccos(df["cos_theta_k"])
    df["dTheta"] = df["theta_K"] - df["theta_M"]

    df['dY'] = df.y2 - df.y1
    df['dS'] = numpy.sqrt( (df.x2 - df.x1)*(df.x2 - df.x1) + (df.z2 - df.z1)* (df.z2 - df.z1))

    return df


def get_calculated_fields_train(df):

    df = get_calculated_fields(df)
    origins = df[['X0', 'Y0', 'Z0']].as_matrix()
    df["DCA"] = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1','z1', 'E2', 'x2', 'y2', 'z2']].as_matrix(), origins)
    df.ix[df.DCA < 0,'DCA'] = BIG_DCA
    df["DCAk"] = energy_matcher.get_dca_for_events(df[['E1', 'x1', 'y1','z1', 'E2k', 'x2', 'y2', 'z2']].as_matrix(), origins)
    df.ix[df.DCAk < 0,'DCAk'] = BIG_DCA

    return df


def get_calculated_fields_test(df, dca_origin):

    # origins = numpy.tile(dca_origin, df.shape[0]).reshape([df.shape[0], 3])
    df = get_calculated_fields(df)
    df["DCA"] = get_dca(df, dca_origin, np.array(df.theta_M))
    df.ix[df.DCA < 0,'DCA'] = BIG_DCA
    df["DCAk"] = get_dca(df, dca_origin, np.array(df.theta_K))
    df.ix[df.DCAk < 0,'DCAk'] = BIG_DCA

    return df


def label_events_rules(df, energies, cross_sections, output_folder, origins, distance_param=10):
    '''
    Attempts to label the energy and scatter sequence of the events.
    :param event_permutations: A list of data frames, each with one of the possible combinations of event energies and scatter sequences
    :param output_folder:
    :param distance_param: Arbitrary value that is used for event weighting. 10 (1 cm) is probably OK.
    :return: dataframe with the best predicted event energies and scatter sequences.
    '''

    energies.sort()

    p1 = df.copy()
    p2 = df.copy()
    p2.ix[:, ['E1', 'x1', 'y1', 'z1', "E2", 'x2', 'y2', 'z2']] = p1.ix[:, ["E2", "x2", "y2", "z2", "E1", "x1", "y1", "z1"]].values
    E1 = numpy.array(p1.E1)
    E2 = numpy.array(p1.E2)
    e = numpy.array(energies)
    xs = numpy.array(cross_sections)
    p2["Ek"] = energy_matcher.get_energy_labels(numpy.array(p2.E1), numpy.array(p2.E2), energies, cross_sections)
    p1["Ek"] = energy_matcher.get_energy_labels(E1, E2, e, xs)


    p1['E2k'] = p1.Ek - p1.E1
    p2['E2k'] = p2.Ek - p2.E1

    #Remove scatters where both orderings are non-physical
    # and set both p1 and p2 to physical ordering in the case where only one is non-physical
    energy_and_compton_edge_allowed = numpy.logical_or(p1.Ek > 0.0, p2.Ek > 0.0)
    print "Before good event cut length, ", p1.shape[0]
    p1 = p1[energy_and_compton_edge_allowed]
    p2 = p2[energy_and_compton_edge_allowed]
    print "After good event cut length, ", p1.shape[0]

    #Lehner et al. IEEE Transactions Vol. 51, NO.4, Aug. 2004
    # at 1.17 and 1.33 MeV the first Compton interaction deposits more
    # energy > 80% of the time. Therefore, we use BigE ordering.
    mask = (p1.E1 < p2.E1) & (p2.Ek > 0.0)
    p1.ix[mask, :] = p2.ix[mask, :].values

    p1.ix[p1.Ek <= 0.0, :] = p2.ix[p1.Ek <= 0.0, :].values
    # p2.ix[p2.Ek <= 0.0, :] = p1.ix[p2.Ek <= 0.0, :].values

    p1["cos_theta_k"] = 1.0 - 0.511*(1.0/p1.E2k - 1.0/p1.Ek)
    p1.ix[p1.Ek < 0,'cos_theta_k'] = 2.0*math.pi



    # #keep event with scattering angle closest to Theta degrees
    # theta = 10
    # cos_theta = cos(theta*math.pi/180.0)
    # mask = numpy.abs(p1.cos_theta_k - cos_theta) < numpy.abs(p2.cos_theta_k - cos_theta)
    # p1.ix[mask, : ] = p2.ix[mask, : ].values

    p1["DCAk"] = energy_matcher.get_dca_for_events(p1[['E1', 'x1', 'y1','z1', 'E2k', 'x2', 'y2', 'z2']].as_matrix(), origins)
    assert(p1.DCAk.all > 0.0)
    # p1.ix[p1.DCAk < 0,'DCAk'] = BIG_DCA # set invalid events


    utilities.plot_1D(p1.DCAk, 100, "DCAk", "Counts", "DCAkrules", output_folder)
    utilities.plot_1D(p1.Ek, 100, "Ek", "Counts", "Ekrules", output_folder)

    return p1


def filter_events(df, verbose=False):

    min_scatter_energy = 0.06
    max_event_energy = 1.5
    min_event_energy = 1.0
    min_scatter_distance = 1.5

    # field_names = ['E1','x1','y1','z1','E2','x2','y2','z2', 'Ek']
    field_names = list(df)

    a = df.loc[:,field_names].as_matrix()

    n0 = a.shape[0]

    a = a[(a[:,0] > min_scatter_energy) & (a[:,4] > min_scatter_energy)]
    n1 = a.shape[0]

    a = a[(a[:,0] + a[:,4]) > min_event_energy]
    n2 = a.shape[0]

    a = a[(a[:,0] + a[:,4]) < max_event_energy]
    n3 = a.shape[0]

    #get the scatter distance
    ds2 = a[:,1:4] - a[:,5:8]
    ds2 *= ds2
    ds2 = ds2.sum(axis=1)

    a = a[ds2 > (min_scatter_distance * min_scatter_distance)]
    n4 = a.shape[0]

    if verbose:
        print "Initial events: %d" % (n0)
        print "Minimum scatter energy > %.1f: %d (%d)" % (min_scatter_energy, n1, n0 - n1)
        print "Minimum event energy > %.1f: %d (%d)" % (min_event_energy, n2, n1 - n2)
        print "Maximum event energy < %.1f: %d (%d)" % (max_event_energy, n3, n2 - n3)
        print "scatter distance > %.1f: %d (%d)" % (min_scatter_distance, n4, n3 - n4)

    df_new = pandas.DataFrame(data=a, columns=field_names)

    return df_new


def filter_bad_origins(df):
    '''
    Some of the MC events have origins outside of the defined creation volume. We can throw out these events
    rather than trying to fix the MC or figure out why it is generating such event.
    :param df: MC training data
    :return: MC data without the events with bad origins.
    '''
    df = df[(df[:, -5] < 30.1) & (df[:, -5] > -30.1)]
    df = df[(df[:, -4] < 1) & (df[:, -4] > -1)]
    df = df[(df[:, -3] < 130) & (df[:, -3] > -130)]

    return df


def produce_results_select_by_rules(df, output_folder, energies=[1.17, 1.33]):
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.style.use('ggplot')
    matplotlib.rcParams.update({'font.size': 14})
    plt.close()
    ax = plt.subplot(121)
    # plt1 = plot_E1_probabilities(output_folder, energies)

    xe = numpy.linspace(0.0, 1.15, 11500)
    # xe = numpy.linspace(0.0, 0.90, 9000)
    dE = 0.03
    p117 = numpy.array([energy_matcher.probability_E1(1.17, e, dE, 1000) for e in xe])
    p133 = numpy.array([energy_matcher.probability_E1(1.33, e, dE, 1000) for e in xe])

    plt.xlabel(r'Energy Deposited in $\mathregular{1^{st}}$ Scatter (MeV)')
    plt.ylabel(r'Free Electron Cross Section (barn)/%.0f (keV)' % (dE * 1000.0))
    plt.title('Relative Cross Sections')
    # plt.plot(xe, numpy.log2(p117))
    # plt.plot(xe, numpy.log2(p133))
    plt.plot(xe, p117, label=r'$\mathregular{1.17 (MeV)}$', linewidth="2.0")
    plt.plot(xe, p133, label=r'$\mathregular{1.33 (MeV)}$', linewidth="2.0")
    plt.grid(True)
    plt.tight_layout()

    legend = ax.legend(loc='upper left', shadow=False, fancybox=True, title=r'Gamma Energy')
    plt.setp(plt.gca().get_legend().get_title(), fontsize='14')
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')


    ax = plt.subplot(122)
    df.groupby(df.Ek).E1.plot.kde(legend=True, xlim=[0.0,1.33])
    plt.xlabel(r'Energy Deposited in $\mathregular{1^{st}}$ Scatter (MeV)')
    plt.ylabel(r'Number of Events' )
    plt.title('Energy Label Distributions')
    plt.grid(True)
    plt.tight_layout()

    legend = ax.legend(loc='upper left', shadow=False, fancybox=True, title=r'Gamma Energy')
    plt.setp(plt.gca().get_legend().get_title(), fontsize='14')
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')

    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)
    plt.savefig("%s/GammaEnergyLabelsByE1.png" % output_folder)
    # plt.show()

    return


def plot_E1_probabilities(output_folder, energies, plt=None):
    if plt == None:
        import matplotlib.pyplot as plt
    ax = plt.subplot()
    # import matplotlib
    # import matplotlib.pyplot as plt
    # matplotlib.style.use('ggplot')
    # matplotlib.rcParams.update({'font.size': 14})
    # plt.close()
    # ax = plt.subplot(121)

    xe = numpy.linspace(0.0, 1.15, 11500)
    # xe = numpy.linspace(0.0, 0.90, 9000)
    dE = 0.03
    p117 = numpy.array([energy_matcher.probability_E1(1.17, e, dE, 1000) for e in xe])
    p133 = numpy.array([energy_matcher.probability_E1(1.33, e, dE, 1000) for e in xe])

    plt.xlabel(r'Energy Deposited in $\mathregular{1^{st}}$ Scatter (MeV)')
    plt.ylabel(r'Free Electron Cross Section (barn)/%.0f (keV)' % (dE * 1000.0))
    plt.title('Relative Cross Sections')
    # plt.plot(xe, numpy.log2(p117))
    # plt.plot(xe, numpy.log2(p133))
    plt.plot(xe, p117, label=r'$\mathregular{1.17 (MeV)}$', linewidth="2.0")
    plt.plot(xe, p133, label=r'$\mathregular{1.33 (MeV)}$', linewidth="2.0")
    plt.grid(True)
    plt.tight_layout()

    legend = ax.legend(loc='upper left', shadow=False, fancybox=True, title=r'Gamma Energy')
    plt.setp(plt.gca().get_legend().get_title(), fontsize='14')
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
    # plt.savefig("%s/relative_probabilities.png" % output_folder)
    # plt.show()

    return plt


#------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------

def usage():
    print "USAGE: %s [training event file] [training origin] [test event file] [test origin] [output folder]"  % (sys.argv[0])
    print "\nExample:\n %s CC_events1.csv 0.0,0.0,0.0  CC_events2.csv 0.0,0.0,6.0  ~/public_html/em" % (sys.argv[0])
    sys.exit(-1)


def main():
    ENERGIES = numpy.array([1.17, 1.33])
    CROSS_SECTIONS =  numpy.array([0.5, 0.5])

    if len(sys.argv) == 5:
        file_train = sys.argv[1]
        file_test = sys.argv[2]
        print "Processing arg[3]: ", sys.argv

        origin = [float(x) for x in sys.argv[3].split(",")]
    else:
        print "ERROR: Expected 5 parameters. Got %d . . ." % len(sys.argv)
        print "arguments: ", sys.argv
        usage()

    output_folder = sys.argv[4]

    # plot_E1_probabilities(output_folder, ENERGIES)
    df_test = read_csv_events(file_test)
    df_test = filter_events(df_test)
    origins = numpy.tile(origin, df_test.shape[0]).reshape([df_test.shape[0], 3])

    df_labeled_rules = label_events_rules(df_test, ENERGIES, CROSS_SECTIONS, output_folder, origins)
    # produce_results_select_by_rules(df_labeled_rules, output_folder, energies=[1.17, 1.33])

    #df_labled_ideal.sort_values(by=['DCA'], ascending=True, inplace=True)
    outputfile = file_test.replace(".csv", "_labelled_RULES.csv")
    df_labeled_rules.to_csv(outputfile,
                       columns=['E1', 'x1', 'y1', 'z1', 'E2', 'x2', 'y2',
                                'z2'], float_format="%0.5f", index=False,
                            header=False)

    outputfile = file_test.replace(".csv", "_labelled_KNOWN.csv")
    df_labeled_rules.to_csv(outputfile,
            columns=['E1', 'x1', 'y1', 'z1', 'E2k', 'x2', 'y2','z2'],
            float_format="%0.5f", index=False, header=False)
    outputfile = file_test.replace(".csv", "_labelled_SMALLE.csv")
    df_labeled_rules.to_csv(outputfile,
            columns=['E2', 'x2', 'y2','z2','E1', 'x1', 'y1', 'z1'],
            float_format="%0.5f", index=False, header=False)

    event_permutations = get_event_permutations(df_test)
    event_permutations = [get_calculated_fields_test(X, origin) for X in event_permutations]

    df_labled_ideal = label_events_ideal(event_permutations, output_folder)
    df_labled_ideal.sort_values(by=['DCA'], ascending=True, inplace=True)
    outputfile = file_test.replace(".csv", "_labelled_IDEAL.csv")
    df_labled_ideal.to_csv(outputfile,
                       columns=['E1', 'x1', 'y1', 'z1', 'E2k', 'x2', 'y2',
                                'z2'], index=False, header=False)

    print "Reading training file %s . . ." % sys.argv[1]
    train = read_training_events(file_train)
    train = filter_bad_origins(train)
    train = energy_matcher.get_dca_for_training_events(train)
    df_train = pandas.DataFrame(numpy.array(train), columns=energy_matcher.column_names_train_enum().keys())

    print energy_matcher.column_names_train_enum().keys()
    df_train = filter_events(df_train)
    # df_bad_events = make_bad_events_sample(df_train)
    # df_train = make_training_sample_weirdlyGood(df_train, df_bad_events)
    df_train = get_calculated_fields_train(df_train)
    df_train = make_training_sample(df_train)

    print df_train[df_train.DCAk > 800].head()

    random_forest_file_name = "%s/%s" % (output_folder, os.path.basename(file_train).replace(".txt","") + "_RF.dat")
    randomforest = train_random_forest(df_train, output_folder, random_forest_file_name, rebuild_forest = False)
    # make_efficiency_plots(df_train, output_folder)

    df_selected = label_events_random_forest(randomforest, event_permutations, output_folder)
    print "Sorting by PDCA . . ."
    df_selected.sort_values(by=['PDCA'], ascending=True, inplace=True)
    outputfile = file_test.replace(".csv", "_labelled_PDCA.csv")
    df_selected.to_csv(outputfile, columns=['E1', 'x1', 'y1', 'z1', 'E2k', 'x2', 'y2', 'z2'], index=False, header=False)
    df_selected.sort_values(by=['dE'], ascending=True, inplace=True)
    outputfile = file_test.replace(".csv", "_labelled_dE.csv")
    df_selected.to_csv(outputfile, columns=['E1', 'x1', 'y1', 'z1', 'E2k', 'x2', 'y2', 'z2'], index=False, header=False)
    df_selected.sort_values(by=['dTheta'], ascending=True, inplace=True)
    outputfile = file_test.replace(".csv", "_labelled_dTheta.csv")
    df_selected.to_csv(outputfile, columns=['E1', 'x1', 'y1', 'z1', 'E2k', 'x2', 'y2', 'z2'], index=False, header=False)
    # print "Correlation: ", df_selected[['dTheta', 'dE', 'PDCA', 'DCAk', 'DCA']].corr(method='spearman')

    sys.exit()
    df_test = get_calculated_fields_test(df_test, origin)
    df_test = put_scatters_in_order(df_test)
    matrix_train = df_train.as_matrix()
    test_events = filter_events(df_test.as_matrix())
    df_test = energy_matcher.run_known_origin_algorithm(test_events, np.array(otest), np.array([1.17, 1.33]), np.array([0.5, 0.5]), 0.1)

    print "Applying random forest . . ."
    df_trained.to_csv("%s/trained_all.csv" % output_folder)


    print "Saving em . . ."
    df_trained.to_csv("%s/trained_em.csv" % output_folder, columns=['E1', 'x1', 'y1', 'z1', 'E2', 'x2', 'y2', 'z2'], index=False)

    print "Saving Ek . . ."
    df_trained['E2Corr'] = df_trained['Ek'] - df_trained['E1']
    df_trained.to_csv("%s/trained_ek.csv" % output_folder, columns=['E1', 'x1', 'y1', 'z1', 'E2Corr', 'x2', 'y2', 'z2'], \
                      float_format="%0.4f", index=False)

    return


if __name__ == "__main__":

    # cProfile.run("main()", sort=2)
    main()





