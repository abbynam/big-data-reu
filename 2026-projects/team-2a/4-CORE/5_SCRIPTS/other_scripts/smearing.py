import sys
import numpy as np
import collections
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import math
import scipy

import pandas

sys.path.insert(0,'/home/dsmackin/projects/CORE/scripts')

import pyximport; pyximport.install()
import energy_matcher
import utilities
import smearing
import event_selector


#----------- FUNCTIONS ---------------------------------------
def combine_approx_equal_values(vals, cutoff=0.5):
    vals = np.array(vals)
    diff = vals[1:] - vals[:-1]
    mask = (diff > cutoff)
    mask = np.insert(mask, 0, True)

    return vals[mask]


def get_discrete_values(vals, counts_needed=100):
    '''
    The Polaris CC system uses discrete positions for single pixel events. For this
    approximation we are ignoring the continuous positions produced by averaging multiple
    pixels. This function counts the number of times each position in vals occurs. If the value
    occurs more than counts_needed times, then the poisition is added to the list of discrete
    position values.
    '''
    vd = {}
    for v in vals:
        vd[v] = vd.get(v, 0) + 1
    v_discrete = dict((k, v) for k, v in vd.items() if v >= counts_needed)

    d = v_discrete.keys()
    d.sort()

    return combine_approx_equal_values(d, cutoff=0.5)


# DISCRETE_VALUES_X = get_discrete_values(df_arr[0].x2, 100)
# DISCRETE_VALUES_Z = get_discrete_values(df_arr[0].z2, 100)

# These values were extracted using get_discrete_values and then saved here
DISCRETE_VALUES_X = np.array([-20.1, -18.38, -16.66, -14.94, -13.22, -11.5, -9.78, -8.06,
                              -6.34, -4.62, -2.9, 2.9, 4.62, 6.34, 8.06, 9.78,
                              11.5, 13.22, 14.94, 16.66, 18.38, 20.09])

DISCRETE_VALUES_Z = np.array([-117.91, -116.19, -114.47, -112.75, -111.03, -109.31, -107.59,
                              -105.87, -104.15, -102.43, -100.71, -94.9, -93.19, -91.47,
                              -89.75, -88.03, -86.31, -84.59, -82.87, -81.15, -79.43,
                              -77.71, -53.14, -51.42, -49.7, -47.98, -46.26, -44.54,
                              -42.82, -41.1, -39.38, -37.66, -35.94, -30.13, -28.42,
                              -26.7, -24.98, -23.26, -21.54, -19.82, -18.1, -16.38,
                              -14.66, -12.94, 11.63, 13.35, 15.07, 16.79, 18.51,
                              20.23, 21.95, 23.67, 25.39, 27.11, 28.83, 34.64,
                              36.35, 38.07, 39.79, 41.51, 43.23, 44.95, 46.67,
                              48.39, 50.11, 51.83, 76.4, 78.12, 79.84, 81.56,
                              83.28, 85., 86.72, 88.44, 90.16, 91.88, 93.6,
                              99.4, 101.12, 102.84, 104.56, 106.28, 108., 109.72,
                              111.44, 113.16, 114.88, 116.6])

def smear_laterial(vals, h):
    vals += np.random.normal(size=len(vals), scale=h, loc=0.0)
    return vals

def smear_longitudinal(vals, h):
    y_min = min(vals)
    y_max = max(vals)
    vals += np.random.normal(size=len(vals), scale=h, loc=0.0)
    vals = np.minimum(vals, y_max)
    vals = np.maximum(vals, y_min)

    return vals

def discretize_values(vals, dvals):
    m = len(vals)
    n = len(dvals)
    d_mat = np.tile(dvals,m).reshape(m,n)
    v_mat = np.tile(vals, n).reshape(m,n, order='F')

    diff = np.abs(d_mat - v_mat)
    indices = np.argmin(diff, axis=1)
    vals = np.take(dvals, indices)

    return vals


def get_smeared_positions(df, h_lateral, h_longitudinal):
    df.x1 = discretize_values(smear_laterial(df.x1, h_lateral), DISCRETE_VALUES_X)
    df.x2 = discretize_values(smear_laterial(df.x2, h_lateral), DISCRETE_VALUES_X)
    df.z1 = discretize_values(smear_laterial(df.z1, h_lateral), DISCRETE_VALUES_Z)
    df.z2 = discretize_values(smear_laterial(df.z2, h_lateral), DISCRETE_VALUES_Z)

    df.y1 = smear_longitudinal(df.y1, h_longitudinal)
    df.y2 = smear_longitudinal(df.y2, h_longitudinal)

    return df


def get_smeared_energy(energy, g1, g2, h1, h2, h3):
    # For explanation of energy uncertainty, see Mackin et al.
    # "The effects of Doppler broadening and detector resolution on the performance of three-stage Compton cameras"

    # Note: the smearing from the paper does not match data, therefore we developed
    # an empirical method.

    epsilon = 1E-10
    N = len(energy)
    if g1 + g2 >= 1.0:
        denom = (g1 + g2 + epsilon)
        g1 = g1 / denom
        g2 = g2 / denom
    N_gauss = int(g1 * N)
    N_uniform = int(g2 * N)
    N_exp = N - N_gauss - N_uniform

    x1 = np.random.normal(size=N_gauss, scale=h1, loc=0.0)
    x2 = np.random.laplace(size=N_exp, scale=h2, loc=0.0)
    #     x2 = np.random.standard_cauchy(size=N_cauchy) * h2

    x3 = h3 * np.random.uniform(low=-1.3, high=0.7, size=N_uniform)
    x = np.concatenate((x1, x2, x3))
    np.random.shuffle(x)

    e_smeared = energy + x * 0.1

    return e_smeared


def smear_MC(df0, g1, g2, h1, h2, h3, h_lateral, h_longitudinal):
    df = df0.copy()
    df = get_smeared_positions(df, h_lateral, h_longitudinal)
    df.E1 = get_smeared_energy(df.E1, g1, g2, h1, h2, h3)
    df.E2 = get_smeared_energy(df.E2, g1, g2, h1, h2, h3)

    df['E'] = df.E1 + df.E2

    return df


def cost_function(n, n0):
    '''Difference metric for n and n0 based on ChiSq. Intended to
        be used for optimization problems, initially for determining
        the best smearing parameters to make Monte Carlo data more
        characteristic of CC data.

    Args:
        n (numpy arr float): bin values for observed data normalized to 1.
        n0 (numpy arr float): bin values for expected data normalized to 1.

    Yields:
        float: cost

    '''

    numerator = (n - n0)
    numerator *= numerator

    mask = n0 > 0.0
    cost = (numerator[mask]/n0[mask]).sum()

    return cost


def get_dca(df, point):
    '''
    Smarter get_dca method that uses dot products instead of point rotations
    to determine the DCA.

    This should be faster and more straightforward.
    '''

    v1_arr = np.tile(point, len(df)).reshape([len(df), len(point)])
    v2_arr = np.transpose(np.array([df.x2, df.y2, df.z2]))
    scat1 = np.transpose(np.array([df.x1, df.y1, df.z1]))

    v1_arr -= scat1
    v2_arr -= scat1

    v1_mag = np.linalg.norm(v1_arr, axis=1, keepdims=False)

    v1_arr = v1_arr / np.linalg.norm(v1_arr, axis=1, keepdims=True)
    v2_arr = v2_arr / np.linalg.norm(v2_arr, axis=1, keepdims=True)

    dot_arr = v1_arr[:, 0] * v2_arr[:, 0] + v1_arr[:, 1] * v2_arr[:, 1] + v1_arr[:, 2] * v2_arr[:, 2]
    theta = math.pi - np.arccos(dot_arr)

    y = v1_mag * np.cos(theta)
    radius = y * np.tan(df.theta_K)
    distance_from_cone_axis = y * np.tan(theta)

    df['DCAnew'] = np.abs((distance_from_cone_axis - radius) * np.cos(df.theta_K))

    return df


def get_smears_dict(df, center, bin_width=0.05, fields=['E1c', "dE1"]):
    df = calculate_mismeasurements(df.copy(), center)
    events = np.array(df[fields])

    energies_bin_num = np.round(events[:, 0] / bin_width)
    bin_energies = np.unique(energies_bin_num) * bin_width
    nbins = int(max(energies_bin_num + 1))

    uncertainties = [np.append(events[energies_bin_num == i, 1], [0.0]) for i in xrange(nbins)]

    return collections.OrderedDict(zip((bin_energies * 1000).astype(int), uncertainties))


def get_smears_dict_E2(df, center, bin_width=0.05):
    return get_smears_dict(df, center, bin_width=0.05, fields=['E2c', "dE2"])


def get_smear(x, smears_dict, scalar=0.5):
    '''
    Produce a smear factor for the first scatter energy based on the smearing from CC data events.

    CC_train is a 2D vector of 1st scatter energies and measured energy errors calculated by determining the energy needed to
    make the DCA for the cone equal zero.
    '''

    bin_energies = np.array(smears_dict.keys()) * 0.001
    assert (len(bin_energies) > 1)

    bin_width = (bin_energies[1] - bin_energies[0])
    energies_bin = np.round(x / bin_width) * bin_width
    #     print np.unique(energies_bin)
    energies_bin[energies_bin > max(bin_energies)] = max(bin_energies)
    #     smear_values = [np.random.choice(smears_dict[int(x*1000)]) for x in energies_bin]

    runif = np.random.random(len(energies_bin))
    rgauss = np.random.normal(loc=0.0, scale=0.05, size=len(energies_bin))
    smear_values = energies_bin * 0.0
    smear_dict_index = (energies_bin * 1000 + 0.5).astype(int)

    ## @TODO iterate over the dictionary array options and use slicing instead
    no_smear_count = 0
    for i, x in enumerate(energies_bin):
        # if(i % 10000 == 0): print "Smearing", i, x
        try:
            smears = smears_dict[smear_dict_index[i]]
            index = int(runif[i] * len(smears))
            smear_values[i] = smears[index]
        except KeyError:
            # print "No smear for", i, x
            smear_values[i] = rgauss[i] * x
            no_smear_count += 1
            continue

    # print no_smear_count, "of", len(energies_bin), "scatter were randomly " \
    #                                                "smeared . . ."
    return smear_values * scalar


def calculate_mismeasurements(df, point):
    '''Calculates the amount of energy needed by E1 to make the DCA go to 0. Both
        the corrected energy, E1c, and the energy added dE1 are added to the
        dataframe.
    '''

    v1_arr = np.tile(point, len(df)).reshape([len(df), len(point)])
    v2_arr = np.transpose(np.array([df.x2, df.y2, df.z2]))
    scat1 = np.transpose(np.array([df.x1, df.y1, df.z1]))

    v1_arr -= scat1
    v2_arr -= scat1

    v1_arr = v1_arr / np.linalg.norm(v1_arr, axis=1, keepdims=True)
    v2_arr = v2_arr / np.linalg.norm(v2_arr, axis=1, keepdims=True)

    dot_arr = v1_arr[:, 0] * v2_arr[:, 0] + \
              v1_arr[:, 1] * v2_arr[:, 1] + \
              v1_arr[:, 2] * v2_arr[:, 2]
    cos_theta = np.cos(math.pi - np.arccos(dot_arr))

    e1c = np.array(df.Ek*df.Ek*(1.0 - cos_theta)/(df.Ek*(1.0 - cos_theta)+0.511))

    df["E1c"] = e1c
    df['E2c'] = df.Ek - df.E1c
    df['dE1'] = df.E1 - df.E1c
    df['dE2'] = df.E2 - df.E2c
    df['thetac'] = np.arccos(cos_theta)
    df['dThetac'] = df.thetac - df.theta_K

    return df


def get_dE1(df, point):
    '''Calculates
    Calculates the amount of energy needed by E1 to make the DCA go to 0.

    Both the corrected energy, E1c, and the energy added dE1 are added to the dataframe.
    '''

    v1_arr = np.tile(point, len(df)).reshape([len(df), len(point)])
    v2_arr = np.transpose(np.array([df.x2, df.y2, df.z2]))
    scat1 = np.transpose(np.array([df.x1, df.y1, df.z1]))

    v1_arr -= scat1
    v2_arr -= scat1

    v1_arr = v1_arr / np.linalg.norm(v1_arr, axis=1, keepdims=True)
    v2_arr = v2_arr / np.linalg.norm(v2_arr, axis=1, keepdims=True)

    dot_arr = v1_arr[:, 0] * v2_arr[:, 0] + v1_arr[:, 1] * v2_arr[:, 1] + v1_arr[:, 2] * v2_arr[:, 2]
    cos_theta = np.cos(math.pi - np.arccos(dot_arr))

    e1c = np.array(df.Ek * df.Ek * (1.0 - cos_theta) / (df.Ek * (1.0 - cos_theta) + 0.511))

    df["E1c"] = e1c
    df['E2c'] = df.Ek - df.E1c
    df['dE1'] = df.E1 - df.E1c
    df['dE2'] = df.E2 - df.E2c
    df['thetac'] = np.arccos(cos_theta)
    df['dThetac'] = df.thetac - df.theta_K

    return df


def bin_data_normed(x, num_bins, x_range):
    assert(len(x_range) == 2)
    assert(x_range[1] - x_range[0] > 0)

    N = len(x)
    range_min = x_range[0] + 0.0
    range_max = x_range[1] * 1.0001 #make bin edge bigger than largest x.
    range_length = range_max - range_min + 1.0E-10

    bins = np.linspace(range_min, range_max, num_bins, endpoint=True)
    x_vals = bins * range_length / num_bins + range_min

    step_factor = float(num_bins)/range_length

    xshifted = x - range_min
    hist = np.zeros(num_bins, dtype='float')

    bin_nums = np.int16(xshifted * step_factor)
    bin_nums[bin_nums == num_bins] += -1
    bin_nums = bin_nums[(bin_nums > 0.0) & (bin_nums < num_bins)]

    def f(n):
        hist[n] += 1
    [f(bin_num) for bin_num in bin_nums]
#
#     for i in range(len(xshifted)):
#         bin_num = int(xshifted[i] * step_factor)
#         if bin_num < 0.0 or bin_num > num_bins: #
# #             print "\n\nError filling hist, invalid bin number\n---------\n", i, bin_num, xshifted[i], step_factor, \
# #                                     range_length, range_min, range_max, num_bins
#             continue
#         if bin_num == num_bins: bin_num += -1
#
#         hist[bin_num] += 1

    hist = 100.0*hist/hist.max()

#     print "LENGTHS:", len(bins), len(x_vals), len(hist)
    return bins, x_vals, hist


def label_ideal(df, origin):
    df = event_selector.filter_events(df)
    df['E'] = df.E1 + df.E2
    event_permutations = event_selector.get_event_permutations(df)
    event_permutations = [event_selector.get_calculated_fields_test(X, origin) for X in event_permutations]
    df_labelled_ideal = event_selector.label_events_ideal(event_permutations)
    df_labelled_ideal = df_labelled_ideal.sort_values(by=['DCAk'], ascending=True)

    return df_labelled_ideal