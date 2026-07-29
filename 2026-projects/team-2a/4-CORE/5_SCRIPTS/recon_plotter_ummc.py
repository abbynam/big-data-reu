#! /usr/bin/python
#  revision date: 13 August 2018

import sys, os
import numpy
import scipy
import scipy.interpolate
import scipy.optimize as opt
import matplotlib
import matplotlib.pyplot as plt
#import multiprocessing   # commented out multiprocessing, crashes in OSX

import Intensity3D


USAGE = '''
USAGE: %s [Path to 3D dose file] [output path]

The following example generates plots using the image intensity information from the file output.dat.

   python %s output.dat outputfolder
'''

#------------------------------------------------------------------
# SETTINGS / PARAMETERS
#------------------------------------------------------------------

### Expected Source Positions ###

#==== 190805 - ALL RUNS (Na-22)  ====#
#EXPECTED_SOURCE_LOCATION = [-1.0, 0.0, 11.0]  # all runs

#==== 220615 - ALL RUNS (66 MeV proton beam)  ====#
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 0.0]  # change to 0, 0, 0 since there is a -11 in x translation during coordinate_transformation

#==== 240502 - ALL RUNS (Co-60 test run)  ====#
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 0.0]  # all runs

#==== 250121 - ALL RUNS (Co-60 test run)  ====#
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 0.0]  # run 1 & 4
#EXPECTED_SOURCE_LOCATION = [12.0, 12.0, 0.0]  # run 2

#==== 250127 - ALL RUNS (Co-60 rotating source run)  ====#
#EXPECTED_SOURCE_LOCATION = [20.0, 11.0, 5.0]  # runs 1 & 2  (source was stationary)
#EXPECTED_SOURCE_LOCATION = [-22.0, 7.0, 5.0]  # run 3  <-- based on results from reconstructions
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 5.0]  # run 4  <-- based on results from reconstructions
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 5.0]  # run 5  <-- based on results from reconstructions
#EXPECTED_SOURCE_LOCATION = [14.0, -18.0, 5.0]  # run 7  (source was stationary)

#==== UMB - Three RUNS (Cs-137 test run)  ====#
#EXPECTED_SOURCE_LOCATION = [150.0, 0.0, 0.0]  # left
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 0.0]  # center
#EXPECTED_SOURCE_LOCATION = [-150.0, 0.0, 0.0]  # right

#==== UMB - 250725 (Cs-137 test run)  ====#
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 0.0]  # only one source position

#==== UMMC - 250828 (Cs-137 point source measurements)  ====#
#EXPECTED_SOURCE_LOCATION = [120.0, 0.0, 0.0]     # run 1
#EXPECTED_SOURCE_LOCATION = [270.0, 0.0, 0.0]     # run 2
#EXPECTED_SOURCE_LOCATION = [120.0, 0.0, -141.0]  # run 3
#EXPECTED_SOURCE_LOCATION = [120.0, 0.0,  140.0]  # run 4

#==== UMMC - 250926 (Cs-137 point source measurements)  ====#
#EXPECTED_SOURCE_LOCATION = [120.0, 0.0, 0.0]   # run 1 & 4
#EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 0.0]     # run 2 & 5
#EXPECTED_SOURCE_LOCATION = [-130.0, 0.0, 0.0]  # run 3 & 6

#==== MPTC - 251206 (proton beam measurements)  ====#
EXPECTED_SOURCE_LOCATION = [0.0, 0.0, 0.0]   # run 18 & 28


#------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------

def plot_generic_planes(obj3D, num_planes_per_dimension = 5):
    ranges = obj3D.get_ranges()

    assert(len(ranges) == 3)
    for i,r in enumerate(ranges):
        step = (r[1] - r[0])/(num_planes_per_dimension + 1)
        depths = numpy.linspace(r[0] + 0.5*step, r[1] - 0.5*step, num_planes_per_dimension)
        for depth in depths:
            plot = obj3D.get_plot(i, depth)


def plot_y0_profiles(obj3D, output_folder):
    try:
        plane = obj3D.getIntensityPlane(1, 0.0)
    except ValueError as e:
        # Skip this plot if the requested plane is out of bounds and continue
        print(f"[plot_y0_profiles] Skipping y=0 profile plot: {e}")
        return None

    k,i = numpy.unravel_index(plane.argmax(), plane.shape)

    x, tmp, z = obj3D.get_coordinates(i, 0, k)

    plot_profiles(obj3D, output_folder, dimension=0, fixed_coordinates = [0.0, z])
    plot_profiles(obj3D, output_folder, dimension=2, fixed_coordinates = [x, 0.0])


def plot_profiles(obj3D, output_folder, dimension=2, fixed_coordinates = [0.0, 0.0], normalize=False):  # change default value of normalize from True to False [swp - 231101]
    try:
        points, values = obj3D.get_profile(dimension, fixed_coordinates)
    except ValueError as ve:
        # Provide more context about which profile request failed
        print(f"[plot_profiles] ValueError when getting profile for dimension={dimension}, fixed_coordinates={fixed_coordinates}: {ve}, skipping plot.")
        return None

    if normalize == True:
        values = values/max(values)

    dimension_labels = ["x", "y", "z"]

    plt.clf()

    ax = plt.subplot(111)
    plt.plot(points, values, label="0 mm", lw=1)

    #plt.title("150 MeV Proton Beam", fontsize=11)
    legend = ax.legend(loc='lower left', title=r'Shift [mm]', fontsize=6)
    ax.grid(True,linestyle='-',which='both',color='0.750')

    plt.setp(ax.get_xticklabels(), rotation=30, fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    plt.setp(legend.get_title(), fontsize=6)
    plt.xlabel("%s [mm]" % dimension_labels[dimension], fontsize=10)

    plt.ylabel("Gamma Emission", fontsize=10)
    plt.gcf().set_size_inches(5,3.3)
    plt.gcf().subplots_adjust(bottom=0.15)

    plot_name = "%s/profile_%s_%.1f_%.1f.png" % (output_folder, dimension_labels[dimension], fixed_coordinates[0], fixed_coordinates[1])
    plt.savefig(plot_name)

    return values


def plot_profile_expected(obj3D, output_folder, dimension = 2, fixed_coordinates = [0.0, 0.0, 0.0], actual_coordinates = [0.0, 0.0, 0.0], normalize = False):   # change default value of normalize from True to False [swp - 231101]

    other_dims = numpy.delete ( ([0,1,2]) , dimension ) # create list of other dimensions
    other_coords = [ fixed_coordinates[ other_dims[0] ], fixed_coordinates[ other_dims[1] ] ]

    try:
        points, values = obj3D.get_profile(dimension, other_coords)
    except ValueError as ve:
        print(f"[plot_profile_expected] ValueError when getting profile for dimension={dimension}, other_coords={other_coords}: {ve}")
        return None

    if normalize == True:
        values = values/max(values)

    dimension_labels = ["x", "y", "z"]

    # make plot
    plt.clf()
    ax = plt.subplot(111)
    # plot 1D profile points for given coordinates
    plt.plot(points, values, color = 'k', label = "Estimated = {:.2f} mm".format(fixed_coordinates[dimension]))
    plt.title('{}-profile at {} = {:.2f} mm and {} = {:.2f} mm'.format( dimension_labels[dimension], dimension_labels[ other_dims[0] ], fixed_coordinates[ other_dims[0] ] , dimension_labels[ other_dims[1] ], fixed_coordinates[ other_dims[1] ] ), fontsize = 10)
    plt.xlabel("%s direction (mm)" % dimension_labels[dimension], fontsize = 10)
    plt.ylabel("Gamma Emission", fontsize = 10)
    # add line for actual value
    ax.axvline(actual_coordinates[dimension], color='b', linestyle = '--', label = "Actual = {:.2f} mm".format(actual_coordinates[dimension])) # max position in data
    # make the plot pretty (size, grid, axis, legend, tc...)
    plt.gcf().set_size_inches(6, 6)
    plt.gcf().subplots_adjust(bottom=0.15)
    ax.grid(True, linestyle = '-', which = 'both', color = '0.750')
    plt.setp(ax.get_xticklabels(), fontsize = 8)
    plt.setp(ax.get_yticklabels(), fontsize = 8)
    ax.set_ylim(-0.04, 1.04)
    plt.legend(loc='upper left', fontsize = 8)
    plt.tight_layout()
    # save plot to file
    plot_name = "%s/profile_%s_%.1f_%.1f_maximum.png" % (output_folder, dimension_labels[dimension], fixed_coordinates[ other_dims[0] ], fixed_coordinates[ other_dims[1] ])
    plt.savefig(plot_name)

    return values


def plot_profile_gaussian(obj3D, output_folder, dimension = 2, fixed_coordinates = [0.0, 0.0, 0.0], actual_coordinates = [0.0, 0.0, 0.0], normalize = True):   # change default value of normalize from True to False [swp - 231101], then changed it back

    other_dims = numpy.delete ( ([0,1,2]) , dimension ) # create list of other dimensions
    other_coords = [ fixed_coordinates[ other_dims[0] ], fixed_coordinates[ other_dims[1] ] ]

    try:
        points, values = obj3D.get_profile(dimension, other_coords)
    except ValueError as ve:
        print(f"[plot_profile_gaussian] ValueError when getting profile for dimension={dimension}, other_coords={other_coords}: {ve}")
        return None

    if normalize == True:
        values = values/max(values)

    dimension_labels = ["x", "y", "z"]

    #  Create ability to kill runaway y-offsets by checking if y-offset way outside of range of values (try less than -1)
    #   - added 3-parameter gaussian fit (four parameters gave big problems for wide octane recons -> y-offsets of -560 mm)
    # one-line gaussian functions (four parameters: p[0]==amplitude, p[1]==mean, p[2]==stdev, p[3]==y-offset)
    gauss_fit = lambda p, x: p[0] * (1.0 / numpy.sqrt(2 * numpy.pi * (p[2]**2))) * numpy.exp(-(x - p[1])**2 / (2 * p[2]**2)) + p[3] # 1D Gaussian func
    e_gauss_fit = lambda p, x, y: (gauss_fit(p, x) - y) # 1D Gaussian fit
    # gaussian fit
    v0 = [0.0, fixed_coordinates[dimension], 1.0, 0.0] # initial guesses for Gaussian Fit.
    out = opt.leastsq(e_gauss_fit, v0[:], args = (points, values), maxfev = 100000, full_output = 1) # Gauss Fit

    #  Check y-offset value (out[0][4])
    if out[0][3] < -1.0:
        print ('    !! Warning !! / y-offset < -1.0 / value = {} / Using 3-parameter fit (fixing y-offset to zero)'.format(out[0][3]))
        #  use 3-parameter gaussian fit (four parameters gave big problems for wide octane recons -> y-offsets of -560 mm)
        # one-line gaussian functions (three parameters: p[0]==amplitude, p[1]==mean, p[2]==stdev)
        gauss_fit = lambda p, x: p[0] * (1.0 / numpy.sqrt(2 * numpy.pi * (p[2]**2))) * numpy.exp(-(x - p[1])**2 / (2 * p[2]**2)) # 1D Gaussian func
        e_gauss_fit = lambda p, x, y: (gauss_fit(p, x) - y) # 1D Gaussian fit
        # gaussian fit
        v0 = [0.0, fixed_coordinates[dimension], 1.0] # initial guesses for Gaussian Fit.
        out = opt.leastsq(e_gauss_fit, v0[:], args = (points, values), maxfev = 100000, full_output = 1) # Gauss Fit

    # Store gaussian output (v==fit parameters, covar==covariance matrix)
    v = out[0];  covar = out[1]
    # creates data points to plot gaussian fit
    points_gauss = numpy.arange(min(points), max(points), points[1] - points[0])
    values_gauss = gauss_fit(v, points_gauss) # runs the gaussian fit
    # make plot
    plt.clf()
    ax = plt.subplot(111)
    # plot 1D profile points for given coordinates
    plt.plot(points, values, color = 'k', label = "Estimated = {:.2f} mm".format(fixed_coordinates[dimension]))
    plt.title('{}-profile at {} = {:.2f} mm and {} = {:.2f} mm'.format( dimension_labels[dimension], dimension_labels[ other_dims[0] ], fixed_coordinates[ other_dims[0] ] , dimension_labels[ other_dims[1] ], fixed_coordinates[ other_dims[1] ] ), fontsize = 10 )
    plt.xlabel("%s direction (mm)" % dimension_labels[dimension], fontsize = 10)
    plt.ylabel("Gamma Emission", fontsize = 10)
    # add line for actual value
    ax.axvline(actual_coordinates[dimension], color='b', linestyle = '--', label = "Actual = {:.2f} mm".format(actual_coordinates[dimension])) # max position in data
    # plot gaussian data
    gaussian_peak_pos = points_gauss[ numpy.where(values_gauss == numpy.max(values_gauss))[0]][0]
    ax.plot(points_gauss, values_gauss, 'r-', label = 'Gaussian = {:.2f} mm'.format(gaussian_peak_pos), lw = 1) #fitted spectrum
    # plot FWHM
    FWHM = 2 * numpy.sqrt(2 * numpy.log(2)) * v[2]
    ax.axvspan( v[1] - FWHM/2, v[1] + FWHM/2, facecolor='g', alpha=0.2, label = 'FWHM = {:.2f} mm'.format(FWHM) )
    # make the plot pretty (size, grid, axis, legend, tc...)
    plt.gcf().set_size_inches(6, 6)
    plt.gcf().subplots_adjust(bottom=0.15)
    ax.grid(True, linestyle = '-', which = 'both', color = '0.750')
    plt.setp(ax.get_xticklabels(), fontsize = 8)
    plt.setp(ax.get_yticklabels(), fontsize = 8)
    ax.set_ylim(-0.04, 1.04)
    plt.legend(loc='upper left', fontsize = 8)
    plt.tight_layout()
    # save plot to file
    plot_name = "%s/profile_%s_%.1f_%.1f_gaussian.png" % (output_folder, dimension_labels[dimension], fixed_coordinates[ other_dims[0] ], fixed_coordinates[ other_dims[1] ])
    plt.savefig(plot_name)
#
#     # print data to screen
#     print ("p[0], a: ", v[0])
#     print ("p[1], mu: ", v[1])
#     print ("p[2], sigma: ", v[2])
#     print ("p[3], y-offset: ", v[3])
#     print ' FWHM (mm): {}'.format(FWHM)
#
    return points_gauss, values_gauss, FWHM


def plot_maximum_planes(obj3D):
    max_point = obj3D.get_maximum()
    plot = obj3D.get_plot(0, max_point[0])
    plot = obj3D.get_plot(1, max_point[1])
    plot = obj3D.get_plot(2, max_point[2])


def plot_maximum_planes_ignore_edge(obj3D, max_point, npx):
    #max_point = obj3D.get_maximum()
    plot = obj3D.get_plot_ignore_edge(0, max_point[0], npx)
    plot = obj3D.get_plot_ignore_edge(1, max_point[1], npx)
    plot = obj3D.get_plot_ignore_edge(2, max_point[2], npx)


def plot_max_profiles(my3DObj, output_folder_path):
    max_point = my3DObj.get_maximum()
    plot_profiles(my3DObj, output_folder_path, 0, [max_point[1], max_point[2]])
    plot_profiles(my3DObj, output_folder_path, 1, [max_point[0], max_point[2]])
    plot_profiles(my3DObj, output_folder_path, 2, [max_point[0], max_point[1]])


def plot_max_profiles_ignore_edge(my3DObj, output_folder_path, max_point):
    #max_point = my3DObj.get_maximum()
    plot_profiles(my3DObj, output_folder_path, 0, [max_point[1], max_point[2]])
    plot_profiles(my3DObj, output_folder_path, 1, [max_point[0], max_point[2]])
    plot_profiles(my3DObj, output_folder_path, 2, [max_point[0], max_point[1]])



def get_maximum_ignore_edges(my3DObj, npx = 1):

    # use masked array to ignore edges (indices of 0 or Max - 1)
    #  source: https://stackoverflow.com/questions/14611250/mask-specific-columns-of-a-numpy-array

    data = my3DObj.values
    # getting number of bins in each direction
    divis = numpy.array(my3DObj.num_bins)
    # creating mask index
    mask_index_min = numpy.arange(0, npx)
    mask_index_max = numpy.arange(numpy.max(divis) - npx, numpy.max(divis))
    mask_index = numpy.array( [mask_index_min, mask_index_max] )
    # create 3d array of zeros with same dimensions as data
    mask = numpy.zeros_like(data)
    # fill the edges with ones
    mask[mask_index, :, :] = 1  # does first and last planes
    mask[:, mask_index, :] = 1  # does rows
    mask[:, :, mask_index] = 1  # does columns
    # create masked data array
    mask_data = numpy.ma.masked_array(data, mask)
    # find array position of maximum value in masked data array
    i, j, k = numpy.unravel_index(mask_data.argmax(), mask_data.shape)
    # convert array location into position
    point = (my3DObj.bin_centers[0][i], my3DObj.bin_centers[1][j], my3DObj.bin_centers[2][k])

    return point


if __name__ == "__main__":

    def print_usage():
        print (USAGE % (sys.argv[0], sys.argv[0]))
        sys.exit(100)

    if len(sys.argv) == 3:

        data_file_path = sys.argv[1]
        output_folder_path = sys.argv[2]
    else:
        print_usage()

    #pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)


    my3DObj = Intensity3D.Intensity3D(data_file_path, output_folder_path)
    # my3DObj.get_plot_range(2, 0.0, -10, 10, "altered_range", True)
    #pool.apply_async(plot_y0_profiles, args=(my3DObj, output_folder_path))
    plot_y0_profiles(my3DObj, output_folder_path)

    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [0.0, 0.0])).get()
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [0.0, 0.0])).get()
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [0.0, 0.0])).get()
    plot_profiles(my3DObj, output_folder_path, 2, [0.0, 0.0])
    plot_profiles(my3DObj, output_folder_path, 1, [0.0, 0.0])
    plot_profiles(my3DObj, output_folder_path, 0, [0.0, 0.0])
    #pool.apply_async(plot_max_profiles, args=(my3DObj, output_folder_path)).get()
    #pool.apply_async(plot_maximum_planes, args=(my3DObj,)).get()
    #pool.apply_async(plot_generic_planes, args=(my3DObj,)).get()
    plot_generic_planes(my3DObj)


#==== MAX ====#

    # Output position of maximum value to screen
    max_pos = my3DObj.get_maximum()
    # convert numpy scalars to native floats for nicer printing and downstream formatting
    max_pos = (float(max_pos[0]), float(max_pos[1]), float(max_pos[2]))
    print ('\n  Position of maximum point in reconstruction (mm): {}'.format(max_pos))

    # Output position of maximum value (ignoring edges) to screen (only works for arrays with equal lengths)
    width_of_edge_in_pixels = 5  # default value is 1
    max_pos_ignore_edge = get_maximum_ignore_edges(my3DObj, width_of_edge_in_pixels)
    # convert numpy scalars to native floats for nicer printing
    max_pos_ignore_edge = (float(max_pos_ignore_edge[0]), float(max_pos_ignore_edge[1]), float(max_pos_ignore_edge[2]))
    print ('  Position of maximum point (ignoring edges with pixel width = {}) in reconstruction (mm): {}'.format(width_of_edge_in_pixels, max_pos_ignore_edge))

    # Should kick in for MLEM reconstructions
    if max_pos != max_pos_ignore_edge:
        max_pos = max_pos_ignore_edge
        # plot maximum profiles and heatmaps
        plot_max_profiles_ignore_edge(my3DObj, output_folder_path, max_pos)
        plot_maximum_planes_ignore_edge(my3DObj, max_pos, width_of_edge_in_pixels)
    else:
        # plot maximum profiles and heatmaps
        plot_max_profiles(my3DObj, output_folder_path)
        plot_maximum_planes(my3DObj)

    # save 1D profile data -> max_value
    print ('\n  Saving 1D profiles through maximum source position (mm): {}'.format(max_pos))
    # x-dir
    X, Values = my3DObj.get_profile(dimension=0, fixed_coordinates = [max_pos[1], max_pos[2]])
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_x_%.1f_%.1f_maximum.csv" % (output_folder_path, max_pos[1], max_pos[2])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))
    # y-dir
    X, Values = my3DObj.get_profile(dimension=1, fixed_coordinates = [max_pos[0], max_pos[2]])
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_y_%.1f_%.1f_maximum.csv" % (output_folder_path, max_pos[0], max_pos[2])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))
    # z-dir
    X, Values = my3DObj.get_profile(dimension=2, fixed_coordinates = [max_pos[0], max_pos[1]])
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_z_%.1f_%.1f_maximum.csv" % (output_folder_path, max_pos[0], max_pos[1])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))


#==== 180315 - RUN 1 (Co-60) & RUN 8 (Cs-137)  ====#
    #exp_pos = [15.0, -5.0, 11.0]  # run 1 & 8
    #exp_pos = [5.0, -5.0, 11.0]  # runs 2 & 7
    #exp_pos = [-5.0, -5.0, 11.0]  # runs 3 & 6
    #exp_pos = [-15.0, -5.0, 11.0]  # runs 4 & 5
    fwhm = []

#==== 180810 - RUN 1 (Co-60)  ====#
    #exp_pos = [73.0, -55.0, 11.0]  # run 1

#==== 190615 - ALL RUNS (Co-60, Cs-137, N1-22)  ====#
    #exp_pos = [-3.0, -6.0, 11.0]  # all runs

#==== 190805 - ALL RUNS (Na-22)  ====#
    #exp_pos = [-1.0, 0.0, 11.0]  # all runs

#==== 220615 - ALL RUNS (66 MeV proton beam)  ====#
    #exp_pos = [11.0, 0.0, 0.0]  # all runs, I think
    #exp_pos = [0.0, 0.0, 0.0]  # change to 0, 0, 0 since there is a -11 in x translation during coordinate_transformation

#==== 240502 - ALL RUNS (Co-60 test run)  ====#
    #exp_pos = [0.0, 0.0, 0.0]  # all runs
    exp_pos = EXPECTED_SOURCE_LOCATION  # all runs

    # plots maximum 1D profiles with expected points
    plot_profile_expected(my3DObj, output_folder_path, 0, max_pos, exp_pos)
    plot_profile_expected(my3DObj, output_folder_path, 1, max_pos, exp_pos)
    plot_profile_expected(my3DObj, output_folder_path, 2, max_pos, exp_pos)

    # plots maximum 1D profiles with gaussian fit and FWHM
    print ('\n  Saving 1D profiles (gaussian fits) through maximum source position (mm): {}'.format(max_pos))
    # x-dir
    X, Values, F = plot_profile_gaussian(my3DObj, output_folder_path, 0, max_pos, exp_pos)
    gpx = X[ numpy.where(Values == numpy.max(Values))[0] ][0]
    fwhm.append(F)
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_x_%.1f_%.1f_gaussian.csv" % (output_folder_path, max_pos[1], max_pos[2])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))
    X, Values, F = plot_profile_gaussian(my3DObj, output_folder_path, 1, max_pos, exp_pos)
    gpy = X[ numpy.where(Values == numpy.max(Values))[0] ][0]
    fwhm.append(F)
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_y_%.1f_%.1f_gaussian.csv" % (output_folder_path, max_pos[0], max_pos[2])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))
    X, Values, F = plot_profile_gaussian(my3DObj, output_folder_path, 2, max_pos, exp_pos)
    gpz = X[ numpy.where(Values == numpy.max(Values))[0] ][0]
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_z_%.1f_%.1f_gaussian.csv" % (output_folder_path, max_pos[0], max_pos[1])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))
    # saving max gaussian points and fwhm
    gau_pos = [gpx, gpy, gpz]
    fwhm.append(F)

    print ('  Maximum point in reconstruction from gaussian fits (mm): ({:.2f}, {:.2f}, {:.2f})'.format(gpx, gpy, gpz))

    # plot 1d profiles through expected source location
    print ('\n  Creating 1D profiles through expected source position (mm): {}'.format(exp_pos))
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [exp_pos[1], exp_pos[2]])).get() # YZ
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [exp_pos[0], exp_pos[2]])).get() # XZ
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [exp_pos[0], exp_pos[1]])).get() # XY
    plot_profiles(my3DObj, output_folder_path, 0, [exp_pos[1], exp_pos[2]])
    plot_profiles(my3DObj, output_folder_path, 1, [exp_pos[0], exp_pos[2]])
    plot_profiles(my3DObj, output_folder_path, 2, [exp_pos[0], exp_pos[1]])

    # save 1D profile data
    print ('  Saving 1D profiles through expected source position (mm): {}'.format(exp_pos))
    # x-dir
    X, Values = my3DObj.get_profile(dimension=0, fixed_coordinates = [exp_pos[1], exp_pos[2]])
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_x_%.1f_%.1f_actual.csv" % (output_folder_path, exp_pos[1], exp_pos[2])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))
    # y-dir
    X, Values = my3DObj.get_profile(dimension=1, fixed_coordinates = [exp_pos[0], exp_pos[2]])
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_y_%.1f_%.1f_actual.csv" % (output_folder_path, exp_pos[0], exp_pos[2])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))
    # z-dir
    X, Values = my3DObj.get_profile(dimension=2, fixed_coordinates = [exp_pos[0], exp_pos[1]])
    output_array = numpy.stack((X, Values), axis=-1)
    output_file = "%s/profile_z_%.1f_%.1f_actual.csv" % (output_folder_path, exp_pos[0], exp_pos[1])
    numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
    print ('   -> output_file: {}'.format(output_file))

    # plot 2d plane slices through expected source location
    print ('  Creating 2D profiles through expected source position (mm): {}'.format(exp_pos))
    my3DObj.get_plot(plane_index=0, depth=exp_pos[0], title="Known X Depth", is_normalized=False)
    my3DObj.get_plot(plane_index=1, depth=exp_pos[1], title="Known Y Depth", is_normalized=False)
    my3DObj.get_plot(plane_index=2, depth=exp_pos[2], title="Known Z Depth", is_normalized=False)

    # calculate the difference between the maximum and expected positions
    difference = numpy.subtract( max_pos, exp_pos )
    magn = numpy.linalg.norm(difference)

    # calculate the difference between the gaussian and expected positions
    difference_g = numpy.subtract( gau_pos, exp_pos )
    magn_g = numpy.linalg.norm(difference_g)

    # calculate combined FWHM
    avg_fwhm = numpy.average(fwhm)

    # calculate resolution of reconstruction
    x_start, x_end = min(my3DObj.edges[0]), max(my3DObj.edges[0])
    y_start, y_end = min(my3DObj.edges[1]), max(my3DObj.edges[1])
    z_start, z_end = min(my3DObj.edges[2]), max(my3DObj.edges[2])
    divis = my3DObj.num_bins
    limits = numpy.array ( ( [x_start, x_end] , [y_start, y_end] , [z_start, z_end] ) , dtype = float )
    deltas = [ (limit[1] - limit[0])/div for limit, div in zip(limits, divis) ]

    # print distance calculations to screen
    print( '\nActual Location of source:     {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm '.format( exp_pos[0] , exp_pos[1] , exp_pos[2] ) )
    print( 'Estimated Location of source:  {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm '.format( max_pos[0] , max_pos[1] , max_pos[2] ) )
    print( 'Gaussian Location of source:   {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm '.format( gau_pos[0] , gau_pos[1] , gau_pos[2] ) )
    print( '                                    _______________________________' )
    print( 'Difference (Estim - Actual):   {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm '.format( difference[0] , difference[1] , difference[2] ) )
    print( 'Magnitude  (Estim - Actual):   {:06.2f} mm '.format( magn ) )
    print( '                                    _______________________________' )
    print( 'Difference (Gauss - Actual):   {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm '.format( difference_g[0] , difference_g[1] , difference_g[2] ) )
    print( 'Magnitude  (Gauss - Actual):   {:06.2f} mm '.format( magn_g ) )
    print( '                                    _______________________________' )
    print( 'Full Width, Half Maximums:     {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm '.format( fwhm[0] , fwhm[1] , fwhm[2] ) )
    print( 'Average FWHM:                  {:06.2f} mm '.format( avg_fwhm ) )
    print( '\nResolution (voxel-size): {:.2f} x {:.2f} x {:.2f} mm'.format( deltas[0] , deltas[1] , deltas[2] ) )

    # save distance calculations to text file
    txt = open(output_folder_path + '/Distance_Calculations_FWHM.txt','w')
    txt.write( 'Actual Location of source:     {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm \n'.format( exp_pos[0] , exp_pos[1] , exp_pos[2] ) )
    txt.write( 'Estimated Location of source:  {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm \n'.format( max_pos[0] , max_pos[1] , max_pos[2] ) )
    txt.write( 'Gaussian Location of source:   {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm \n'.format( gau_pos[0] , gau_pos[1] , gau_pos[2] ) )
    txt.write( '                                    _______________________________\n' )
    txt.write( 'Difference (Estim - Actual):   {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm \n'.format( difference[0] , difference[1] , difference[2] ) )
    txt.write( 'Magnitude  (Estim - Actual):   {:06.2f} mm \n'.format( magn ) )
    txt.write( '                                    _______________________________\n' )
    txt.write( 'Difference (Gauss - Actual):   {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm \n'.format( difference_g[0] , difference_g[1] , difference_g[2] ) )
    txt.write( 'Magnitude  (Gauss - Actual):   {:06.2f} mm \n'.format( magn_g ) )
    txt.write( '                                    _______________________________\n' )
    txt.write( 'Full Width, Half Maximums:     {:06.2f} mm   ,   {:06.2f} mm   ,   {:06.2f} mm \n'.format( fwhm[0] , fwhm[1] , fwhm[2] ) )
    txt.write( 'Average FWHM:                  {:06.2f} mm \n'.format( avg_fwhm ) )
    txt.write( '\nResolution (voxel-size): {:.2f} x {:.2f} x {:.2f} mm\n'.format( deltas[0] , deltas[1] , deltas[2] ) )
    txt.close()


#==== 180425 - RUN 1 (Na-22) ====#
#
#     # plot 1d profiles through source [-1.5, 0.0, +12 mm]
#     exp_pos = [-3.5, 0.0, 12.0]
#     print "  Creating 1D profiles through expected source position for run 1 (mm):", exp_pos
#     pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [exp_pos[1], exp_pos[2]])).get() # YZ
#     pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [exp_pos[0], exp_pos[2]])).get() # XZ
#     pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [exp_pos[0], exp_pos[1]])).get() # XY
#
#     # save 1D profile data
#     print "  Saving 1D profiles through expected source position for run 1 (mm):", exp_pos
#     # x-dir
#     X, Values = my3DObj.get_profile(dimension=0, fixed_coordinates = [exp_pos[1], exp_pos[2]])
#     output_array = numpy.stack((X, Values), axis=-1)
#     output_file = "%s/profile_x_%.1f_%.1f.csv" % (output_folder_path, exp_pos[1], exp_pos[2])
#     numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
#     print "   -> output_file:", output_file
#     # y-dir
#     X, Values = my3DObj.get_profile(dimension=1, fixed_coordinates = [exp_pos[0], exp_pos[2]])
#     output_array = numpy.stack((X, Values), axis=-1)
#     output_file = "%s/profile_y_%.1f_%.1f.csv" % (output_folder_path, exp_pos[0], exp_pos[2])
#     numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
#     print "   -> output_file:", output_file
#     # z-dir
#     X, Values = my3DObj.get_profile(dimension=2, fixed_coordinates = [exp_pos[0], exp_pos[1]])
#     output_array = numpy.stack((X, Values), axis=-1)
#     output_file = "%s/profile_z_%.1f_%.1f.csv" % (output_folder_path, exp_pos[0], exp_pos[1])
#     numpy.savetxt(output_file, output_array, fmt='%1.4f', delimiter=",")
#     print "   -> output_file:", output_file
#
#     # plot 2d plane slices through source [-15, -5, +11 mm]
#     print "  Creating 2D profiles through expected source position for run 1 (mm):", exp_pos
#     my3DObj.get_plot(plane_index=0, depth=exp_pos[0], title="Known X Depth", is_normalized=False)
#     my3DObj.get_plot(plane_index=1, depth=exp_pos[1], title="Known Y Depth", is_normalized=False)
#     my3DObj.get_plot(plane_index=2, depth=exp_pos[2], title="Known Z Depth", is_normalized=False)
#
#
#    # plot 1d profiles through source [-15, -5, +11 mm] -> Run 4 & 5
#    print("  Creating 1D profiles through (-15, -5, +11 mm) -> Co60 source location / run 4 & 5")
#    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [-5.0, 11.0])).get()  # YZ
#    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [-15.0, 11.0])).get() # XZ
#    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [-15.0, -5.0])).get() # XY

#    # plot 1d profiles through source [+10, +10, -5 mm] -> Run 5
#    print("  Creating 1D profiles through (+10, +10, -5 mm) -> Cs137 source location / run 5")
#    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [10.0, -5.0])).get()  # YZ
#    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [10.0, -5.0])).get()  # XZ
#    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [10.0, 10.0])).get()  # XY

#    # plot 2d plane slices through source [-21, +6, +11 mm]
#    print("  Creating 2D profiles through (-21, +6, +11 mm)")
#    my3DObj.get_plot(plane_index=0, depth=-21, title="Known X Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=1, depth=6, title="Known Y Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=2, depth=11, title="Known Z Depth", is_normalized=False)

#    # plot 2d plane slices through source [0, +5, +11 mm] -> Run 1
#    print("  Creating 2D profiles through (0, +5, +11 mm) -> Co60 source location / run 1")
#    my3DObj.get_plot(plane_index=0, depth=0, title="Known X Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=1, depth=5, title="Known Y Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=2, depth=11, title="Known Z Depth", is_normalized=False)

#    # plot 2d plane slices through source [15, +5, +11 mm] -> Run 2
#    print("  Creating 2D profiles through (15, +5, +11 mm) -> Co60 source location / run 2")
#    my3DObj.get_plot(plane_index=0, depth=15, title="Known X Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=1, depth=5, title="Known Y Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=2, depth=11, title="Known Z Depth", is_normalized=False)



#    # plot 2d plane slices through source [-15, -5, +11 mm] -> Run 4 & 5
#    print("  Creating 2D profiles through (-15, -5, +11 mm) -> Co60 source location / run 4 & 5")
#    my3DObj.get_plot(plane_index=0, depth=-15, title="Known X Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=1, depth=-5, title="Known Y Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=2, depth=11, title="Known Z Depth", is_normalized=False)

#    # plot 2d plane slices through source [+10, +10, -5 mm] -> Run 5
#    print("  Creating 2D profiles through (+10, +10, -5 mm) -> Cs137 source location / run 5")
#    my3DObj.get_plot(plane_index=0, depth=10, title="Known X Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=1, depth=10, title="Known Y Depth", is_normalized=False)
#    my3DObj.get_plot(plane_index=2, depth=-5, title="Known Z Depth", is_normalized=False)

    # pull 1d profile for X at Y = 6, Z = 11
    #X, Values = my3DObj.get_profile(dimension=0, fixed_coordinates = [6.0, 11.0])
    #Y, Values = my3DObj.get_profile(dimension=1, fixed_coordinates = [-21.0, 11.0])
    #Z, Values = my3DObj.get_profile(dimension=2, fixed_coordinates = [-21.0, 6.0])
