#! /usr/bin/python

import sys, os
import numpy
import scipy
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing

import Intensity3D

USAGE = '''
USAGE: %s [Path to 3D dose file] [output path]

The following example generates plots using the image intensity information from the file output.dat.

   python %s output.dat outputfolder
'''


def plot_generic_planes(obj3D, num_planes_per_dimension = 5):
    ranges = obj3D.get_ranges()
    assert(len(ranges) == 3)
    for i,r in enumerate(ranges):
        step = (r[1] - r[0])/(num_planes_per_dimension + 1)
        depths = numpy.linspace(r[0] + 0.5*step, r[1] - 0.5*step, num_planes_per_dimension)
        for depth in depths:
            plot = obj3D.get_plot(i, depth)


def plot_y0_profiles(obj3D, output_folder):
    plane = obj3D.getIntensityPlane(1, 0.0)
    k,i = numpy.unravel_index(plane.argmax(), plane.shape)

    x, tmp, z = obj3D.get_coordinates(i, 0, k)

    plot_profiles(obj3D, output_folder, dimension=0, fixed_coordinates = [0.0, z])
    plot_profiles(obj3D, output_folder, dimension=2, fixed_coordinates = [x, 0.0])



def plot_profiles(obj3D, output_folder, dimension=2, fixed_coordinates = [0.0, 0.0], normalize=False):
    try:
        points, values = obj3D.get_profile(dimension, fixed_coordinates)
    except ValueError, ve:
        print str(ve)
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

    plt.setp(ax.get_xticklabels(), rotation='30', fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    plt.setp(legend.get_title(), fontsize=6)
    plt.xlabel("%s [mm]" % dimension_labels[dimension], fontsize=10)

    plt.ylabel("Gamma Emission", fontsize=10)
    plt.gcf().set_size_inches(5,3.3)
    plt.gcf().subplots_adjust(bottom=0.15)

    plot_name = "%s/profile_%s_%.1f_%.1f.png" % (output_folder, dimension_labels[dimension], fixed_coordinates[0], fixed_coordinates[1])
    plt.savefig(plot_name)

    return values


def plot_maximum_planes(obj3D):
    max_point = obj3D.get_maximum()
    print("MAXPOINT: %.2f, %.2f, %.2f" % (max_point[0], max_point[1], max_point[2]))
    plot = obj3D.get_plot(0, max_point[0])
    plot = obj3D.get_plot(1, max_point[1])
    plot = obj3D.get_plot(2, max_point[2])


def plot_max_profiles(my3DObj, output_folder_path):
    max_point = my3DObj.get_maximum()
    plot_profiles(my3DObj, output_folder_path, 0, [max_point[1], max_point[2]])
    plot_profiles(my3DObj, output_folder_path, 1, [max_point[0], max_point[2]])
    plot_profiles(my3DObj, output_folder_path, 2, [max_point[0], max_point[1]])



if __name__ == "__main__":

    def print_usage():
        print USAGE % (sys.argv[0], sys.argv[0])
        sys.exit(100)

    if len(sys.argv) == 3:

        data_file_path = sys.argv[1]
        output_folder_path = sys.argv[2]
    else:
        print_usage()

    pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)


    my3DObj = Intensity3D.Intensity3D(data_file_path, output_folder_path)

    my3DObj.get_plot(1, 0.0)
    my3DObj.get_plot(2, 0.0)
    my3DObj.get_plot(2, -38)
    my3DObj.get_plot(0, 50)
    my3DObj.get_plot(1, 25)
    for x in (10.0, 30.0, 50.0, 75.0, 80.0, 90.0):
        print(x)
        my3DObj.get_plot(0,x)

    # my3DObj.get_plot_range(2, 0.0, -10, 10, "altered_range", True)
    pool.apply_async(plot_y0_profiles, args=(my3DObj, output_folder_path))
    #pool.apply_async(my3DObj.get_plot, args=(1,0.0)).get()
    #pool.apply_async(my3DObj.get_plot, args=(2,0.0)).get()
    #pool.apply_async(my3DObj.get_plot, args=(1,3.3)).get()

    pool.apply_async(plot_max_profiles, args=(my3DObj, output_folder_path)).get()
    pool.apply_async(plot_maximum_planes, args=(my3DObj,)).get()
    #pool.apply_async(plot_generic_planes, args=(my3DObj,)).get()
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [6.0, 11.0])).get()
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [-21.0, 11.0])).get()
    #pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [-21.0, 6.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [0.0, 0.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [0.0, 0.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [0.0, 0.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [0.0, 0.0])).get()
    

        
    my3DObj.get_plot(1, 0)
    my3DObj.get_plot(2, 0)
    #my3DObj.get_plot(0, 300.0)
    #my3DObj.get_plot(0, 400.0)
    #my3DObj.get_plot(0, 350.0)
