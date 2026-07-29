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
USAGE: %s [Path to CORE output file] [output path]

This script converts the CORE output to the dimensions specified in the script constants (originally 512x512x1200)

   python %s output.dat outputfolder 
'''
X_RANGE = [-100, 412, 512]
Y_RANGE = [-200, 312, 512]
Z_RANGE = [50, 150, 100]


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


    # my3DObj.get_plot_range(2, 0.0, -10, 10, "altered_range", True)
    pool.apply_async(plot_y0_profiles, args=(my3DObj, output_folder_path))

    pool.apply_async(plot_max_profiles, args=(my3DObj, output_folder_path)).get()
    pool.apply_async(plot_maximum_planes, args=(my3DObj,)).get()
    pool.apply_async(plot_generic_planes, args=(my3DObj,)).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [6.0, 11.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [-21.0, 11.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [-21.0, 6.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 2, [0.0, 0.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 1, [0.0, 0.0])).get()
    pool.apply_async(plot_profiles, args=(my3DObj, output_folder_path, 0, [0.0, 0.0])).get()
    my3DObj.get_plot(0, -16.0)
    #my3DObj.get_plot(0, -18.0)
    my3DObj.get_plot(0, -21.0)
    #my3DObj.get_plot(0, -23.0)
    #my3DObj.get_plot(1, 5.0)
    my3DObj.get_plot(1, 6.0)
    my3DObj.get_plot(1, 7.0)
    #my3DObj.get_plot(1, 8.0)
    #my3DObj.get_plot(2, 8.0)
    my3DObj.get_plot(2, 11.0)
    #my3DObj.get_plot(2, 14.0)
    #my3DObj.get_plot(2, 16.0)
