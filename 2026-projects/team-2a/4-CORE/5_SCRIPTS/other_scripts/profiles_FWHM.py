#! /usr/bin/python

import sys, os
import numpy
import scipy
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt

import Intensity3D

#---------- CONSTANTS -----------------

USAGE = '''
USAGE: %s [Path to 3D dose file] [output path]

The following example generates plots using the image intensity information from the file output.dat.

   python %s output.dat outputfolder
'''

#---------- PLOTTING FUNCTIONS -----------------
def get_FWHM(x, y):

    x_max_index = numpy.argmax(y)
    half_max = 0.5 * max(y)

    x_right = 0
    x_left = 0

    x_right = x_max_index
    for i in range(int(x_max_index), 10000):
        if y[i] < half_max:
            x_right = x[i]
            break
        else:
            print i, y[i], half_max

    for i in range(0, 10000):
        if y[x_max_index - i] < half_max:
            x_left = x[x_max_index - i]
            break
        else:
            print i, y[x_max_index - i], half_max

    fwhm = x_right - x_left

    return fwhm


def plot_profiles(obj3Ds, output_folder, dimension=1):

    plt.clf()
    ax = plt.subplot(111)

    labels = [r'$\Delta$E', r'$\Delta$$\Theta$', 'window', 'basic']
    for i,obj in enumerate(obj3Ds):
      
        #fixed coordinates are the values for the coorinated not being profiled. 
        #  If dimension is 2 then (0,0) produces a profile along the z axis.
        x, intensities = obj.get_profile(dimension=dimension, fixed_coordinates=[0.0,0.0])
        fwhm = get_FWHM(x, intensities)
        intensities = intensities/max(intensities)
        plt.plot(x, intensities, label="%s: (%.1f mm)" % (labels[i], fwhm), lw=1)


    plt.title("Co-60 Source Profile", fontsize=11)
    legend = ax.legend(loc='upper left', title=r'Filter: (FWHM)', fontsize=8)
    ax.grid(True,linestyle='-',which='both',color='0.750')

    plt.setp(ax.get_xticklabels(), fontsize=9)
    plt.setp(ax.get_yticklabels(), fontsize=9)
    plt.setp(legend.get_title(), fontsize=9)
    plt.xlabel("position [mm]", fontsize=10)

    plt.ylabel("Fraction of Maximum", fontsize=10)
    plt.gcf().set_size_inches(9,6)
    plt.gcf().subplots_adjust(bottom=0.15)

    plot_name = "%s/profile_FWHM.png" % (output_folder)
    plt.savefig(plot_name)

    return True

if __name__ == "__main__":
    OUTPUT_FOLDER = '/home/dsmackin/public_html/dtheta'
    # def print_usage():
    #     print USAGE % (sys.argv[0], sys.argv[0])
    #     sys.exit(100)
    #
    # if len(sys.argv) == 3:
    #
    #     data_file_path = sys.argv[1]
    #     output_folder_path = sys.argv[2]
    # else:
    #     print_usage()

    # rawObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/out/Co60c_raw/output.dat")
    # windowObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/out/Co60c_window/output.dat")
    # thetaObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/out/Co60c_dthetam/output.dat")
    # energyObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/out/Co60c_dE/output.dat")
    #
    # rawObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-MAXIMUM_NUMBER_CONES-3000__TYPE-m__MIN_GAMMA_ENERGY-0.0__MAX_GAMMA_ENERGY-10.0/output.dat")
    # windowObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-MAXIMUM_NUMBER_CONES-3000__TYPE-m__MIN_GAMMA_ENERGY-1.0__MAX_GAMMA_ENERGY-1.5/output.dat")
    # thetaObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-MAXIMUM_NUMBER_CONES-3000__TYPE-m_dTheta__MIN_GAMMA_ENERGY-0.0__MAX_GAMMA_ENERGY-10.0/output.dat")
    # energyObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-MAXIMUM_NUMBER_CONES-3000__TYPE-m_dE__MIN_GAMMA_ENERGY-0.0__MAX_GAMMA_ENERGY-10.0/output.dat")

    rawObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-m__CONE_LENGTH_CORRECTION-1.0__MAXIMUM_NUMBER_CONES-10000/output.dat")
    windowObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-m__CONE_LENGTH_CORRECTION-0.90__MAXIMUM_NUMBER_CONES-10000/output.dat")
    thetaObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-m_dTheta__CONE_LENGTH_CORRECTION-1.0__MAXIMUM_NUMBER_CONES-10000/output.dat")
    energyObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/dtheta/dtheta-TYPE-m_dE__CONE_LENGTH_CORRECTION-0.90__MAXIMUM_NUMBER_CONES-10000/output.dat")


    plot_profiles([energyObj, thetaObj, windowObj, rawObj], OUTPUT_FOLDER, 2)

