#! /usr/bin/python

import sys, os
import numpy
import scipy
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt

import ROOT

sys.path.insert(0, './scripts')
import Intensity3D


def plot_profiles(obj3Ds, output_folder, dimension=2):

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
    plt.gcf().set_size_inches(5,3.3)
    plt.gcf().subplots_adjust(bottom=0.15)

    plot_name = "%s/profile_FWHM.png" % (output_folder)
    plt.savefig(plot_name)

    return True

def read_dose_file(filepath):
    with open(filepath) as f:
        content = f.readlines()

    X = numpy.zeros(len(content) - 1)
    Y = X + 0
    for i, line in enumerate(content[1:]):
        X[i], Y[i] = (float(v) for v in line.strip().split(","))

    return X, Y

def read_PG_origins_from_root(rootfile):

    f = ROOT.TFile(rootfile)
    myTree = f.Get("tripleGammas")
    z_vals = numpy.zeros(myTree.GetEntries())

    for i in range(myTree.GetEntries()):
        myTree.GetEntry(i)
        z_vals[i] = myTree.origin_z

    return z_vals

def kde(x, data, bandwidth=10):
    invbandwidth = 1.0/bandwidth

    y = invbandwidth * len(data) * numpy.sum( numpy.exp(-(data - x)*(data - x)*invbandwidth*invbandwidth))
    # for v in data:
    #     y +=
    # Y = numpy.array([numpy.sum( numpy.exp(-(data - x)*invbandwidth)) for x in X])
    # Y /= max(Y)

    return y


def kde_curve(X, data, bandwidth=10):
    invbandwidth = 1.0/bandwidth
    Y = X * 0.0
    for i,x in enumerate(X):
        Y[i] = kde(x, data, bandwidth)
    # Y = numpy.array([numpy.sum( numpy.exp(-(data - x)*invbandwidth)) for x in X])
    Y /= max(Y)

    return Y


def plot_dose(depth, dose, origins, obj3Ds, outputfolder, dimension=2):
    plt.clf()
    ax = plt.subplot(111)
    depth *= 10.0
    dose /= numpy.max(dose)

    temp = numpy.array([163.0, 161.0, 159.0, 157.0, 155.0, 154.0, 153.0])
    newdepths = numpy.concatenate((temp, depth))

    pg_density = kde_curve(newdepths, origins, bandwidth=8)
    # pg_density = depth*0.00 +1

    # labels = [r'$\Delta$E', r'$\Delta$$\Theta$', 'window', 'basic']
    labels = [r'Single CC', r'Orthogonal', r'Parallel']
    for i,obj in enumerate(obj3Ds):

        #fixed coordinates are the values for the coorinated not being profiled.
        #  If dimension is 2 then (0,0) produces a profile along the z axis.
        x, intensities = obj.get_profile(dimension=dimension, fixed_coordinates=[0.0,0.0])
        intensities = intensities/max(intensities)
        plt.plot(x, intensities, label="%s" % (labels[i]), lw=2)

    plt.plot(depth, dose, label="Dose", lw=2, color='black')
    plt.plot(newdepths, pg_density, label="PG Origin", lw=1, color='black', linestyle = 'dashed')

    plt.title("200 MeV Proton Beam", fontsize=18)
    legend = ax.legend(loc='upper left', title=r'Configurations', fontsize=10)
    ax.grid(True,linestyle='-',which='both',color='0.750')

    plt.setp(ax.get_xticklabels(), fontsize=14)
    plt.setp(ax.get_yticklabels(), fontsize=14)
    plt.setp(legend.get_title(), fontsize=10)
    plt.xlabel("Depth [mm]", fontsize=16)

    plt.ylabel("Fraction of Maximum", fontsize=16)
    plt.gcf().set_size_inches(10,6)
    plt.gcf().subplots_adjust(bottom=0.15)

    plot_name = "%s/profile_dose.png" % (outputfolder)
    plt.savefig(plot_name)

    return True

#------------------------   MAIN -----------------------------------------
USAGE = '''
USAGE: %s [Path to 3D dose file] [output path]

The following example generates plots using the image intensity information from the file output.dat.

   python %s output.dat outputfolder
'''


def main():

    output_folder = '/home/dsmackin/public_html/3D'
    dosefile = '/y_drive/CCData/MC_spot/DoseFiles/200MeV-water_1DdepthDose.txt'
    rootfile = '/y_drive/CCData/MC_spot/200_270.root'

    pg_origin_z = read_PG_origins_from_root(rootfile)
    depth, dose = read_dose_file(dosefile)

    # singleObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/batch/dose-TYPE-200_270_direct__IMAGE_ALGORITHM-OCTANE__MAXIMUM_NUMBER_CONES-20000/output.dat")
    # orthoObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/batch/dose-TYPE-200_ortho__IMAGE_ALGORITHM-OCTANE__MAXIMUM_NUMBER_CONES-20000/output.dat")
    # parallelObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/batch/dose-TYPE-200_parallel__IMAGE_ALGORITHM-OCTANE__MAXIMUM_NUMBER_CONES-20000/output.dat")

    singleObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/batch/dose-TYPE-200_270_direct__IMAGE_ALGORITHM-SOE__MAXIMUM_NUMBER_CONES-20000/output.dat")
    orthoObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/batch/dose-TYPE-200_ortho__IMAGE_ALGORITHM-SOE__MAXIMUM_NUMBER_CONES-20000/output.dat")
    parallelObj = Intensity3D.Intensity3D("/home/dsmackin/projects/CORE/batch/dose-TYPE-200_parallel__IMAGE_ALGORITHM-SOE__MAXIMUM_NUMBER_CONES-20000/output.dat")

    # plot_dose(depth, dose, pg_origin_z, [], output_folder)
    plot_dose(depth, dose, pg_origin_z, [singleObj, orthoObj, parallelObj], output_folder)

if __name__ == "__main__":
    import cProfile
    cProfile.run("main()")
    # main()