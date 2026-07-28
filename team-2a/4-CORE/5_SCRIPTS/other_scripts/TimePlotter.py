#!/usr/bin/python
"""
Description:
This script creates timing plots for the PolarisJ2 detector. The intent is to show that we can
discern the cyclotron frequency in the data.
"""
__author__ = "Dennis Mackin <dsmackin@mdanderson.org>"
__date__ = "Aug. 26, 2016"
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
import pyximport; pyximport.install()
import utilities
import StringIO
import numpy

def get_records_from_text(text):
    df = pandas.read_csv(StringIO.StringIO(text), sep=',', names=['scatters', 'detector', 'energy', 'x', 'y', 'z', 'time'])

    return df


def convert_text_to_csv(text):
    # text = re.sub(r'[;,\t ]*[;,\t ]', ",", text) #re is 10x slower!
    return text.replace(" ", ",").replace("\t", ",").replace(",,,,", ",").replace(",,,", ",").replace(",,", ",").replace("\r", "")


def read_cc_data_file(file_path):

    with open(file_path, "r") as f:
        for rec in f:
            # rec = rec.strip()
            if rec[0] != '0' and rec[0] != "#":
                break
        text = rec + "".join(f)

    print "converting the file to CSV, could take several minutes . . ."
    text = convert_text_to_csv(text)
    print "converting to dataframe . . ."
    df = get_records_from_text(text)

    return df


def plot_timing_info(df, output_folder, bins=1000, length_in_seconds=1.0, offset_in_seconds=0):

    print "# Records: %s" % (len(df))
    detector = df.detector.as_matrix().astype(int)
    timestamps = df.time.as_matrix().astype(float)
    print "min(%s), max(%s)" % (numpy.min(timestamps), numpy.max(timestamps))
    timestamps -= numpy.min(timestamps)
    print "min(%s), max(%s)" % (numpy.min(timestamps), numpy.max(timestamps))
    #convert from 10's of nanoseconds to nanoseconds then to seconds
    timestamps *= 10.0 * 1E-9
    df = pandas.DataFrame({'detector':detector, 'timestamp': timestamps})
    df = df[df.timestamp >= offset_in_seconds]
    df = df[df.timestamp < offset_in_seconds + length_in_seconds]

    detector_numbers = numpy.unique(detector)

    for i in detector_numbers:

        #Select based on the detector. Verify that the
        # detector and timestamps are still lined up correctly
        assert(len(detector) == len(timestamps))
        t_d = df[df.detector == i].timestamp

        num_events = len(t_d)
        print "Detector %d: %d events" % (i, len(t_d))
        if num_events == 0:
            t_d.append(-1)
            t_d.append(-2)
        utilities.plot_time(t_d, bins, 'counts per %.1f ms' % (1000.0 * length_in_seconds / float(bins)), "Timing Detector %02d" % i, "timing_D%02d" % i, output_folder)

    print numpy.min(timestamps)



#------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------
def usage():
    print "USAGE: %s [event file] [output folder]" % (sys.argv[0])
    print "\nExample:\n %s CC_events.csv ~/public_html/scratch" % (sys.argv[0])
    sys.exit(-1)


def main():

    if len(sys.argv) != 3:
        usage()

    input_file = sys.argv[1]
    output_folder = sys.argv[2]

    print "Reading in CC events from %s . . ." % input_file
    df_cc = read_cc_data_file(input_file)

    print "making time plots . . ."
    plot_timing_info(df_cc, output_folder, 1000, length_in_seconds=100.0, offset_in_seconds=1.0)



if __name__ == "__main__":

    #cProfile.run("main()", sort=2)
    main()

