'''
08/31/2018

The following code converts the reconData of the *.root file obtained from the POLARIS3.
into a *.csv file. It requires "uproot" preinstalled on the computer and may also require lz4.
If not already installed, these packages can be installed using:

$pip install uproot --user
$pip install lz4 --user

Author :: Jerimy Polf / S W Peterson

Additionally, to read the strings from Ttree branches/leaves, uproot should be upgraded to
version 3.0.0b2 (currently a beta version). Upgrading can be done using:

$pip install -U "uproot>=3.0.0b2"

Modified: V R SHARMA
Modified: H LEWIS 2025-06-20
'''

import uproot
import awkward as ak
import numpy as np
import os
import argparse
import sys

# command line args
parser = argparse.ArgumentParser()
parser.add_argument("folder", type=str)
args = parser.parse_args()
folder_path = args.folder

if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a valid directory.")
        sys.exit(1)

print(f"Received folder: {folder_path}")

MU = 100 #---->  This is a monitor unit case [Eliminated for current project]

# change me
NAME_PREFIX = "compiled"

# load the data
DATA_FILE = folder_path + "/root/" + NAME_PREFIX + ".root"
TREE_NAME = "scatterData"

# make csv folder if needed
if not os.path.isdir(folder_path + "/csv/"):
        os.makedirs(folder_path + "/csv/")

# output csv
OUT_FILE = folder_path + '/csv/' + NAME_PREFIX

MyFile = uproot.open(DATA_FILE)
MyTree = MyFile[TREE_NAME]

#Form an array for each leaf
gevent = MyTree["event"].array()
gtrack = MyTree["track"].array() #[7 4 4 ... 6 4 7]
gstep = MyTree["step"].array()
gdetPos = MyTree["detPos"].array()
gdetSize = MyTree["detSize"].array()
gpos_x = MyTree["pos_x"].array()
gpos_y = MyTree["pos_y"].array()
gpos_z = MyTree["pos_z"].array()
ggammaEnergy = MyTree["gammaEnergy"].array()
genergyDeposited = MyTree["energyDeposited"].array()
gscatAng = MyTree["scatAng"].array()
gorigin_x = MyTree["origin_x"].array()
gorigin_y = MyTree["origin_y"].array()
gorigin_z = MyTree["origin_z"].array()
gorigin_energy = MyTree["origin_energy"].array()

#Array of each branch (??)
gdetector = MyTree["detector"].array()
gprocess = MyTree["process"].array()
gorigin_process = MyTree["origin_process"].array()
gorigin_volume = MyTree["origin_volume"].array()

N = len(gevent)
print('Number of events found: %i: ' % (N))

