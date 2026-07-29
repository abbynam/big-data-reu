import numpy as np
import sys
import shutil
import os

path = sys.argv[1]
name = sys.argv[2]

original_train = np.load(path+name+"/train/X.npy")
original_test = np.load(path+name+"/test/X.npy")

reshaped_train = original_train.reshape(-1, 3, 5)
reshaped_test = original_test.reshape(-1,3,5)

os.makedirs(path+name+"_3x5/train/", exist_ok=True)
os.makedirs(path+name+"_3x5/test/", exist_ok=True)

np.save(path+name+"_3x5/train/X.npy", reshaped_train)
np.save(path+name+"_3x5/test/X.npy", reshaped_test)

shutil.copy(path+name+"/train/y.npy", path+name+"_3x5/train/y.npy")
shutil.copy(path+name+"/test/y.npy", path+name+"_3x5/test/y.npy")
