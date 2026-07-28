#!/usr/bin/env python3

'''
This file should take in 12-column data and reorient it two different ways,
from either 90 or 270 degree gantry angle to -y, under bench data
(the inverse operation should be the same)

Input: csv
[e1, x1, y1, z1, e2, x2, y2, z2, e3, x3, y3, z3]

Output 90:
x -> y
y -> -x
[e1, y1, -x1, z1, e2, y2, -x2, z2, e3, y3, -x3, z3]

Output 270:
x -> -y
y -> x
[e1, -y1, x1, z1, e2, -y2, x2, z2, e3, -y3, x3, z3]

'''

import sys
import numpy as np
import pandas as pd


def repair_data(data, labels):

    X = np.load(data)
    print(f"Data rows: {X.shape[0]}")

    #y = pd.read_csv(labels, header=None).to_numpy().flatten()
    y = np.load(labels).flatten()
    print(f"Label rows: {y.shape[0]}")

    if X.shape[0] != y.shape[0]:
        print("Data rows != Label rows\nExiting...")
        sys.exit(1)

    # define the output struct
    processed_rows = []

    # define sorting rules
    RULES = {
        0:  (0, 1, 2),
        1:  (0, 2, 1),
        2:  (1, 0, 2),
        3:  (1, 2, 0),
        4:  (2, 0, 1),
        5:  (2, 1, 0),
        6:  (0, 1, None),
        7:  (1, 0, None),
        8:  (0, 2, None),
        9:  (2, 0, None),
        10: (1, 2, None),
        11: (2, 1, None),
    }

    # define counter
    COUNT = {
        0:  0,
        1:  0,
        2:  0,
        3:  0,
        4:  0,
        5:  0,
        6:  0,
        7:  0,
        8:  0,
        9:  0,
        10: 0,
        11: 0,
        12: 0,
    }

    # define a blank interaction
    BLANK = ["", "", "", ""]

    # loop through all the rows
    for row, label in zip(X,y):

        # increment counter
        COUNT[label]+=1

        # if the interactions are three unrelated singles, just skip that row
        if label ==  12:
            continue

        # grab the three interactions
        scatters=[
            list(row[0:4]),
            list(row[5:9]),
            list(row[10:14])
        ]

        # define output row
        reordered=[]

        # select rule
        rule = RULES[label]

        # follow the indices from the rule
        for idx in rule:

            # add blank row
            if idx is None:
                reordered.extend(BLANK)

            # add interactions in the new order
            else:
                reordered.extend(scatters[idx])

        # add reordered row to output
        processed_rows.append(reordered)

    print(f"0 : {COUNT[0 ]}")
    print(f"1 : {COUNT[1 ]}")
    print(f"2 : {COUNT[2 ]}")
    print(f"3 : {COUNT[3 ]}")
    print(f"4 : {COUNT[4 ]}")
    print(f"5 : {COUNT[5 ]}")
    print(f"6 : {COUNT[6 ]}")
    print(f"7 : {COUNT[7 ]}")
    print(f"8 : {COUNT[8 ]}")
    print(f"9 : {COUNT[9 ]}")
    print(f"10: {COUNT[10]}")
    print(f"11: {COUNT[11]}")
    total_count = COUNT[0]+COUNT[1]+COUNT[2]+COUNT[3]+COUNT[4]+COUNT[5]+COUNT[6]+COUNT[7]+COUNT[8]+COUNT[9]+COUNT[10]+COUNT[11]+COUNT[12]
    print(f"12: {COUNT[12]} (removed) Remaining: {total_count-COUNT[12]}")
    print(f"Total: {total_count}")

    # output to csv
    df = pd.DataFrame(processed_rows)
    # df.to_csv(
    #     "repaired.csv",
    #     index=False,
    #     header=False
    # )
    return df

def reorient_data(csv, gantry_angle):
    '''
    oreint  for 90G we're changing +X to +Y and +Y to -X and then for 270G changing +X to -Y and +Y to +X
    reverse orient for -90g we're changing +X to -Y and +Y to +X and then for - 270 changing +X to +Y and +Y to -X
    '''
    print(f"csv input: {csv}")
    #print(f"csv input: {csv.split('/')}")
    #print(f"csv input: {csv.split('/')[-1]}")
    #print(f"csv input: {csv.split('/')[-1].split('.')[0]}")

    data = pd.read_csv(csv, header=None)

    print(f"Data rows: {data.shape[0]}")

    # Iteratng through all x's
    for i in range(1, 12, 4):
        if gantry_angle == '90':
            neg_x = -data.iloc[:,i].to_numpy()
            data.iloc[:,i] = data.iloc[:, i+1].to_numpy() # x takes on the values of y
            data.iloc[:,i+1] = neg_x # y takes on the values of negative x

        if gantry_angle == '270':
            neg_y = -data.iloc[:,i+1].to_numpy()
            data.iloc[:,i+1] = data.iloc[:,i].to_numpy() # y takes on the values of x
            data.iloc[:,i] = neg_y # x takes on the values of y

        if gantry_angle == '-90':
            neg_y = -data.iloc[:,i+1].to_numpy()
            data.iloc[:,i+1] = data.iloc[:,i]
            data.iloc[:,i] = neg_y

        if gantry_angle == '-270':
            neg_x = -data.iloc[:,i].to_numpy()
            data.iloc[:,i] = data.iloc[:,i+1]
            data.iloc[:,i+1] = neg_x

    out = f"{csv.split('/')[-1].split('.')[0]}_reoriented{gantry_angle}.csv"

    data.to_csv(
       out,
       index=False,
       header=False
    )

    print(f"csv output: {out}")


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage:")
        print(f"python {sys.argv[0]} data.csv gantry_angle") # input.csv integer:90/270
        sys.exit(1)

    csv = sys.argv[1]
    gantry_angle = sys.argv[2]
    reorient_data(csv, gantry_angle)
