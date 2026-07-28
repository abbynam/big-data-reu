#!/usr/bin/env python3

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
    df.to_csv(
        "repaired.csv",
        index=False,
        header=False
    )

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage:")
        print(f"python {sys.argv[0]} data.npy labels.npy")
        sys.exit(1)

    data_npy = sys.argv[1]
    labels_txt = sys.argv[2]

    repair_data(data_npy, labels_txt)

