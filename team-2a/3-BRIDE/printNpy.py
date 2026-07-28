#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd


def print_npy(input_file, line_count):
    arr = np.load(input_file)

    print(f"Size: {arr.shape}")
    print(arr[:line_count])  # first 5 rows/elements


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage:")
        print("python printNpy.py input.npy [line_count]")
        sys.exit(1)

    input_txt = sys.argv[1]
    # default line count
    line_count = 15

    if len(sys.argv) == 3:
            line_count = int(sys.argv[2])

    print_npy(input_txt, line_count)

