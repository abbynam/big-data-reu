#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd


def euclidean_distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))


def process_csv(csv_path, output_path="X.npy"):
    # Read CSV without headers
    df = pd.read_csv(csv_path)

    if df.shape[1] != 12:
        raise ValueError(
            f"Expected 12 columns, got {df.shape[1]}"
        )

    processed_rows = []

    for _, row in df.iterrows():

        (
            e1, x1, y1, z1,
            e2, x2, y2, z2,
            e3, x3, y3, z3
        ) = row.tolist()

        p1 = (x1, y1, z1)
        p2 = (x2, y2, z2)
        p3 = (x3, y3, z3)

        # Pairwise Euclidean distances
        d12 = euclidean_distance(p1, p2)
        d23 = euclidean_distance(p2, p3)
        d13 = euclidean_distance(p1, p3)

        # Match observed 5-column repeating structure
        processed = [
            e1, x1, y1, z1, d12,
            e2, x2, y2, z2, d23,
            e3, x3, y3, z3, d13
        ]

        processed_rows.append(processed)

    X = np.array(processed_rows, dtype=np.float32)
    Y = np.zeros(X.shape[0], dtype=np.int64)

    np.save(output_path, X)
    np.save("y.npy", Y)

    print(f"Saved {output_path}")
    print(f"Shape: {X.shape}")

    print(f"Saved y.npy")
    print(f"Shape: {Y.shape}")

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage:")
        print("python csvToNpy.py input.csv [output.npy]")
        sys.exit(1)

    input_csv = sys.argv[1]

    if len(sys.argv) >= 3:
        output_npy = sys.argv[2]
    else:
        output_npy = "X.npy"

    process_csv(input_csv, output_npy)
