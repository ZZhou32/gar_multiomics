#!/usr/bin/env python3
import sys
from collections import defaultdict

def main():
    matrix = defaultdict(dict)
    samples = set()

    for line_num, line in enumerate(sys.stdin, 1):
        line = line.strip()
        if not line:
            continue

        parts = line.split('\t')
        if len(parts) != 3:
            print(f"Skipping malformed line {line_num}: {line}", file=sys.stderr)
            continue

        try:
            sample1, sample2, value_str = parts
            # Remove "jaccard:" prefix if present
            value_str = value_str.replace("jaccard:", "")
            value = float(value_str)
            
            matrix[sample1][sample2] = value
            matrix[sample2][sample1] = value  # Ensure symmetry
            samples.update([sample1, sample2])
        except ValueError:
            print(f"Invalid number in line {line_num}: {value_str}", file=sys.stderr)
            continue

    if not samples:
        print("Error: No valid data processed", file=sys.stderr)
        sys.exit(1)

    samples = sorted(samples)
    print("\t" + "\t".join(samples))
    for s1 in samples:
        print("\t".join([s1] + [str(matrix[s1].get(s2, "NA")) for s2 in samples]))

if __name__ == "__main__":
    main()