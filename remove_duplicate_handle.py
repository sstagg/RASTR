#!/usr/bin/env python
import cupy as cp
import pandas as pd
from starparse import StarFile
import sys

def process_chunk(chunk, repeats):
    for _, row in chunk.iterrows():
        image_number = int(row.iloc[0]) + 1
        if image_number in repeats:
            continue

        correlation_values = cp.array(row.iloc[1:])
        mean_correlation = correlation_values.mean()
        rmsd = cp.sqrt(cp.mean((correlation_values - mean_correlation)**2))
        threshold = mean_correlation + 10 * rmsd

        indexes = cp.where(correlation_values > threshold)[0].get() + 1
        new_repeats = [int(i) for i in indexes if int(i) > image_number and int(i) not in repeats]
        
        for i in new_repeats:
            print(image_number, i)
        repeats.update(new_repeats)

        print(image_number, end='\r')

def main():
    if len(sys.argv) != 4:
        print("Usage: script.py <csv_file> <star_file> <output_file>")
        return

    chunk_size = 1000  # Adjust based on your system's memory
    repeats = set()

    # Process CSV in chunks
    for chunk in pd.read_csv(sys.argv[1], chunksize=chunk_size):
        process_chunk(chunk, repeats)

    # Save results
    pd.DataFrame(list(repeats)).to_csv(sys.argv[3], index=False)
    print(f"\nTotal repeats: {len(repeats)}")

if __name__ == '__main__':
    main()
