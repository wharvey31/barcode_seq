#!/bin/env python3
from __future__ import annotations
import argparse
import gzip
import os, sys
from pytrie import SortedStringTrie
import warnings
import pandas as pd


def opener(filename) -> object:
    # Reads in first two bytes of the input file and looks for gzip bytes
    # Returns gzip.open if found, else open
    f = open(filename, "rb")
    if f.read(2) == b"\x1f\x8b":
        f.close()
        return gzip.open
    else:
        f.close()
        return open


def read_input(input_file: str) -> tuple[dict, list]:
    # Reads the input file and returns a list and distionary of the reads
    # Check if file exists
    try:
        assert os.path.isfile(input_file)
    except:
        print(f"ERROR: {input_file} does not exist!")
        sys.exit(1)
    open_func = opener(input_file)
    read_counts = {}
    read_list = []
    with open_func(input_file, "rt") as file_in:
        for i, line in enumerate(file_in):
            if i % 4 == 1:
                if i not in read_counts:
                    read_counts[line.rstrip()] = {}
                    read_list.append(line.rstrip())
                else:
                    continue
    return read_counts, read_list


def find_anchor(cutoff: float, read_list: list) -> tuple[list, list, int, list]:
    # Detects where the anchor sequence starts by checking frequency of bases at each position
    base_count = [[], [], [], []]
    for i in range(len(read_list[0])):
        base_dict = {"A": 0, "G": 0, "T": 0, "C": 0}
        for read in read_list:
            base_dict[read[i]] += 1
        for j, base in enumerate(base_dict):
            base_count[j].append(base_dict[base] / len(read_list))
    # Computes value of most frequent base at each position
    max_val = [max([y[x] for y in base_count]) for x in range(len(base_count[0]))]
    # Finds where the anchor starts based on last value which is below proportional cutoff
    # Adds one to end so indexing is properly evaluated
    anchor_start = (
        max([x for x, val in enumerate(max_val) if val < args.anchor_cutoff]) + 1
    )
    # Truncates reads at anchor position
    map_reads = [x[0:anchor_start] for x in read_list]
    return base_count, max_val, anchor_start, map_reads


def find_hamming(
    map_reads: list, read_list: list, read_counts: dict, cells: int
) -> tuple[int, list, list]:
    # Creates a trie of all sequences
    trie = SortedStringTrie.fromkeys(map_reads)
    map_counts = {}
    # Iterates through potential sequences and looks for hamming-1 substiutions in trie
    for string in read_list:
        map_string = string[0:anchor_start]
        for i in range(len(map_string)):
            for letter in "GATC":
                if map_string[i] != letter:
                    neighbor = map_string[:i] + letter + map_string[i + 1 :]
                    # Adds dict entry for each sequence if hamming-1 found in trie
                    if neighbor in trie and map_string not in map_counts:
                        map_counts[map_string] = {}
                        map_counts[map_string][neighbor] = {}
                    elif neighbor in trie and neighbor not in map_counts[map_string]:
                        map_counts[map_string][neighbor] = {}
    # Computes hamming-1 degree of all reads
    hamming_degree = [len(map_counts[x]) for x in map_counts]
    # Zips and sorts potential sequences by their degree
    degree_zip = zip(map_counts, hamming_degree)
    degree_sorted = sorted(degree_zip, key=lambda x: x[1])[::-1]
    name_list, degree_list = zip(*degree_sorted)
    # Computes degree cutoff based on nyumber of sequenced cells
    max_degree = degree_list[100 - 1]
    return max_degree, name_list, degree_list


def write_output(output: str, degree_list: list, name_list: list, max_degree: int):
    # Writes to the output file
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as outfile:
        for i, degree in enumerate(degree_list):
            if degree >= max_degree:
                outfile.write(f"{name_list[i]}\n")
            else:
                break


def plot_degree_anchor(
    degree_list: list,
    cells: int,
    max_val: list,
    cutoff: int,
    base_count: list,
    anchor_start: int,
):
    try:
        import matplotlib.pyplot as plt
    except:
        warnings.warn("Unable to load matplotlib. Still exiting successfully")
        sys.exit(0)
    os.makedirs("plots", exist_ok=True)
    # Plot histogram with Hamming-1 degree of mapseq nodes
    fig, ax = plt.subplots()
    ax.hist(degree_list, bins=len(set(degree_list)), color="tab:blue", align="right")
    ax.axvline(degree_list[cells], linestyle="--", color="tab:red")
    ax.set_xlabel("Hamming-1 Degree")
    ax.set_ylabel("Barcodes")
    plt.title("Hamming-1 Neighbors")
    plt.savefig("plots/hamming_hist.png")
    # Plot base frequency
    fig, ax = plt.subplots()
    color_dict = {"A": "black", "G": "tab:red", "T": "tab:green", "C": "tab:blue"}
    for i, base in enumerate(color_dict):
        ax.plot(
            range(1, len(base_count[i]) + 1),
            base_count[i],
            color=color_dict[base],
            alpha=0.7,
            label=base,
        )
    ax.plot(range(1, len(max_val) + 1), max_val, color="tab:orange", linestyle="--")
    # Plot lines to indicate where anchor starts and anchor cutoff
    ax.axhline(cutoff, linestyle="--", color="tab:purple")
    ax.axvline(anchor_start + 1, linestyle="--", color="tab:purple")
    ax.set_xlabel("Read Position (bp)")
    ax.set_ylabel("Proportion of Base")
    ax.legend()
    plt.title("Sequence Content of Reads")
    plt.savefig("plots/read_content.png", bbox_inches="tight")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Python script to find Hamming-1 neighbors and degree of MapSeq Reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input FASTQ file name (can be gzipped or unzipped)",
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-c",
        "--cells",
        type=int,
        required=True,
        help="Number or cells sequenced",
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output file name",
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-a",
        "--anchor_cutoff",
        type=float,
        required=False,
        default=0.8,
        help="Frequency of most common base cutoff",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        required=False,
        help="Whether or not to produce anchor sequence plot and hamming degree plot",
    )

    args = parser.parse_args()

    print(f"Reading {args.input}...")
    read_counts, read_list = read_input(args.input)
    print(f"Detecting anchor position start with frequency: {args.anchor_cutoff}")
    base_count, max_val, anchor_start, map_reads = find_anchor(
        args.anchor_cutoff, read_list
    )
    print(f"...Detected anchor start position of {anchor_start+1} (1-based)")
    print(f"Searching for Hamming-1 neighbors")
    max_degree, name_list, degree_list = find_hamming(
        map_reads, read_list, read_counts, args.cells
    )
    print(f"Creating output: {args.output}")
    write_output(args.output, degree_list, name_list, max_degree)

    if args.plot:
        print(f"Plotting")
        plot_degree_anchor(
            degree_list,
            args.cells,
            max_val,
            args.anchor_cutoff,
            base_count,
            anchor_start,
        )
