#!/bin/env python3

import argparse
import gzip
import os, sys
import pandas as pd
import numpy as np
import tracemalloc
import matplotlib.pyplot as plt
from pytrie import SortedStringTrie


def opener(filename):
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'):
        f.close()
        return gzip.open
    else:
        f.close()
        return open


# parser = argparse.ArgumentParser() 

# parser.add_argument("--input", "-i", type=str, required=True)

# args = parser.parse_args()

# tracemalloc.start()

args = pd.Series(['/home/wharvey/Downloads/10_5_candidates.fastq.gz', 1, 10000], index=['input', 'hamming', 'cells'])

try:
    assert os.path.isfile(args.input)
except:
    print(f'ERROR: {args.input} does not exist!')
    sys.exit(1)

open_func = opener(args.input)

read_counts = {}
read_list = []

with open_func(args.input, 'rt') as file_in:
    for i, line in enumerate(file_in):
        if i % 4 == 1:
            if i not in read_counts:
                read_counts[line.rstrip()] = {}
                read_list.append(line.rstrip())
            else:
                continue

result = most_common_substring(read_list[0:100])

# degree_match = [0] * len(read_list)

# check_list = []


# distance = args.hamming

# for j, string in enumerate(read_counts):
#     print(j)
#     # Create a Trie containing all strings that are within the desired Hamming distance
#     trie_hamm = SortedStringTrie()
#     queue = deque([(string, 0)])
#     while queue:
#         string, d = queue.popleft()
#         if d == distance:
#             trie_hamm[string] = True
#         else:
#             for neighbor in generate_neighbors(string):
#                 queue.append((neighbor, d+1))
#     for key in trie_hamm.keys():
#         if key in trie and key != string:
#             read_counts[string][key] = {}


trie = SortedStringTrie.fromkeys(read_list)

# for i, string in enumerate(read_list[0:1]):
for j, string in enumerate(read_list):
    for i in range(len(string)):
        for letter in 'GATC':
            if string[i] != letter:
                neighbor = string[:i] + letter + string[i+1:]
                if neighbor in trie and neighbor not in read_counts[string]:
                    read_counts[string][neighbor] = {}


hamming_degree = [len(set(read_counts[x])) for x in read_list]

degree_zip = zip(read_list, hamming_degree)

degree_sorted = sorted(degree_zip, key = lambda x: x[1])[::-1]

name_list, degree_list = zip(*degree_sorted)





fig, ax = plt.subplots()

ax.hist(degree_list, bins=54, color='tab:blue')
ax.axvline(degree_list[args.cells], linestyle='--', color='tab:red')


plt.show()