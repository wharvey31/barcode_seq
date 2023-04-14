# barcode_seq

Program for determining hamming-1 neighbors and degree of MapSeq reads


## Requirments
```
python >= 3.7
pytrie
Optional: matplotlib
```


## Usage

```
usage: find_barcodes.py [-h] -i INPUT -c CELLS -o OUTPUT [-a ANCHOR_CUTOFF] [-p]

Python script to find Hamming-1 neighbors and degree of MapSeq Reads

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTQ file name (can be gzipped or unzipped)
  -c CELLS, --cells CELLS
                        Number or cells sequenced
  -o OUTPUT, --output OUTPUT
                        Output file name
  -a ANCHOR_CUTOFF, --anchor_cutoff ANCHOR_CUTOFF
                        Frequency of most common base cutoff (default: 0.8)
  -p, --plot            Whether or not to produce anchor sequence plot and hamming degree plot (default: False)
```

## Optional Plots

### Hamming Histogram
```{figure}
plots/hamming_hist.png
```

### Read Content
```{figure}
plots/read_content.png
```