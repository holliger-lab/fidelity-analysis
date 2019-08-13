#!/bin/bash

mkdir bc_split
echo "Splitting barcodes from fastq."
python fasta_demux.py -i ../../HS035_qfilt_q30_p90.fastq -b pol_barcodes.txt -o bc_split/.

