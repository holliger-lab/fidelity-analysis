### Fasta_demux.py
### Searches through HTS data and extracts barcoded sequences.

import re, sys
import argparse

parser = argparse.ArgumentParser(description='fasta_demux: Searches through HTS fasta file and extracts sequences containing one or more barcodes.')
parser.add_argument('-i', help='Input fasta file.', type=argparse.FileType('r'), required=True)
parser.add_argument('-b', help='Barcodes file.', type=argparse.FileType('r'), required=True)
parser.add_argument('-o', help='Output location', required=True)
parser.add_argument('-bc', help='Remove the barcode from the output sequence.', default=False, action='store_true')

args = parser.parse_args()


def readFasta(fastaFile):
    fh = open(fastaFile, 'r')
    for line in fh:
        header = ""
        seq = ""
        if line[0] == '>':
            header = line.rstrip()[1:]
            if sys.version_info[0] < 3:
                seq = fh.next().rstrip()
            else:
                seq = fh.readline().rstrip()
        yield [header, seq]
    fh.close()


def readFastq(fastqFile, tile_id=None):
    if fastqFile.split('.')[-1] == "gz":
        fh = gzip.open(fastqFile, 'rt')
    else:
        fh = open(fastqFile, 'r')

    for line in fh:
        if line[0] == '@':
            header = line.rstrip()
            seq = ''
            if sys.version_info[0] < 3:
                line = fh.next()
            else:
                line = fh.readline()
            while line[0] != '+':
                seq += line.rstrip()
                if sys.version_info[0] < 3:
                    line = fh.next()
                else:
                    line = fh.readline()
            if tile_id == None:
                yield [header, seq]
            else:
                header_split = header.split(':')
                lane = int(header_split[3])
                tile_id_fasta = "%d_%s" % (lane, header_split[4])
                if tile_id == tile_id_fasta:
                    yield [header, seq]
            if sys.version_info[0] < 3:
                line = fh.next()
            else:
                fh.readline()
    fh.close()


def readFile(file_path):
    if file_path.split('.')[-1] == "fasta":
        return readFasta(file_path)
    elif file_path.split('.')[-1] == "fastq":
        return readFastq(file_path)

## Parse the barcodes file.
barcodes = []
for barcode in args.b:
    if barcode[0] != '#':
        barcodes.append(barcode.split())

# Open file handles for each barcode and append to bc_fh list.
bc_fh = []
for barcode in barcodes:
    output_location = args.o
    output_fasta = barcode[0] + ".fasta"
    fh = open(output_location + "/" + output_fasta, 'w')
    bc_fh.append(fh)


total_count = 0
bc_count = [0]*len(barcodes)
# Move through fasta file
for entry in readFile(args.i.name):

    # For each barcode, find in range.
    for i, barcode in enumerate(barcodes):
        barcode_seq = barcode[1]
        start = int(barcode[2])
        end = int(barcode[3])
        seq = entry[1][start:end]
        # print(barcode)

        #process the barcode for N or W
        barcode_seq_proc = ""
        for nt in barcode_seq:
            if nt == 'N':
                n_new = "(A|T|C|G)"
                barcode_seq_proc += n_new
            elif nt == 'W':
                n_new = "(A|T)"
                barcode_seq_proc += n_new
            elif nt == 'S':
                n_new = "(C|G)"
                barcode_seq_proc += n_new
            else:
                barcode_seq_proc += nt

        pattern = r"(" + barcode_seq_proc + ")"

        search = re.search(pattern, seq)
        if search:
            # If here, we have a match. Write to file.
            if args.bc:
                end_idx = search.span()[1]
                seq_mod = entry[1][end_idx:]
                bc_fh[i].write(">"+entry[0]+"\n")
                bc_fh[i].write(seq_mod+"\n")
            else:
                bc_fh[i].write(">"+entry[0]+"\n")
                bc_fh[i].write(entry[1]+"\n")
            bc_count[i] += 1

    total_count += 1


print("Read %d sequences." % total_count)
print("Wrote:")
for i, barcode in enumerate(barcodes):
    print("%s  %i" % (barcode[0], bc_count[i]))
    bc_fh[i].close()


