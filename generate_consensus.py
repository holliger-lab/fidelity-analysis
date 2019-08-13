import re, sys, os
import argparse
from collections import defaultdict

## Script configuration
pre_NNN_bc = ["GTGTGCCTGG", 15, 35] ## Sequence of a constant region just prior to the NNN barcode. Numbers are roughly where to find this sequence.
pre_template_bc = ["TGGACTGTGG", 57, 80] ## Sequence of a constant region just prior to the template. Numbers are roughly where to find this sequence.
NNN_length = 25 ## Length of the NNN barcode
template_length = 19 ## Length of the template within this read. This must also match the length of the template.
n_min_reads = 3 ## Minimum number of reads required for each unique barcode.
consensus_threshold = 0.66 ## For each unique barcode, sequencing errors need to be removed by consensus. This is the threshold for determing the correct base.


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

def calcConsensus(seqList, threshold=0.66):
    consensus = ''
    con_len = len(seqList[0])


    ## Move through the sequences, one base at a time.
    for i in range(con_len):

        nt_dict = {}
        nt_count = 0

        for entry in seqList:
            # Move through all of the sequences in the seqList.
            if entry[i] not in nt_dict:
                nt_dict[entry[i]] = 1
            else:
                nt_dict[entry[i]] += 1

            nt_count += 1

        max_nt = []
        max_count = 0

        for nt in nt_dict:
            if nt_dict[nt] > max_count:
                max_nt = [nt]
                max_count = nt_dict[nt]
            elif nt_dict[nt] == max_count:
                max_nt.append(nt)

        if (float(max_count) / float(nt_count)) >= threshold:
            consensus += max_nt[0]
        else:
            ## Consensus is ambiguous, abort this list.
            return None

    return consensus



parser = argparse.ArgumentParser(description='process_fidelity_data.py: Process split barcode fasta files to write consensus sequences.')
parser.add_argument('-i', help='Input fasta file.', required=True)
parser.add_argument('-o', help='Output file', required=True)

args = parser.parse_args()


bc_split_file = args.i
consensus_out_file = args.o

sequence_dict = defaultdict(list)
fasta_data = readFasta(bc_split_file)


counter = 0
for read in fasta_data:

    seq = read[1]

    NNN_seq = ""
    template_seq = ""

    ## Extract NNN
    pattern_NNN = r"(" + pre_NNN_bc[0] + ")"
    pat_match_NNN = re.compile(pattern_NNN)

    search_NNN = re.search(pat_match_NNN, seq[pre_NNN_bc[1]:pre_NNN_bc[2]])

    if search_NNN:
        NNN_start = pre_NNN_bc[1]+search_NNN.span()[1]
        NNN_seq = seq[NNN_start:NNN_start+NNN_length]

    ## Extract template bits
    pattern_template = r"(" + pre_template_bc[0] + ")"
    pat_match_template = re.compile(pattern_template)

    search_template = re.search(pat_match_template, seq[pre_template_bc[1]:pre_template_bc[2]])

    if search_template:
        template_start = pre_template_bc[1]+search_template.span()[1]
        template_seq = seq[template_start:template_start+template_length]


    if NNN_seq != "" and template_seq != "":
        if len(template_seq) != template_length:
            continue

        if NNN_seq not in sequence_dict:
            sequence_dict[NNN_seq] = [template_seq]
        else:
            sequence_dict[NNN_seq].append(template_seq)

        counter += 1



fh = open(consensus_out_file, 'w')

counter_all = 0
counter_three = 0
for barcode in sequence_dict.items():
    if len(barcode[1]) >= n_min_reads:
        barcode_seq = barcode[0]
        seq_list = barcode[1]

        consensus_seq = calcConsensus(seq_list, consensus_threshold)
        if consensus_seq == None:
            continue

        ## Correct all reads with the same barcode.
        # for i in range(len(seq_list)):
        #     fh.write(">%d:%s\n" % (counter_all, barcode_seq))
        #     fh.write(consensus_seq + "\n")
        #     counter_all += 1

        # Write out the corrected consensus sequence for this barcode.
        fh.write(">%d:%s\n" % (counter_all, barcode_seq))
        fh.write(consensus_seq + "\n")
        counter_all += 1

        counter_three += 1

fh.close()

print("Total unique: %d\n" % len(sequence_dict))
print("Total unique w/ at least %d reads: %d\n" % (n_min_reads, counter_three))



