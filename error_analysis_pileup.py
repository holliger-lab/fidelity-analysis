import argparse
import numpy as np

parser = argparse.ArgumentParser(description='error_analysis: Calculates error rates and indels of HTS data from a pileup file.')
parser.add_argument('-pileup', help='pileup file. REQUIRED.', required=True)

args = parser.parse_args()



## Main ####
encode_dict = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    'N': 4,
    '+': 5,
    '-': 6,
}

count_array = []
total_error = []
total_matches = 0
total_mismatches = 0
total_indels = 0

ref_table = [[0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0]]

ref_sequence = ""

pileup = open(args.pileup, 'r')

header = "Position\tRef nt\tReads\tA\tC\tG\tT\tN\tIns\tDel\t\tA%\tC%\tG%\tT%\tN%\tIns%\tDel%\tError"

print(header)

for line in pileup:
    pileup_line = line.strip().split("\t")
    position = int(pileup_line[1])
    ref_nt = pileup_line[2]
    ref_sequence += ref_nt
    total_reads = int(pileup_line[3])
    read_data = pileup_line[4]

    count_match = read_data.count('.') + read_data.count(',')
    count_A = read_data.count('A') + read_data.count('a')
    count_C = read_data.count('C') + read_data.count('c')
    count_G = read_data.count('G') + read_data.count('g')
    count_T = read_data.count('T') + read_data.count('t')
    count_N = read_data.count('N') + read_data.count('n')
    count_ins = read_data.count('+')
    count_del = read_data.count('-')

    ## Populate the total matches/mismatches/indels variables
    total_matches += count_match
    total_mismatches += (count_A+count_C+count_G+count_T+count_N)
    total_indels += (count_ins+count_del)



    ## [A, C, G, T, N, +, -]
    position_counts = [count_A, count_C, count_G, count_T, count_N, count_ins, count_del]
    position_counts[encode_dict[ref_nt]] += count_match
    count_array.append(position_counts)

    # Calculate ref vs nt counts
    ## [A, C, G, T, N]
    position_counts_nts = [count_A, count_C, count_G, count_T, count_N]
    position_counts_nts[encode_dict[ref_nt]] += count_match
    ref_table[encode_dict[ref_nt]] = [sum(x) for x in zip(ref_table[encode_dict[ref_nt]], position_counts_nts)]

    ## Compute errors
    position_error_array = []

    for i in range(7):
        position_error_array.append(position_counts[i]/total_reads*100)

    position_total_error = (total_reads-count_match)/total_reads
    total_error.append(position_total_error)
    position_error_array.append(position_total_error)


    ## Print the position data
    print("%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.2E" % (
        position, ref_nt, total_reads,
        position_counts[0], position_counts[1], position_counts[2], position_counts[3],
        position_counts[4], position_counts[5], position_counts[6],
        position_error_array[0], position_error_array[1], position_error_array[2],
        position_error_array[3], position_error_array[4], position_error_array[5],
        position_error_array[6], position_total_error
        ))


print("\nMean total error\t%.2E" % np.mean(total_error))

## Calculate reference match vs reference mismatch
print("\n")
print('\tA\tC\tG\tT\tN')
print("ref A\t%d\t%d\t%d\t%d\t%d" % (ref_table[0][0], ref_table[0][1], ref_table[0][2], ref_table[0][3], ref_table[0][4]))
print("ref C\t%d\t%d\t%d\t%d\t%d" % (ref_table[1][0], ref_table[1][1], ref_table[1][2], ref_table[1][3], ref_table[1][4]))
print("ref G\t%d\t%d\t%d\t%d\t%d" % (ref_table[2][0], ref_table[2][1], ref_table[2][2], ref_table[2][3], ref_table[2][4]))
print("ref T\t%d\t%d\t%d\t%d\t%d" % (ref_table[3][0], ref_table[3][1], ref_table[3][2], ref_table[3][3], ref_table[3][4]))

print("\n")
print("Total matches\t%.2E" % total_matches)
print("Total mismatches\t%d" % total_mismatches)
print("Total indels\t%d" % total_indels)
print("Error rate\t%.2E" % ((total_mismatches+total_indels)/(total_mismatches+total_indels+total_matches)))

pileup.close()

