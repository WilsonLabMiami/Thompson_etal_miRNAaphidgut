
# USAGE: python3 ./elimMatched.py [APHID SAM FILE] [BUCHNERA SAM FILE] \
#                                 [PLANT SAM FILE] [READS TO BE FILETERED] \
#                                 [DESIRED OUTPUT READS FILE]

    
from sys import argv
import csv

matched_dict = {}
genomes_dict = {}
aphid = set()
buchnera = set()
plant = set()


primary = ['Buchnera', 'Myzus']
secondary = ['plants', 'other_bacteria', 'viruses']
genomes_dict = {}
genomes_dict_collapsed = {}



with open(argv[1], newline='') as f:
    csvreader = csv.reader(f, delimiter = '\t')
    for row in csvreader:
        if (len(row) >= 14) and (('XM:i:0' or 'XM:i:1') in str(row)):
            readnum = row[0]
            refgen = row[2]
            seq = row[9]
            refgen = 'aphid'
            aphid.add(readnum)
            matched_dict.setdefault(readnum, []).append(refgen)
        else:
            pass

with open(argv[2], newline='') as f:
    csvreader = csv.reader(f, delimiter = '\t')
    for row in csvreader:
        if (len(row) >= 14) and (('XM:i:0' or 'XM:i:1') in str(row)):
            readnum = row[0]
            refgen = row[2]
            seq = row[9]
            refgen = 'buchnera'
            buchnera.add(readnum)
            matched_dict.setdefault(readnum, []).append(refgen)
        else:
            pass

with open(argv[3], newline='') as f:
    csvreader = csv.reader(f, delimiter = '\t')
    for row in csvreader:
        if (len(row) >= 14) and (('XM:i:0' or 'XM:i:1') in str(row)):
            readnum = row[0]
            refgen = row[2]
            seq = row[9]
            refgen = 'plant'
            plant.add(readnum)
            matched_dict.setdefault(readnum, []).append(refgen)
        else:
            pass

mapped_to_a_b_p = aphid | buchnera | plant

reads_in = open(argv[4], 'r')
# Full reads file
reads_out = open(argv[5], 'w')
# File that will contain only unmatched reads
for header_line in reads_in:
    # each header_line should be a FASTA header with format '>#'
    seq_line = next(reads_in)
    # seq_line should be FASTA seq
    readnum = header_line[1:-1]
    # removes '>' from begining of header line leaving read number

    if readnum not in matched_dict:
        reads_out.write(header_line)
        reads_out.write(seq_line)
        # writes to reads_out in FASTA format for unmatched entries in reads_in
    else:
        pass

reads_in.close()
reads_out.close()
