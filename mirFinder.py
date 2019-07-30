
import csv

m_s_dict = {}
nt_vars_set = set()
mirs = set()
reads_dict = {}
combined_reads_dict = {}
lowers = []

file = [input('/path/to/miRNA_precursors.csv')]

with open(file, 'r') as f:
    csvreader = csv.reader(f, delimiter = ',')
    for row in csvreader:
        lowers.append([row[3], row[0]])
        lowers.append([row[5], row[0]])
        
        cluster = row[0]
        per_map = row[1]
        m_len = row[2]
        mature = row[3].upper().replace('U', 'T')
        m_reads = row[4]
        star = row[5].upper().replace('U', 'T')
        s_reads = row[6]
        s_len = row[7]
        m_1nt = row[8].upper().replace('U', 'T')
        s_1nt = row[9].upper().replace('U', 'T')
        precursor = row[10]
        sec_struct = row[11]
        location = row[12]
        strand = row[13]
        sseqid = row[17]
        
        m_s_dict.setdefault(cluster, []).append(mature)
        m_s_dict.setdefault(cluster, []).append(star)
        m_s_dict.setdefault(cluster, []).append(m_1nt)
        m_s_dict.setdefault(cluster, []).append(s_1nt)
        m_s_dict.setdefault(cluster, []).append(precursor)
        m_s_dict.setdefault(cluster, []).append(sec_struct)
        m_s_dict.setdefault(cluster, []).append(sseqid)
        m_s_dict.setdefault(cluster, []).append(location)
        m_s_dict.setdefault(cluster, []).append(strand)
        nt_vars_set.add(m_1nt)
        nt_vars_set.add(s_1nt)
        mirs.add(mature)
        mirs.add(star)

nt_vars_string = ''.join(nt_vars_set)

#%%

reads_files = []
# Need to place the plant-only read files from each dataset here

G006_reads_dict = {}
G002_reads_dict = {}
BTIRed_reads_dict = {}

with open(reads_files[0], 'r') as f:
    csvreader = csv.reader(f, delimiter = '\t')
    for row in csvreader:
        header = row[0]
        read_num = header.lstrip('>')
        seq = next(csvreader)[0]
        if seq in mirs:
            G006_reads_dict.setdefault(seq, []).append(read_num)
        else:
            pass

with open(reads_files[1], 'r') as f:
    csvreader = csv.reader(f, delimiter = '\t')
    for row in csvreader:
        header = row[0]
        read_num = header.lstrip('>')
        seq = next(csvreader)[0]
        if seq in mirs:
            G002_reads_dict.setdefault(seq, []).append(read_num)
        else:
            pass

with open(reads_files[2], 'r') as f:
    csvreader = csv.reader(f, delimiter = '\t')
    for row in csvreader:
        header = row[0]
        read_num = header.lstrip('>')
        seq = next(csvreader)[0]
        if seq in mirs:
            BTIRed_reads_dict.setdefault(seq, []).append(read_num)
        else:
            pass

reads_dict = {}
for reads_file in reads_files:
    with open(reads_file, 'r') as f:
        csvreader = csv.reader(f, delimiter = '\t')
        for row in csvreader:
            header = row[0]
            read_num = header.lstrip('>')
            seq = next(csvreader)[0]
            if seq in mirs:
                reads_dict.setdefault(seq, []).append(read_num)
            else:
                pass

#%%

total = 0
lines = []
for k in reads_dict.keys():
    try:
        G006_len = len(G006_reads_dict[k])
    except KeyError:
        G006_len = 0
    try:
        G002_len = len(G002_reads_dict[k])
    except KeyError:
        G002_len = 0
    try:
        BTIRed_len = len(BTIRed_reads_dict[k])
    except KeyError:
        BTIRed_len = 0
    total += len(reads_dict[k])
    line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(k, len(k), G006_len, G002_len, 
                                                 BTIRed_len, 
                                                 len(reads_dict[k]))
    lines.append(line)

# lines now contains a printable display of all miRNAs found in your datasets
# you may need to generate fasta files and csvs for futher analysis

#%%

print_lines = []
with open(input('/path/to/file_with_miRNAs_found_in_your_datasets.fasta'), \
          'r') as f:
    file_lines = f.readlines()

    for line in lines:
        line_seq = line.split('\t')[0]
        len_pls = len(print_lines)
        for i in range(0, len(file_lines), 2):
            name = file_lines[i]
            seq = file_lines[i + 1].rstrip('\n')
            name = name.lstrip('>').rstrip('\n')
            if line_seq == seq:
                line = "{0}\t{1}".format(name, line)
                print_lines.append(line)
        if len(print_lines) == len_pls:
            name = 'Not in Gut'
            line = "{0}\t{1}".format(name, line)
            print_lines.append(line)

#%%

prev_mirs = {}
with open(input('/path/to/file_with_miRNAs_to_check_yours_against.csv'), \
          'r') as f:
    # can be used when you collect miRNA sequences from other papers to see
    # if any show up in your dataset
    csvreader = csv.reader(f, delimiter = ',')
    for row in csvreader:
        p_name = row[0]
        p_seq = row[1]
        p_mir = row[2]
        prev_mirs.setdefault(p_name, []).append([p_mir, p_seq])

species_s = sorted(prev_mirs.keys())
i = 0
for line in print_lines:
    seq = line.split('\t')[1]
    found_list = [''] * len(species_s)
    for species in species_s:
        species_list = list(zip(*prev_mirs[species]))
        if seq in species_list[1]:
            sp_ind = species_s.index(species)
            seq_ind = species_list[1].index(seq)
            found_list[sp_ind] = species_list[0][seq_ind]
        else:
            pass
        print_lines[i] = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(line, 
                         found_list[0], found_list[1], found_list[2],
                         found_list[3], found_list[4], found_list[5])
    i += 1

with open('mirFinder_Bac.out', 'w') as f:
    f.write("Homolog\tmiRNA Sequence\tLength\tExact Matches In: G006\tG002\tUS"
            "DA\tTotal\tPreviously Found In: {0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"\
            .format(species_s[0], species_s[1], species_s[2], species_s[3], 
                    species_s[4], species_s[5]))
    for line in print_lines:
        f.write(line + '\n')
        print(line)
    print('Total: {0} reads'.format(total))
    f.write('Total: {0} reads\n'.format(total))
