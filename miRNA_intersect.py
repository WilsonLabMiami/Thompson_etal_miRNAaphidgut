
import csv

strand_list = []
miranda_list = []
pita_list = []
hybrid_list = []
total_list = []
m_dict = {}
p_dict = {}
h_dict = {}

# check below for necessary input file formating or adjust to your format
with open('exact_mirs_targets.out.1l_parsed.su') as f:  # miranda
    csvreader = csv.reader(f, delimiter = '\t')
    for row in csvreader:
        start = (int(row[5].split(' ')[0]))
        end = (int(row[5].split(' ')[1]))
        mirna = row[0].lstrip('>')
        target = row[1]
        miranda_list.append([mirna+'\t'+target, set(list(range(start, end)))])
        m_dict.setdefault(mirna+'\t'+target, []).append(set(list(range(start, end))))
        

with open('exact_mirs_targets_pita_results.sort.uniq') as f:  # pita
    csvreader = csv.reader(f, delimiter = '\t')
    header = next(csvreader)
    for row in csvreader:
        mirna = row[1].split(' ')[0]
        target = row[0]
        start = (int(row[3]))
        end = (int(row[2]))
        pita_list.append([mirna+'\t'+target, set(list(range(start, end)))])
        p_dict.setdefault(mirna+'\t'+target, []).append(set(list(range(start, end))))

with open('exact_mirs_RNAhybrid_cat.out') as f:  # RNAhybrid
    csvreader = csv.reader(f, delimiter = ':')
    header = next(csvreader)
    for row in csvreader:
        mirna = row[3]
        target = ':'.join(row[0:2])
        start = int(row[7])
        end = (start + len(row[8]))
        hybrid_list.append([mirna+'\t'+target, set(list(range(start, end)))])
        h_dict.setdefault(mirna+'\t'+target, []).append(set(list(range(start, end))))

m_set = set()
p_set = set()
h_set = set()
for m in miranda_list:
    m_set.add(m[0])
for p in pita_list:
    p_set.add(p[0])
for h in hybrid_list:
    h_set.add(h[0])
total_set = m_set & p_set & h_set

for entry in total_set:
    if len(m_dict[entry]) > 1 or len(p_dict[entry]) > 1 or len(h_dict[entry]) > 1:
        b = set()
        for value in m_dict[entry]:
            b = b | value
            m_dict[entry] = b
        b = set()
        for value in p_dict[entry]:
            b = b | value
            p_dict[entry] = b
        b = set()
        for value in h_dict[entry]:
            b = b | value
            h_dict[entry] = b
        if len(set(m_dict[entry] & p_dict[entry] & h_dict[entry])) > 0:
            seed = sorted(set(m_dict[entry] & p_dict[entry] & h_dict[entry]))
            if (1 + seed[-1] - seed[0]) > 8:
                print('ERROR: Fix entry below manually')
                # multiple seeds overlap and must be manualy split apart in
                # the output file
                print(seed)
                print(entry+'\t'+str(seed[0])+'\t'+str(seed[-1]))
            else:
                print(entry+'\t'+str(seed[0])+'\t'+str(seed[-1]))
    else:
        if len(set(m_dict[entry][0] & p_dict[entry][0] & h_dict[entry][0])) > 0:
            seed = sorted(set(m_dict[entry][0] & p_dict[entry][0] & h_dict[entry][0]))
            print(entry+'\t'+str(seed[0])+'\t'+str(seed[-1]))


