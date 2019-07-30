
import csv
from subprocess import call

names = []
seqs = []

def VARNA_fig(name, row, struct):
    title = '{0} ({1})'.format(name.replace('/', ','), row[0])
    name = '{0}_{1}'.format(name.replace('/', ','), row[0])
    mature = row[3]
    star = row[5]
    precursor = row[10]
    seq = precursor
    if any(m.islower() for m in mature) or any(s.islower() for s in star):
        pass
    else:
        m_start = precursor.index(mature)
        m_end = m_start + len(mature)
        s_start = precursor.index(star)
        s_end = s_start + len(star)
        mir_place = round(((m_end - m_start)/2), 0)
        mir_s_place = round((s_end - s_start)/2)
        m_start = m_start + 1
        s_start = s_start + 1
        if m_start > s_start:
            call(['java', '-cp', input('/path/to/VARNA.jar'), 
                  'fr.orsay.lri.varna.applications.VARNAcmd', 
                  '-annotations', 'M:type=B,size=24,color=#000000,anchor='+str(int(mir_place) + m_start)
                  +';*:type=B,size=24,color=#000000,anchor='+str(int(mir_s_place) + s_start), 
                  '-title', title, '-titlesize', '12', 
                  '-sequenceDBN', seq, 
                  '-structureDBN', struct, 
                  '-highlightRegion', str(m_start)+'-'+str(m_end)+':radius=16,fill=#BCFFDD,outline=#6ED86E;'
                  +str(s_start)+'-'+str(s_end)+':radius=16,fill=#FF9999,outline=#FF3333;'
                  +str(s_end + 1)+'-'+str(m_start - 1)+':radius=16,fill=#FFF1A6,outline=#FFD700', 
                  '-o', name + '_plant-precursors_structure_RNAfold.png'])
        elif s_start > m_start:
            call(['java', '-cp', input('/path/to/VARNA.jar'), 
                  'fr.orsay.lri.varna.applications.VARNAcmd', 
                  '-annotations', 'M:type=B,size=24,color=#000000,anchor='+str(int(mir_place) + m_start)
                  +';*:type=B,size=24,color=#000000,anchor='+str(int(mir_s_place) + s_start), 
                  '-title', title, '-titlesize', '12', 
                  '-sequenceDBN', seq, 
                  '-structureDBN', struct, 
                  '-highlightRegion', str(m_start)+'-'+str(m_end)+':radius=16,fill=#BCFFDD,outline=#6ED86E;'
                  +str(s_start)+'-'+str(s_end)+':radius=16,fill=#FF9999,outline=#FF3333;'
                  +str(m_end + 1)+'-'+str(s_start - 1)+':radius=16,fill=#FFF1A6,outline=#FFD700', 
                  '-o', name + '_plant-precursors_structure_RNAfold.png'])
        else:
            pass
    
    return title

with open(input('/path/to/found_miRNAs.fasta'), 'r') as f:
    for name in f:
        name = name.rstrip('\n').lstrip('>')
        names.append(name)
        seq = next(f).rstrip('\n').replace('T', 'U')
        seqs.append(seq)

struct_names = []
structs = []
with open(input('/path/to/precursors_for_found_miRNAs.fass'), 'r') as f:
    # fass files are fasta files where the nucleotide sequence is followed by
    # the sequence for the structure of the precursor in dot/bracket notation
    # Structures can be obtained from ShortStack though the trimming rounds 
    # may not allow VARNA to match up all of the brackets.
    # Feed RNAfold the sequences is a good alternative. It will produce a fass
    # file as the stdout
    for name in f:
        name = name.rstrip('\n').lstrip('>')
        seq = next(f).rstrip('\n').replace('T', 'U')
        struct = next(f).rstrip('\n')
        struct_names.append(name)
        structs.append(struct)

rows = []
with open(input('/path/to/file_with_all_precursor_information.csv'), 'r') as f:
    csvreader = csv.reader(f, delimiter = ',')
    header = next(csvreader)
    for row in csvreader:
        mature = row[3]
        star = row[5]
        if mature in seqs:
            seqs_ind = seqs.index(mature)
            name = names[seqs_ind]            
            struct_name = '{0}_{1}'.format(name, row[0])
            struct_ind = struct_names.index(struct_name)
            struct = structs[struct_ind]
            title = VARNA_fig(name, row, struct)
        elif star in seqs:
            seqs_ind = seqs.index(star)
            name = names[seqs_ind]
            struct_name = '{0}_{1}'.format(name, row[0])
            struct_ind = struct_names.index(struct_name)
            struct = structs[struct_ind]
            title = VARNA_fig(name, row, struct)
        else:
            pass
        
