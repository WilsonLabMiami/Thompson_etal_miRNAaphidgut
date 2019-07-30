from os import listdir
import csv
import re
import matplotlib.pyplot as plt

map_folder = input('/path/to/bamAnalysis_output_alignment_files.txt')

txt_files = []
map_folder_files = listdir(map_folder)
for file in map_folder_files:
    if file.startswith('bol') and file.endswith('_reads.txt') and file not in ['summary_reads.txt']:
        txt_files.append(file)
    else:
        pass

for file in txt_files:
    with open(map_folder + file) as f:
        csvreader = csv.reader(f, delimiter = ' ')
        name = next(csvreader)[0].replace(',', '/')
        precursor = ref_seq = next(csvreader)[0].replace('U', 'T')
        mature = re.split(r'([AGCTU]+)', next(csvreader)[0])[1]
        star = re.split(r'([AGCTU]+)', next(csvreader)[0])[1]
        empty = next(csvreader)
        
        m_start = precursor.index(mature) + 1
        m_end = m_start + len(mature) - 1
        s_start = precursor.index(star) + 1
        s_end = s_start + len(star) - 1
        end = len(precursor)
        
        read_seqs = []
        read_counts = []
        reads_expanded = []
        mat_count = 0
        for row in csvreader:
            if row[0][0] == '.':
                read_seq = re.split(r'([AGCTU]+)', row[0])[1]
                read_count = int(row[2].lstrip('a='))
                start_ind = ref_seq.index(read_seq) + 1
                end_ind = start_ind + len(read_seq) - 1
                if start_ind == m_start and end_ind == m_end and read_count > 300:
                    read_count = 1
                else:
                    pass
                read_seqs.append(read_seq)
                read_counts.append(read_count)
                reads_expanded += [read_seq] * read_count

#%%
        inds = []
        inds_2 = []
        lines = []
        a_ind = 0
        a = []
        for read in reads_expanded:
            start_ind = ref_seq.index(read) + 1
            end_ind = start_ind + len(read) - 1
            inds.append([start_ind, end_ind, read])
        len_inds = len(inds)
        while len(inds_2) < len_inds:
            line = []
            split_list = list(zip(*inds))

            first_line_entry = inds[split_list[1].index(min(split_list[1]))]
            line.append(first_line_entry)
            last_entry_index = inds.index(line[-1])
            inds = inds[:last_entry_index] + inds[last_entry_index + 1:]
            split_list = list(zip(*inds))
            if len(split_list) == 0:
                break
            else:
                pass
            while any(split_list[0][j] > (line[-1][1] + 1) for j in range(0, len(split_list[0]))):
                a = []
                for entry in inds:
                    if entry[0] > (1 + line[-1][1]):
                        a.append(entry)
                    else:
                        pass
                if len(a) > 0:
                    a = list(zip(*a))
                    a_ind = a[1].index(min(a[1]))
                    line_entry = [a[0][a_ind], a[1][a_ind], a[2][a_ind]]
                    line.append(line_entry)
                    last_entry_index = inds.index(line[-1])
                    inds = inds[:last_entry_index] + inds[last_entry_index + 1:]
                    split_list = list(zip(*inds))
                else:
                    pass
                if len(split_list) == 0:
                    break
                else:
                    pass

            for entry in line:
                inds_2.append(entry)
            lines.append(line)
            if len(split_list) == 0:
                break
            else:
                pass
            
#%%
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        i = mat_count = 0
        darkdarkgreen = '#015b00'
        lw = 2
        for line in lines:
            i += 1
            for entry in line:
                x = list(range(entry[0], entry[1] + 1))
                y = [i] * len(x)
                if entry[0] == m_start and entry[1] == m_end:
                    plt.plot(x, y, color = darkdarkgreen, linewidth = lw)
                    mat_count += 1
                elif entry[0] == s_start and entry[1] == s_end:
                    plt.plot(x, y, color = 'darkred', linewidth = lw)
                else:
                    plt.plot(x, y, color = 'black', linewidth = lw)

#%%
    alpha = 0.25
    ymin = 0
    if m_start > s_start:
        plt.axvspan(0.5, s_start - 0.5, color = 'black', 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(s_start - 0.5, s_end + 0.5, color = 'firebrick', 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(s_end + 0.5, m_start - 0.5, color = 'gold', 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(m_start - 0.5, m_end + 0.5, color = darkdarkgreen, 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(m_end + 0.5, end + 0.5, color = 'black', 
                    alpha = alpha, zorder = 1, ymin = ymin)
    elif s_start > m_start:
        plt.axvspan(0.5, m_start - 0.5, color = 'black', 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(m_start - 0.5, m_end + 0.5, color = darkdarkgreen, 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(m_end + 0.5, s_start - 0.5, color = 'gold', 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(s_start - 0.5, s_end + 0.5, color = 'firebrick', 
                    alpha = alpha, zorder = 1, ymin = ymin)
        plt.axvspan(s_end + 0.5, end + 0.5, color = 'black', 
                    alpha = alpha, zorder = 1, ymin = ymin)
    plt.xlim(1, len(ref_seq))
    plt.ylim(0.25, i + 0.75)
    plt.yticks([], [])
    plt.xticks(range(0, len(precursor), 25)[1:], 
               range(0, len(precursor), 25)[1:], fontsize = 8)
    ax.spines['right'].set_visible(False)  # Removes right axis
    ax.spines['left'].set_visible(False)  # Removes left axis
    ax.spines['top'].set_visible(False)  # Removes right axis
    ax.yaxis.set_ticks_position('none')  # Keeps vertical ticks hidden
    plt.xlabel('\nPosition in Precursor')
    plt.title('{0}\n(Cluster_{1})'.format(name.split('_Cluster_')[0],
                                          name.split('_Cluster_')[1]))
    
    name = name.replace('/', ',')
    
    fig.set_size_inches(3.25, (0.075 * i))
    try:
        if len(reads_expanded) > 50:
            plt.savefig(map_folder + '{0}_individual-reads.svg'.format(name), 
                        bbox_inches = 'tight', format = 'svg', dpi = 300)
        else:
            plt.savefig(map_folder + '{0}_individual-reads.svg'.format(name), 
                        bbox_inches = 'tight', format = 'svg', dpi = 300)
    except RuntimeError:
        print('RuntimeError for', file)
    plt.show()
    plt.close('all')