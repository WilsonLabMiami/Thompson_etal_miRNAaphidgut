import csv
from Bio import SeqIO
import re
import pysam  # only Mac or Linux
import numpy as np
import matplotlib.pyplot as plt

docs_folder = input('/path/to/bam_alignment_files_and_indexes')
map_folder = input('/path/to/folder_for_alignment_visualizations')

m_s_dict = {}
with open(docs_folder + input('/path/to/precursor_info_file_only_including_'\
                              'those_with_miRNAs_in_your_dataset'), 'r') as f:
    # check below for file formatting needs or modify below to match your file
    csvreader = csv.reader(f, delimiter = ',')
    for row in csvreader:
        name = row[0]
        mature = row[3]
        star = row[5]
        precursor = row[10]
        m_s_dict[name] = {'mature' : mature, 'star' : star, 
                          'precursor' : precursor}

record_dict = SeqIO.to_dict(SeqIO.parse(docs_folder + 
                                        input('/path/to/precursors.fasta'), 
                                        'fasta'))

bam_ext = '_Gut_plant-precursors_bowtie1_a_T_mm0_gap0.sorted.bam'
# file extension of bam files
BTIRed_bam = pysam.AlignmentFile(docs_folder + 'BTIRed' + bam_ext)
G002_bam = pysam.AlignmentFile(docs_folder + 'G002' + bam_ext)
G006_bam = pysam.AlignmentFile(docs_folder + 'G006' + bam_ext)

bams = [G006_bam, G002_bam, BTIRed_bam]

G006_not_m_s = open(docs_folder + 'G006_not_m_s_reads.fa', 'w')
# will store all reads that are not assigned to the mature or star sequences
# for later analysis
G002_not_m_s = open(docs_folder + 'G002_not_m_s_reads.fa', 'w')
BTIRed_not_m_s = open(docs_folder + 'BTIRed_not_m_s_reads.fa', 'w')

matures_summary = stars_summary = mature_tails_summary = star_tails_summary = \
loops_summary = mat_tail_overs_summary = star_tail_overs_summary = \
loop_mat_overs_summary = loop_star_overs_summary = others = 0

mat_count_sum = star_count_sum = loop_count_sum = mat_tail_count_sum = \
star_tail_count_sum = i_sum = 0

for name in record_dict.keys():
    matures = stars = mature_tails = star_tails = loops = mat_tail_overs = \
    star_tail_overs = loop_mat_overs = loop_star_overs = 0
    mature = m_s_dict[name]['mature']
    star = m_s_dict[name]['star']
    precursor = m_s_dict[name]['precursor']
    if any(m.islower() for m in mature) or any(s.islower() for s in star):
        pass
    else:
        if '/' in name:
            name = name.replace('/', ',')
        else:
            pass
        with open(map_folder + '{0}_reads.txt'.format(name), 'w') as f:
            f.write('{0}\n'.format(name))
            f.write('{0} Precursor\n'.format(precursor))
            m_start = precursor.index(mature)
            m_end = m_start + len(mature)
            s_start = precursor.index(star)
            s_end = s_start + len(star)
            end = len(precursor)
            
            name = name.replace(',', '/')
            read_counts = [0] * len(record_dict[name].seq)
            positions = list(range(1, len(read_counts) + 1))
            
            f.write('{0}{1}{2} Mature\n'.format('-' * m_start, 
                                            mature.replace('U', 'T'), 
                                            '-' * (end - m_end)))
            f.write('{0}{1}{2} Star\n\n'.format('-' * s_start, 
                                             star.replace('U', 'T'), 
                                             '-' * (end - s_end)))
            
            read_lines = []
            for bam in bams:
                aphid_line = str(bam.filename).lstrip("b'" + docs_folder).split(bam_ext + "'")[0]
                for read in bam:
                    read_ref_name = read.reference_name.split('_Cluster_')[-1]
                    name_name = name.split('_Cluster_')[-1]
                    if read_ref_name == name_name:
                        readnum = read.qname
                        read_positions = read.positions
                        record = record_dict[name]
                        read_seq = read.seq
                        ref_seq = record.seq
                        read_line = ''
                        for i in range(0, len(ref_seq)):
                            if i in read_positions:
                                read_line += str(ref_seq[i])
                                read_counts[i] += 1
                            else:
                                read_line += '.'
                        for pos in range(0, len(read_line)):
                            if read_line[pos] == '.' and \
                            ((m_start <= pos < m_end) or \
                             (s_start <= pos < s_end)):
                                read_line = read_line[:pos] + '_' + \
                                            read_line[pos+1:]
                        read_lines.append(read_line)

                        mature_positions = positions[m_start : m_end + 1]
                        star_positions = positions[s_start : s_end + 1]
                        if m_start > s_start:
                            loop_positions = positions[s_end + 1 : m_start]
                            mature_tail_positions = positions[m_end + 1 :]
                            star_tail_positions = positions[: s_start]
                        elif s_start > m_start:
                            loop_positions = positions[m_end + 1 : s_start]
                            star_tail_positions = positions[s_end + 1:]
                            mature_tail_positions = positions[: m_start]
                        else:
                            pass
                        i = 0
                        mat_count = star_count = loop_count = mat_tail_count =\
                        star_tail_count = 0
                        for pos in read_positions:
                            i += 1
                            i_sum += 1
                            if pos in mature_positions:
                                mat_count += 1
                                mat_count_sum += 1
                            elif pos in star_positions:
                                star_count += 1
                                star_count_sum += 1
                            elif pos in mature_tail_positions:
                                mat_tail_count += 1
                                mat_tail_count_sum += 1
                            elif pos in star_tail_positions:
                                star_tail_count += 1
                                star_tail_count_sum += 1
                            elif pos in loop_positions:
                                loop_count += 1
                                loop_count_sum += 1
                            else:
                                pass
                        if (mat_tail_count <= 1 and loop_count <= 1) and mat_count > 0:
                            matures += 1
                            matures_summary += 1
                        elif (star_tail_count <= 1 and loop_count <= 1) and star_count > 0:
                            stars += 1
                            stars_summary += 1
                        elif i - mat_tail_count <= 1 and mat_tail_count > 0:
                            mature_tails += 1
                            mature_tails_summary += 1
                            if aphid_line == 'G006':
                                G006_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'G002':
                                G002_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'BTIRed':
                                BTIRed_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            else:
                                pass
                        elif i - star_tail_count <= 1 and star_tail_count > 0:
                            star_tails += 1
                            star_tails_summary += 1
                            if aphid_line == 'G006':
                                G006_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'G002':
                                G002_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'BTIRed':
                                BTIRed_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            else:
                                pass
                        elif (star_count <= 1 and mat_count <= 1) and loop_count > 0:
                            loops += 1
                            loops_summary += 1
                            if aphid_line == 'G006':
                                G006_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'G002':
                                G002_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'BTIRed':
                                BTIRed_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            else:
                                pass
                        else:
                            others += 1
                            if aphid_line == 'G006':
                                G006_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'G002':
                                G002_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            elif aphid_line == 'BTIRed':
                                BTIRed_not_m_s.write(">{0}\n{1}\n".format(readnum, read_seq))
                            else:
                                pass
# Below matplotlib section can be used to generat histogram visualizations
# of read mapping for each precursor.

#            fig = plt.figure(figsize = (15, 7))
#            ax = plt.subplot(1, 1, 1)
#            barlist = plt.bar(positions, read_counts, width = 1.0, 
#                              color = 'black', zorder = 5)
#            
#            for bar in barlist[m_start:m_end]:
#                bar.set_color('darkgreen')
#            for bar in barlist[s_start:s_end]:
#                bar.set_color('firebrick')
#            if m_start > s_start:
#                for bar in barlist[s_end:m_start]:
#                    bar.set_color('gold')
#                plt.axvspan(0, s_start + 0.5, color = 'black', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(s_start + 0.5, s_end + 0.5, color = 'firebrick', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(s_end + 0.5, m_start + 0.5, color = 'gold', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(m_start + 0.5, m_end + 0.5, color = 'darkgreen', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(m_end + 0.5, end, color = 'black', 
#                            alpha = 0.25, zorder = 1)
#            elif s_start > m_start:
#                for bar in barlist[m_end:s_start]:
#                    bar.set_color('gold')
#                plt.axvspan(0, m_start + 0.5, color = 'black', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(m_start + 0.5, m_end + 0.5, color = 'darkgreen', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(m_end + 0.5, s_start + 0.5, color = 'gold', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(s_start + 0.5, s_end + 0.5, color = 'firebrick', 
#                            alpha = 0.25, zorder = 1)
#                plt.axvspan(s_end + 0.5, end, color = 'black', 
#                            alpha = 0.25, zorder = 1)
#            else:
#                pass
#            if max(read_counts) <= 5:
#                ind = np.arange(0, max(read_counts) + 1, 1)
#                plt.yticks(ind[1:], ind[1:])
#            elif 5 < max(read_counts) <= 10:
#                ind = np.arange(0, max(read_counts) + 1, 2)
#                plt.yticks(ind[1:], ind[1:])
#            elif 10 < max(read_counts) <= 50:
#                ind = np.arange(0, max(read_counts) + 1, 5)
#                plt.yticks(ind[1:], ind[1:])
#            elif 50 < max(read_counts) <= 150:
#                ind = np.arange(0, max(read_counts) + 1, 10)
#                plt.yticks(ind[1:], ind[1:])
#            elif 150 < max(read_counts):
#                ind = np.arange(0, max(read_counts) + 1, 25)
#                plt.yticks(ind[1:], ind[1:])
#            else:
#                pass
#            plt.xlim(0, end)
#            plt.xlabel('\nPosition in Precursor')
#            plt.ylabel('Reads Mapped\n')
#            plt.title('{}\n'.format(name))
#            ax.spines['right'].set_visible(False)  # Removes right axis
#            ax.spines['top'].set_visible(False)  # Removes top axis
#            ax.yaxis.set_ticks_position('none')  # Keeps vertical ticks hidden
#            ax.xaxis.set_ticks_position('bottom')
#            if max(read_counts) > 10000:
#                plt.savefig(map_folder + '/top/{0}_plant-precursors_bowtie1_a_'
#                                        'T_mm0_gap0split.svg'.format(name), 
#                            bbox_inches = 'tight', format = 'svg')
#            else:
#                plt.savefig(map_folder + '{0}_plant-precursors_bowtie1_a_T_'
#                                         'mm0_gap0split.svg'.format(name), 
#                            bbox_inches = 'tight', format = 'svg')
#            plt.close('all')
#            
            read_lines_set = set()
            for line in read_lines:
                read_lines_set.add('{0} l={1} a={2}'.format(line, 
                                          len(re.split(r'([AGCTU]+)', line)[1]),
                                                    read_lines.count(line)))
            read_lines_sorted = list(sorted(read_lines_set, key = lambda x: \
                                      (len(re.split(r'([AGCTU]+)', x)[0] + \
                                           re.split(r'([AGCTU]+)', x)[1]),
                                    2 ** -len(re.split(r'([AGCTU]+)', x)[1])), 
                                    reverse = False))

            for line in read_lines_sorted:
                f.write('{0}\n'.format(line))
            f.write("Mature: {0}\nStar: {1}\nMature Tail: {2}\nStar Tail: "
                    "{3}\nLoop: {4}\n".format(round(100*mat_count/i, 2), 
                                              round(100*star_count/i, 2), 
                                              round(100*mat_tail_count/i, 2), 
                                              round(100*star_tail_count/i, 2), 
                                              round(100*loop_count/i, 2)))
            loops_per = (loops*100)/(matures + stars + loops + mature_tails + star_tails)
            if loops_per >= 1:
                print("Had >1% loop reads: {0}\n{1}% Loops ({2})\tMatures {3}\tStars {4}"\
                      .format(name, round(loops_per, 2), loops, matures, stars))
            else:
                pass
            
    BTIRed_bam = pysam.AlignmentFile(docs_folder + 'BTIRed' + bam_ext)
    G002_bam = pysam.AlignmentFile(docs_folder + 'G002' + bam_ext)
    G006_bam = pysam.AlignmentFile(docs_folder + 'G006' + bam_ext)

    bams = [G006_bam, G002_bam, BTIRed_bam]  # close after one read through


total_reads = matures_summary + stars_summary + mature_tails_summary + \
              star_tails_summary + loops_summary + others
with open(map_folder + 'summary_reads.txt'.format(name), 'w') as f:
    f.write("Matures: {0}\nStars: {1}\nMature Tails: {2}\n""Star Tails: {3}\n"
            "Loops: {4}\nOthers: {5}\n\n".format(matures_summary, stars_summary, 
                                  mature_tails_summary, star_tails_summary, 
                                  loops_summary, others))
    f.write("Matures: {0}%\nStars: {1}%\nMature Tails: {2}%\n""Star Tails: {3}%\n"
            "Loops: {4}%\nOthers: {5}%\n\n".format(round(100*matures_summary/total_reads, 2), 
                                  round(100*stars_summary/total_reads, 2), 
                                  round(100*mature_tails_summary/total_reads, 2), 
                                  round(100*star_tails_summary/total_reads, 2), 
                                  round(100*loops_summary/total_reads, 2),
                                  round(100*others/total_reads, 2)))
    f.write("Mature: {0}%\nStar: {1}%\nMature Tail: {2}%\nStar Tail: {3}%\n"
            "Loop: {4}%\n".format(round(100*mat_count_sum/i_sum, 2), 
                                  round(100*star_count_sum/i_sum, 2), 
                                  round(100*mat_tail_count_sum/i_sum, 2), 
                                  round(100*star_tail_count_sum/i_sum, 2), 
                                  round(100*loop_count_sum/i_sum, 2)))

count_plot_list = [100*mature_tails_summary/total_reads, 
                   100*matures_summary/total_reads, 
                   100*loops_summary/total_reads, 
                   100*stars_summary/total_reads, 
                   100*star_tails_summary/total_reads,
                   100*others/total_reads]
nt_plot_list = [100*mat_tail_count_sum/i_sum, 
                100*mat_count_sum/i_sum, 
                100*loop_count_sum/i_sum, 
                100*star_count_sum/i_sum, 
                100*star_tail_count_sum/i_sum, 
                0]
plot_list = np.array(list(zip(count_plot_list, nt_plot_list)))
darkdarkgreen = '#015b00'
colors = ['black', darkdarkgreen, 'gold', 'firebrick', 'black', 'midnightblue']
label_list = ['Mature Tail', 'Mature miRNA', 'Loop', 'Star miRNA', 'Star Tail', 
              'Fragments']
plot_dict = {}

fig = plt.figure(figsize = (5.25, 3.25))
ax = plt.subplot(1, 1, 1)
left = 0
N = [0.05, 0.75]
for entry, color, label in zip(plot_list, colors, label_list):
    print("{0}:\t{1}% (/nt)\t{2}% (/read)".format(label, round(entry[0], 2), round(entry[1], 2)))
    plot_dict[label] = ax.barh(N, entry, height = 0.55, color = color, 
                               left = left, linewidth = 0, zorder = 5)
    left += entry
plt.legend([plot_dict['Mature miRNA'][0], plot_dict['Star miRNA'][0], 
            plot_dict['Loop'][0], plot_dict['Mature Tail'][0], 
            plot_dict['Star Tail'][0], plot_dict['Fragments'][0]], 
           ['Mature miRNA', 'Star miRNA', 'Loop', 'Mature Tail', 'Star Tail', 
            'Fragments'], frameon = False, loc = 'upper left', 
           bbox_to_anchor = (0.15, -0.175), ncol = 2)
plt.xlim(right = 100)
plt.yticks((0.05, 0.75), ('per Read', 'per Nucleotide'))
plt.xlabel('Percentage')
ax.spines['right'].set_visible(False)  # Removes right axis
ax.spines['left'].set_visible(False)  # Removes left axis
ax.yaxis.set_ticks_position('none')  # Keeps vertical ticks hidden
ax.xaxis.set_ticks_position('both')  # Shows x-axis ticks
ax.tick_params(axis = 'x', direction = 'in')  # Directs the x-axis ticks inward
plt.subplots_adjust(top = 0.75, wspace = 0.5)

plt.savefig(map_folder + 'summary_alignment.svg', bbox_inches = 'tight', 
            format = 'svg')
#plt.show()

G006_not_m_s.close()
G002_not_m_s.close()
BTIRed_not_m_s.close()