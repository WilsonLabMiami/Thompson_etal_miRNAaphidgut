# import matplotlib as mpl
# mpl.use('SVG')


import time
start_time = time.time()

from sys import argv
import csv
from numpy.random import choice
from joblib import Parallel, delayed 
import multiprocessing
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import seaborn as sns
import pickle


cogs_dict_g = {}  # holds genome COGs; turns to percentages
cogs_dict_utr_g = {}

cogs_dict = {}  # holds targets COGs; turns to percentages
cogs_dict_utr = {}
cogs_per = []  # used for plotting target composition percentages in order

cogs_dict_final = {}  # holds target - genome percentages with labels
labels_final = []  # holds ordered labels
cogs_per_final = []   # used for plotting target - genome composition percentages in order

probs = []  # holds probability (weights) for each COG in genome
probs_dict = {}  # probability with label
probs_dict_bc = {}  # for boxcox transformation
pvals_raw = []
pvals_dict = {}
ci = {}  # holds confidence intervals


filelist = ['bol4_myseq0.fa.emapper.annotations', 
            'bol4_myseq16000.fa.emapper.annotations', 
            'bol4_myseq24000.fa.emapper.annotations', 
            'bol4_myseq32000.fa.emapper.annotations', 
            'bol4_myseq8000.fa.emapper.annotations', 
            'bol4_myseq12000.fa.emapper.annotations', 
            'bol4_myseq20000.fa.emapper.annotations', 
            'bol4_myseq28000.fa.emapper.annotations', 
            'bol4_myseq4000.fa.emapper.annotations']
            # genome COG assignments

for file in filelist:
    with open(file) as f:
        csvreader = csv.reader(f, delimiter = '\t')
        for row in csvreader:
            cog = row[11]
            name = row[0]

            if cog == "S":
                pass
            else:
                if len(cog) > 1:
                    cogs_list_g = cog.split(', ')  # because entries can match to multiple COGs
                    for entry in cogs_list_g:
                        # entry is cog
                        cogs_dict_g.setdefault(entry, []).append(name)
                else:
                    cogs_dict_g.setdefault(cog, []).append(name)
total_g = 0
for key in cogs_dict_g:
    total_g = total_g + len(cogs_dict_g[key])

with open('oleracea_transcripts_pep_exact_mirs_targets.fa.emapper.annotations') as f:
    # target COG assignments
    csvreader = csv.reader(f, delimiter = '\t')
    header = next(csvreader)
    for row in csvreader:
        cog = row[11]
        name = row[0]
        if cog == 'S':
            pass
        else:
            if len(cog) > 1:
                cogs_list = cog.split(', ')
                for entry in cogs_list:
                    cogs_dict.setdefault(entry, []).append(name)
            elif len(cog) == 1:
                cogs_dict.setdefault(cog, []).append(name)
            else:
                pass
total_t = 0
for key in cogs_dict:
    total_t = total_t + len(cogs_dict[key])

labels = ["RNA processing and modification", 
          "Chromatin structure and dynamics", 
          "Energy production and conversion", 
          "Cell cycle control, cell division,\nchromosome partitioning", 
          "Amino acid transport and metabolism", 
          "Nucleotide transport and metabolism", 
          "Carbohydrate transport and metabolism", 
          "Coenzyme transport and metabolism", 
          "Lipid transport and metabolism", 
          "Translation, ribosomal structure,\nand biogenesis", 
          "Transcription", 
          "Replication, recombination, and repair", 
          "Cell wall/membrane/envelope biogenesis",  
          "Posttranslational modification,\nprotein turnover, chaperones", 
          "Inorganic ion transport\nand metabolism", 
          "Secondary metabolites biosynthesis,\ntransport, and catabolism", 
          "Signal transduction mechanisms", 
          "Intracellular trafficking, secretion,\nand vesicular transport", 
          "Defense mechanisms", 
          "Extracellular structures", 
          "Nuclear structure", 
          "Cytoskeleton"]
#          COG long form labels over letters
#          "Cell motility",

labels_letters = ["A", "B", "C", "D", "E", "F", "G", "H", 
                  "I", "J", "K", "L", "M", "O", "P", 
                  "Q", "T", "U", "V", "W", "Y", "Z"]
                  # COG short form letters in the same order as labels
#                  "N",

for (long_form, short_form) in zip(labels, labels_letters):
    try:
        cogs_dict[long_form] = cogs_dict.pop(short_form)  
        # exchanges letter COG for label
        cogs_dict[long_form] = len(cogs_dict[long_form])/(total_t/100)
        # creates percentage from number of values (proteins) over the total
    except KeyError:
        pass

not_in_targets = ['Extracellular structures', 'Nuclear structure']
# these entries show up in the genome proteins but not targets
# this sorts it out for the figure production
for entry in not_in_targets:
    cogs_dict[entry] = 0


for (long_form, short_form) in zip(labels, labels_letters):
    try:
        cogs_dict_g[long_form] = cogs_dict_g.pop(short_form)
        cogs_dict_g[long_form] = len(cogs_dict_g[long_form])/(total_g/100)
        cogs_dict_final[long_form] = cogs_dict[long_form] - cogs_dict_g[long_form]
    except KeyError:
        pass


#%%
        
mid_time = time.time()
print(round(mid_time - start_time, 2))


sorted_dict_list = sorted(cogs_dict.items(), key=lambda x: (abs(x[1]), abs(cogs_dict_final[x[0]])), reverse = True)
# sorts the dictionary into a list to keep a consistent order
# second sorting key goes by magnitude of percent difference
# secondary sorting also assures the same order each time
for entry in sorted_dict_list:
    label = entry[0]
    percent = entry[1]
    labels_final.append(label)
    probs_dict[label] = []
    cogs_per.append(percent)
    cogs_per_final.append(cogs_dict_final[label])
    # this loop estabilishes many lists for the figure at once, keeping the order the same

for label in labels_final:
    probs.append(cogs_dict_g[label]/sum(cogs_dict_g.values()))
    # total percent is >100 but this needs a ratio, so it's divided by the total percent
    # this generates probabilities for the MC simulation


def choices_list(i):
    choices = []
    counted_dict = {}
    for i in range(0, 576):
        name = list(choice(labels_final, 1, p = probs))[0]
        choices.append(name)
    for label in labels:
        percent = (choices.count(label))/(len(choices)/100)
        counted_dict.setdefault(label, []).append(percent)
    # runs simulations to generate null distributions
    return counted_dict

#num_cores = multiprocessing.cpu_count()
num_cores = 4
counted_dicts = Parallel(n_jobs=num_cores)(delayed(choices_list)(i) for i in range(1, 100001))

for counted_dict in counted_dicts:
    for label in counted_dict.keys():
        probs_dict.setdefault(label, []).extend(counted_dict[label])

with open('probs_dict_exact_mirs_plant_100k_{0}.pkl'.format(argv[1]), 'wb') as dict_file:
    pickle.dump(probs_dict, dict_file, protocol = pickle.HIGHEST_PROTOCOL)
    # argv and this pickle formating is used to allow parallelization across
    # several jobs if submitting to a cluster. Fewer simulations can be run 
    # each time and ammased later.
    # Once you have enough simulations, build a final dictionary and save it
    # with pickle to avoid needing to resimulate

mid_time_2 = time.time()
print(round(mid_time_2 - mid_time, 2))

#%%

for label in labels_final:
    target_per = cogs_dict[label]
    total_sims = len(probs_dict[label])
    if cogs_dict_final[label] < 0:
        lesser_pval = len([i for i in probs_dict[label] if i < target_per])/total_sims
        lesser_ci = sorted(probs_dict[label])[int(0.025*total_sims)]
        ci.setdefault(label, []).append(lesser_ci)
        pvals_raw.append(lesser_pval)
        if lesser_pval < 0.05 and lesser_pval != 0:
            print('{0} (Lesser): {1}\nConfidence interval: {2}\n'.format(label, lesser_pval, lesser_ci))
            pvals_dict[label] = {'lesser', lesser_pval}
        else:
            pass
    elif cogs_dict_final[label] > 0:
        greater_pval = len([i for i in probs_dict[label] if i > target_per])/total_sims
        greater_ci = sorted(probs_dict[label])[int(0.975*total_sims)]
        ci.setdefault(label, []).append(greater_ci)
        pvals_raw.append(greater_pval)
        if greater_pval < 0.05 and greater_pval != 0:
            print('{0} (Greater): {1}\nConfidence interval: {2}\n'.format(label, greater_pval, greater_ci))
            pvals_dict[label] = {'greater', greater_pval}
        else:
            pass
    else:
        pass
    # Note that with empirical p-value determination, there must be a simulation
    # as extreme as, or more extreme than, the target value; otherwise, the
    # p-value must be represented as p < 1/(simulation number)

with open('ci_plant.pkl', 'wb') as dict_file:
    pickle.dump(ci, dict_file, protocol = pickle.HIGHEST_PROTOCOL)
with open('pvals_dict_plant.pkl', 'wb') as dict_file:
    pickle.dump(pvals_dict, dict_file, protocol = pickle.HIGHEST_PROTOCOL)

#%%

sns.set_style('white')  # sets sns background styling to white
vmax = max(cogs_per_final)
vmin = min(cogs_per_final)
mid = 1 - vmax / (vmax + abs(vmin))

top = cm.get_cmap('Blues_r', 256000)
middle_blue = cm.get_cmap('Blues_r', 256000)
middle_red = cm.get_cmap('Reds', 256000)
bottom = cm.get_cmap('Reds', 256000)

newcolors = np.vstack((top(np.linspace(0.1, 0.85, 251000)),
                       middle_blue(np.linspace(0.85, 0.9, 5000)),
                       middle_red(np.linspace(0.05, 0.15, 5000)),
                       bottom(np.linspace(0.15, 0.9, 251000))))
centered_cm = ListedColormap(newcolors, name='Blue_Red')

colors_2 = []
for entry in cogs_per_final:
    colors_2.append(centered_cm((entry - -14.5)/(float(14.5) - -14.5)))

plot = plt.scatter([-14.5, 14.5], [-14.5, 14.5], c = [-14.5, 14.5], 
                   cmap = centered_cm)
# establishes colorbar in scatterplot
plt.close()  # clears scatterplot
 

#####################
# Figure Production #
#####################

fig = plt.figure(figsize = (15, 15))
ind = np.arange(len(labels_final))
ax = plt.subplot(1, 22, (1, 7))
barlist = plt.barh(ind, cogs_per, height = 1, align = 'edge', zorder = 10)
for i, color in zip(range(0, len(barlist)), colors_2):
    barlist[i].set_color('#808080')
    if i%2 == 0:
        ax.axhspan(i, 1 + i, color = 'white', alpha = 0.15, zorder = 1)
    else:
        ax.axhspan(i, 1 + i, color = 'dimgrey', alpha = 0.15, zorder = 1)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
plt.yticks(ind + 0.5, labels_final, rotation = "horizontal", fontsize = 12, 
           multialignment = 'center', zorder = 12)
plt.ylim([0, ind.size])  # Removes whitespace to right side
plt.xlim([30, 0])
plt.xlabel('Percent Composition of Targets', x = 0.5, fontsize = 12)
plt.xticks((0, 5, 10, 15, 20, 25, 30), (0, 5, 10, 15, 20, 25, 30), fontsize = 12)
ax.tick_params(direction = 'out')
ax.spines['right'].set_visible(False)  # Removes right axis
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)  # Removes top axis
ax.yaxis.set_ticks_position('none')  # Keeps vertical ticks hidden
ax.xaxis.set_ticks_position('bottom')  # Keeps horizontal ticks hidden on top


font = {'family': 'serif',
        'color':  'darkgrey',
        'weight': 'bold',
        'size': 24,
        }

ind = np.arange(len(labels_final))
ax = plt.subplot(1, 22, (8, 20))
barlist = plt.barh(ind, cogs_per_final, height = 1, align = 'edge', zorder = 10)
for i, color, label in zip(range(0, len(barlist)), colors_2, labels_final): # entry and ci missing now
    barlist[i].set_color(color)
    if i%2 == 0:
        ax.axhspan(i, 1 + i, color = 'white', alpha = 0.15, zorder = 1)
    else:
        ax.axhspan(i, 1 + i, color = 'dimgrey', alpha = 0.15, zorder = 1)
    if label in pvals_dict.keys():
        ax.axvline(ci[label][0] - cogs_dict_g[label], color = 'darkgrey', 
              ymin = (i)/22, ymax = (1 + i)/22, linewidth = 1.5, zorder = 11)
        plt.text(-13.5, i + 0.15, '*', fontdict = font, zorder = 11)
    else:
        ax.axvline(ci[label][0] - cogs_dict_g[label], color = 'darkgrey', 
              ymin = (i)/22, ymax = (1 + i)/22, linewidth = 1.5, zorder = 11)
plt.ylim([0, ind.size])  # Removes whitespace to right side
plt.xlim([-41, 15])
ax.axvspan(-15, -15.25, color = 'white', zorder = 11)
ax.axvspan(-40.7, -40.95, color = 'white', zorder = 11)
plt.yticks(ind, '')
plt.xlabel('Percent Enrichment Over Genome', fontsize = 12, x = (41/(41+15)))
plt.xticks((-15, -10, -5, 0, 5, 10, 15), (-15, -10, -5, 0, 5, 10, 15), fontsize = 12)
ax.tick_params(direction = 'out')
ax.spines['right'].set_visible(False)  # Removes right axis
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)  # Removes top axis
ax.yaxis.set_ticks_position('none')  # Keeps vertical ticks hidden
ax.xaxis.set_ticks_position('bottom')  # Keeps horizontal ticks hidden on top

ax = plt.subplot(1, 22, (21, 22))
plt.yticks([], '')
plt.xticks([], '')
ax.spines['right'].set_visible(False)  # Removes right axis
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)  # Removes top axis
ax.spines['bottom'].set_visible(False)
ax.yaxis.set_ticks_position('none')  # Keeps vertical ticks hidden
ax.xaxis.set_ticks_position('none')  # Keeps horizontal ticks hidden on top
plt.colorbar(plot, extend = 'both', extendfrac = 0.025)  # plots colorbar from earlier scatterplot


plt.subplots_adjust(wspace = 0.00)
if len(probs_dict[label]) == 1000000:
    plt.savefig('Figure_COG_Center_Legend_plant_1Mil_flipped.svg', 
                bbox_inches = 'tight', format = 'svg', dpi = 300)
else:
    plt.savefig('Figure_COG_Center_Legend_plant_flipped.svg', 
                bbox_inches = 'tight', format = 'svg', dpi = 300)
#plt.show()
plt.cla()
plt.clf()
plt.close("all")

end_time = time.time()
print(round(end_time - start_time, 2))
