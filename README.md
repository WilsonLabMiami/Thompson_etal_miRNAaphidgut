# The green peach aphid gut contains host plant microRNAs identified by comprehensive annotation of Brassica oleracea small RNA data
###### Max C. Thompson, Honglin Feng, Stefan Wuchty, and Alexandra C. C. Wilson

- Code included here was used to format, analyze, and display the data and results present in the above paper.
- Paths within the Python scripts have been mostly altered to input('') phrases that explain what file path needs to go in there. The extensions shown tell what type of file is currently expected.
  - The scripts are not expected to be run while manually inputting these file paths each time (though it would be possible, just painful). Instead, swap them out for strings that have the file paths you will be using in the scripts.


## Initial Workflow

1. Genomes for [Myzus persicae](https://bipaa.genouest.org/sp/myzus_persicae/download/genome/CloneG006_v2/), [Buchnera aphidicola](https://www.ncbi.nlm.nih.gov/assembly/?term=(%22Buchnera+aphidicola+(Myzus+persicae)%22)+AND+(G006%5BInfraspecific+Name%5D+OR+G002%5BInfraspecific+Name%5D+OR+USDA%5BInfraspecific+Name%5D)), and [Brassica oleracea var. oleracea](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/695/525/GCF_000695525.1_BOL/) were manually downloaded from AphidBase and NCBI. 
2. Index the genomes using bowtie2
3. Filter small-RNA-seq data to your specifications
4. Align small-RNA-seq data to the respective genomic libraries using bowtie2
5. Analyze the alignments with respect to your genomes of interest using [`percentage_calculations.py`](percentage_calculations.py)
   - Parses bowtie2 output SAM files and assigns each read alignment to a genome set.
   - Genome sets are compared for overlap and percentages are assigned to each section based on the total reads in the fasta file used in making these alignments.
6. Generate read subset files for futher analysis using [`elimMatched.py`](elimMatched.py) and ['elimUnmatched_plant-only.py'](elimUnmatched_plant-only.py)
   - elimMatched will parse a given SAM file and associate fasta file to generate a new fasta file with only the reads that have not met the criteria for a match
     - This will set up the "unknown" reads to be aligned against the bacterial and viral libraries
   - elimUnmatched_plant-only works int he same manner as above only to output reads with have met the criteria for a match to the plant genome


## Bacterial and Viral Assignment of Unknown Reads
1. Use [`fetchGenomes.py`](fetchGenomes.py) to collect genomes for bacterial and viral libraries
   - Pulls genomes from literature review into downloaded FASTA files.
2. Index the genomes using bowtie2
3. Align small-RNA-seq data to the respective genomic libraries using bowtie2
4. Analyze the alignments with respect to your genomes of interest using [`percentage_calculations.py`](percentage_calculations.py)
5. Produce panels for figures 2 and 3 using [`figure_2-3.py`](figure_2-3.py)
   - This will generate SVGs that can be manipulated in vector graphics software (ie: Illustrator).


## Plant miRNA Analysis
1. Collect SRA datasets for plant small-RNA-seq
2. Filter small-RNA-seq data to your specifications
3. Align small-RNA-seq data to the respective genomic libraries using bowtie2
4. Annotate potential miRNA precursors using [ShortStack](https://github.com/MikeAxtell/ShortStack)
5. Assess homology to existing precursors using local BLAST against [miRBase](http://mirbase.org/ftp.shtml)
6. Collect accepted miRNA precursors and assocaited information into a csv with similar formating to Supplementary Table 3.
7. Use [`mirFinder.py`](mirFinder.py) to search small RNA reads for exact copies of 5' and 3' sequences from annotated plant miRNA precursors
   - mirFinder can also be fed csvs of additional miRNA sequences to see if they are represented in your found miRNAs
###### Target Analysis 
8. Predict targets in both the plant genome and the aphid genome using multiple miRNA target prediction software packages
   - [miRanda](http://www.microrna.org/microrna/getDownloads.do), [PITA](https://genie.weizmann.ac.il/pubs/mir07/mir07_exe.html), and [RNAhybrid](https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid) were used for aphid predictions
   - [psRNAtarget](http://plantgrn.noble.org/psRNATarget/?dowhat=Help) and [TargetFinder](https://github.com/carringtonlab/TargetFinder) were used for plant predictions
   - Compact output should be used and will need to be parsed to match input formating for miRNA_intersect
9. Parse miRNA targets to find target sites where all programs agree on the seed region using [`miRNA_intersect.py`](miRNA_intersect.py) for aphid data and [`miRNA_intersect_plant.py`](miRNA_intersect_plant.py) for plant data
   - Input files may require some modification to match parameters within or parameters can be adjusted
   - Output sites are agreed upon by all packages
10. Collect protein sequences for the genes found to be targeted and assess their function using HMMs via [EggNOG-mapper](http://eggnog-mapper.embl.de/).
    - Repeat this process for the genome of interest as it will be used to generate null distributions.
11. Process target COG assignments in relation to genome COG assignments and produce SVGs for Figures 5 and 6 using [`isoformAnalysis_small.py`](isoformAnalysis_small.py) and [`isoformAnalysis_small_plant.py`](isoformAnalysis_small_plant.py) for the aphid and plant genomes respectively
    - Simulations for the generation of a null distribution can be further parallelized across multiple jobs on a cluster through depositing the probability dictionaries into pickle files each time and later amassing a main dictionary to save.
      - It is recommened that the simulations only be run once and then loaded from a pickle dictionary each time as this is the most time and computatinoally intensive part.
12. Visualize structures for the precursors using Python subprocess calls to VARNA CLI with [`varna_structure_from_fass.py`](varna_structure_from_fass.py)
    - [VARNA](http://varna.lri.fr/index.php?lang=en&page=downloads&css=varna) will need to be downloaded and accessible in the commmand line environment this Python script is run from.
    - A FASS (fasta and secondary structure) file will also be needed with all of the precursors structures. This can come from ShortStack, but the rounds of optimizing and trimming within may leave the brackets not matching up, and VARNA will be unable to use it. Instead [RNAfold](https://www.tbi.univie.ac.at/RNA/#download) can be used to produce this file.
###### Precursor Processing Analysis
13. Align all processed small-RNA-seq data against plant miRNA precursors in step 6 above using bowtie2
14. Filter alignments to include only those with zero mismatches and zero gaps in a BAM file
    - Filtering and file conversion can be accomplish with Bamtools and Samtools respectively
    - Filtered BAM files will also need to be indexed
15. Process precursor alignments using [`bamAnalysis.py`](bamAnalysis.py) to generate text file visualizations of small read alignment for each precursor
    - An example txt file is available for [bol-miR168b-3p (Cluster_2535)](examples/bol-miR168b-3p_Cluster_2535_reads.txt)
    - bamAnalysis will also calculate the percentage of reads assigned to mature, star, loop, 5', and 3' sections and produce and SVG figure that displays this in a bar chart
    - You can also create histogram figures for each precursor based on per nucleotide mapping
    - This uses Python packages only built for Mac/UNIX systems. Linux Subsystem for Windows can be used to run these on a Windows machine
16. Using the text file outputs from bamAnalysis, generate read alignment figures for each precursor using [`read_mapping_figure.py`](read_mapping_figure.py)
    - This uses Python packages only built for Mac/UNIX systems. Linux Subsystem for Windows can be used to run these on a Windows machine

