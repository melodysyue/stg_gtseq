# Sturgeon_gtseq
 

# Data

`./data/baseline/`: 
  - stg.G30.meanDP15.mac3.G70.l50.ID.vcf: vcf file after filtering (328 individuals, 55,112 SNPs) with ID column.
  - stg.G30.meanDP15.mac3.G70.l50.id: list of snp IDs in the filtered vcf file. 
  - popmap_sturgeon_N328.txt: popmap for 328 individuals after filtering.
  
`./data/panel_development/`
  - sumstats_acrosspops.csv: summary statistics for 55,112 SNPs, generated using polyfreqs output files. 
  - N502_primer.csv: primer information of 502 loci in the final optimized GTseq panel.
  - GTseq_sturgeon_primer_probe_final502.txt: primer/probe file for the final GTseq panel with 502 loci. 
  
`./data/gtseq/`: 
  - sturgeon_sampleinfo_updated.csv: sample meta information for 936 samples used in GTseq.
  - parentage.csv: raw parentage relationship information (unconfirmed).
  - true_parentage.txt: true parent-offspring trio relationship (confirmed).
  - minAF01_AvgReadDepth200_429loci_773SNPs_AlleleRatioCategory.list: a list of 773 informative SNPs with good coverage (average read depth >=200, AF>=1%) along with their allele ratio categories (great, okay, tbd).
  - `/scatterPlots_429loci_773SNPs/`: allele ratio plots by three categories (great, okay, tbd). 
  - polyGenResults_singSNP_letter_good.txt: GTscore/polyGen genotype output (4N).
  - polygene input files:
    - polygene_geno_pop_input.txt: input for 180 samples for population analysis (all 773 SNPs).
    - polygene_geno_pop_greatLoci.txt: input for 180 samples for population structure analysis (258 great SNPs only).
    - polygene_trio_geno_input.txt: input for 177 samples used for parentage analysis (all 773 SNPs).
  - polygene output files:
    - o_parentage_paternity.txt: output for likelihood-based paternity analysis;
    - o_relatedness_edited.txt: output for pairwise relatedness analysis.

# Scripts

The following scripts are for reference only. Please adjust accordingly for your computing environment and working directory.

`./scripts/`: 
 - AlleleCount_vcf_4n.py: python script to generate input matrices to run polyfreqs for a ploidy of 4N. Adopted from the original AlleleCount.py from Limborg, Larson, Seeb, & Seeb (2017).
 - polyfreqs_parallel.r & polyfreqs_parallel.slurm: run polyfreqs for all pops in parallel.
 - polyfreqs_postrun.rmd: process polyfreqs output files and generate summary statistics.
 - MI_calculation.r & MI_plot.r: analyses of Mendelian incompatibilities. 
 - paternity_LODplot.r: likelihood-based paternity plot.
 - relatedness.r: pairwise relatedness analysis.
 

