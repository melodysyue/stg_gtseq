########################################################################################################################################
########################################################################################################################################
###### Python 2 script to open vcf files from stacks output & create AlleleCount matrices for TOTAL reads and for REFERENCE reads ######
########################################################################################################################################
########################################################################################################################################

###################################
## Defining input & output files ##
###################################

#Read data:
import os

os.chdir("/scratch-lustre/wlarson/yue/sturgeon")
vcf_stacks_alleles = 'stg.G30.meanDP15.mac3.G70.l50.ID.vcf' 							#Define the name of the vcf input file with allele read counts for all loci

#Define and read in list with 4n loci names:
Dup_4n_loci_file = 'stg.G30.meanDP15.mac3.G70.l50.id' 											#Define file with list of all uplicated (4n) locus names
Dup_4n_loci_list = []              												#Read the file with list of 4n loci and create a list with 4n loci names:
with open (Dup_4n_loci_file, "r") as Dup_4n_loci_infile:
    for line in Dup_4n_loci_infile:
        locus = line.strip()
        Dup_4n_loci_list.append(locus)


## Define output files ##
Total_Reads_4n = 'Total_Reads_4n.txt' 											#Define output file with matrix of Total allelelic read counts for 4n loci
Reference_Reads_4n = 'Reference_Reads_4n.txt' 									#Define output file with matrix of Reference allelelic read counts for 4n loci
STACKS_genotypes_4n = 'STACKS_genotypes_4n.txt'									#Define output file with matrix of STACKS called genotypes for 4n loci


###############################################################################################
## For each locus (RADtag): Parse vcf file info and put allele read counts into dictionaries ##
###############################################################################################

RADtag_TotalAlleleCount_dict = {}
RADtag_ReferenceAlleleCount_dict = {}
RADtag_STACKS_Genotypes_dict = {}
RADtags = []

with open(vcf_stacks_alleles, 'r') as INFILE:
   #Looping over each line in the vcf input file
   for line in INFILE:
        linecontent = line.strip().split('\t')[0::] 												#Read line content and split by tab
     
        if len(linecontent) == 1:
            print 'Info line; skipping this one'            										#Check if line is just an info line; if yes skip this line

        else:

            if linecontent[0] == '#CHROM':                  										#Check if line is the header line - if yes; extract sample names
                individuals = line.strip().split('\t')[9::] 										#Record individual sample names (skipping first 9 info columns)
                print 'Header line; sample names are extracted for downstream analyses'
            
            else:
                #Read next line with genotypes and allele counts for each RADtag
                
                genotypes = line.strip().split('\t')[9::]
                RADtag = linecontent[2]     														#Extract RADtag name
                RADtags.append(RADtag)    															#Append RADtag name to list of all RADtag names
                print 'Processing RADtag: ' + RADtag
            
                List_of_Total_alleles_count = []        											#Define empty list to store Total allele counts for each individual
                List_of_Reference_alleles_count = []    											#Define empty list to store Reference allele counts for each individual

                List_of_STACKS_genotypes = []    													#Define empty list to store STACKS called genotypes


                for i in range(len(genotypes)):
                    #Loop over all individual allele counts for RADtag:
                    allele_counts = genotypes[i].strip().split(':')[2] 								#Parsing out the counts for the two alleles
                    
                    ###Loop over all STACKS genotypes for RADtag:
                    STACKS_genotypes = genotypes[i].strip().split(':')[0] 							#Parsing out the STACKS genotype
                    
                    #If number of reads of individual RADtag is zero (i.e. missing data) - record 0; else record allele counts:
                    if allele_counts.strip().split(',')[0] == '.':
                        Reference_allele_count = 0      											#Record missing data with a count = 0
                        Total_alleles_count = 0         											#Record missing data with a count = 0

                    else:
                        Reference_allele_count = int(allele_counts.strip().split(',')[0])         	#Recording observed number of the reference allele (to be printed in final matrix)
                        alternate_allele_count = int(allele_counts.strip().split(',')[1])         	#Recording observed number of the alternate allele
                        Total_alleles_count = Reference_allele_count + alternate_allele_count     	#Recording total number of observed alleles (to be printed in final matrix)

                    #Append individual allele counts to the RADtag list
                    List_of_Reference_alleles_count.append(Reference_allele_count)
                    List_of_Total_alleles_count.append(Total_alleles_count) 

                    ###Append individual STACKS genotypes to the respective list
                    List_of_STACKS_genotypes.append(STACKS_genotypes)
                        
                #Write Total and Reference allele counts for RADtag out to the respective dictionaries for Total and Reference allele counts
                RADtag_TotalAlleleCount_dict[RADtag] = List_of_Total_alleles_count
                RADtag_ReferenceAlleleCount_dict[RADtag] = List_of_Reference_alleles_count
                
                ###Write STACKS genotypes to the respective dictionary
                RADtag_STACKS_Genotypes_dict[RADtag] = List_of_STACKS_genotypes


#################################################################################################################################
## Create and write matrices with 'Total' and 'Reference' allele counts per individual for use in e.g. the POLYFREQS software ###
#################################################################################################################################

## Create new lists for 'mapped' 2n and 4n loci represented in the .vcf data file ##
RADtags_4n = []
for RADtag in RADtags:
    if RADtag in Dup_4n_loci_list:
        RADtags_4n.append(RADtag)

    else:
        print 'RADtag: ' + RADtag + ' is not polymorphic in the population data set (vcf file)'


## TOTAL ALLELE COUNTS: Write matrices for 4n loci respectively ##
#4n loci:
with open(Total_Reads_4n, 'w')as OUTFILE_4n:
    OUTFILE_4n.write('\t'.join(RADtags_4n))   					#Write RADtag names in the header (first) line
    OUTFILE_4n.write('\n')

    for i in range(len(individuals)):          					#Loop over each individual to parse out allele counts
        ind_Total_Reads_4n = []
        for RADtag in RADtags_4n:
            ind_Total_Reads_4n.append(str(RADtag_TotalAlleleCount_dict[RADtag][i]))

        OUTFILE_4n.write(individuals[i] + '\t')            		#Write name of individual in row 1
        OUTFILE_4n.write('\t'.join(ind_Total_Reads_4n))    		#Write out all Total allele counts for individual i
        OUTFILE_4n.write('\n')


## REFERENCE ALLELE COUNTS: Write matrices for 4n loci##
#4n loci:
with open(Reference_Reads_4n, 'w')as OUTFILE_4n:
    OUTFILE_4n.write('\t'.join(RADtags_4n))    					#Write RADtag names in the header (first) line
    OUTFILE_4n.write('\n')

    for i in range(len(individuals)):          					#Loop over each individual to parse out allele counts
        ind_Reference_Reads_4n = []
        for RADtag in RADtags_4n:
            ind_Reference_Reads_4n.append(str(RADtag_ReferenceAlleleCount_dict[RADtag][i]))

        OUTFILE_4n.write(individuals[i] + '\t')            		#Write name of individual in row 1
        OUTFILE_4n.write('\t'.join(ind_Reference_Reads_4n))    	#Write out all Total allele counts for individual i
        OUTFILE_4n.write('\n')



## STACKS GENOTYPES: Write matrices for 4n loci ##
#4n loci:
with open(STACKS_genotypes_4n, 'w')as OUTFILE_4n:
    OUTFILE_4n.write('\t'.join(RADtags_4n))    					#Write RADtag names in the header (first) line
    OUTFILE_4n.write('\n')

    for i in range(len(individuals)):          					#Loop over each individual to parse out allele counts
        ind_STACKS_genotypes_4n = []
        for RADtag in RADtags_4n:
            ind_STACKS_genotypes_4n.append(str(RADtag_STACKS_Genotypes_dict[RADtag][i]))

        OUTFILE_4n.write(individuals[i] + '\t')            		#Write name of individual in row 1
        OUTFILE_4n.write('\t'.join(ind_STACKS_genotypes_4n))   	#Write out all Total allele counts for individual i
        OUTFILE_4n.write('\n')


#######################################################################
############################### THE END ###############################
#######################################################################
