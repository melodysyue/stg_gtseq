library(magrittr) #to use %>% since Turing doesn't have tidyverse installed;
library(parallel) #to use mclapply and mcmapply
library(polyfreqs) 

##check how many cores do you have
##detectCores()
##on Turing, usually it has 16-32 cores. Let's use 16 just to be safe.And, we happen to have 16 populations to process. So it is perfect;
  
##read the read count files and convert them to matrix as required by polyfreqs; 
tot_files <- list.files(pattern = "*_Total_Reads_4n.txt") %>% 
  mclapply(read.table,row.names=1, header=TRUE, mc.cores=16) %>% 
  mclapply(as.matrix, mc.cores=16)

ref_files <- list.files(pattern = "*_Reference_Reads_4n.txt") %>% 
  mclapply(read.table,row.names=1, header=TRUE, mc.cores=16) %>% 
  mclapply(as.matrix, mc.cores=16)


##generate the list of output file names using pop names 
pops.name <- list.files(pattern = "*_Total_Reads_4n.txt")
mcmc_files <- sub("_Total_Reads_4n.txt", "_mcmc.out", pops.name)
het_obs_files <- sub("_Total_Reads_4n.txt", "_het_obs.txt", pops.name)
het_exp_files <- sub("_Total_Reads_4n.txt", "_het_exp.txt", pops.name)


##create a function with which you can modify the input  with the polyfreqs function;
##note: each input is a list, all the inputs should have the same length. Otherwise the shorter one(s) will get recycled;  
polyfreqs_par <- function(total,reference,output,obs,exp)
{
  p <- polyfreqs(total, reference, ploidy=4, 
                     iter=100000,
                     thin=100,
                     outfile=output)
  write.table(p$het_obs, obs, quote=F, row.names=F, col.names=F)
  write.table(p$het_exp, exp, quote=F, row.names=F, col.names=F)
}


##run polyfreqs in parallel (multi-core mode)
mcmapply(polyfreqs_par, tot_files, ref_files, mcmc_files, het_obs_files, het_exp_files, mc.cores=16)

