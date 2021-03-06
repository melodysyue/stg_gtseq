######----------------------------------------------------------------------------------------------------------------
##calculate Mendelian incompatibilities ("percent mismatch") by locus and by trio

library(tidyverse)
library(gtools)
library(compare)

rm(list=ls())

###load genotype files
geno <- read.table("./data/gtseq/polyGenResults_singleSNP_letter_good.txt", check.names = F, header=T)
dim(geno)
ploidy=4

#check.names=F to avoid leading X and .
geno[1:5,1:5]
#remove , from the genotype
geno <- geno %>% 
  mutate_all(funs(str_replace_all(., ",","")))
geno[1:5,1:5]

#load relationship info
meta <- read.csv("./data/gtseq/sturgeon_sampeinfo_updated.csv", header=T)
meta$Sample <- str_replace_all(meta$Sample, "_", "-")

dam <- c("MM12-P014", "MM12-P029", "MM12-P041", "MM12-P117", "MM12-P173", "MM12-P687",
         "18BLA126", "18BLA150", "18BLA163","19BLA134", "19BLA148")
sire <- c("MM12-P016", "MM12-P027", "MM12-P035", "MM12-P126", "MM12-P174", "MM12-P178", 
          "18BLA111", "18BLA112","18BLA151", "18BLA152","18BLA166","18BLA168", 
          "19BLA115", "19BLA122","19BLA146","19BLA147")

offspring <- meta %>% 
  filter(Parentage.or.Populations.=="Parentage") %>% 
  select(Sample, Female_parent=Female.Parent, Male_parent=Male.Parent) %>% 
  filter(Female_parent != "Parent" & Female_parent !="Parents") %>% 
  pull(Sample) %>% 
  as.character()

#how many trios were retained?

good.ind <- colnames(geno)

dam.good <- good.ind[good.ind %in% dam] #4
sire.good <- good.ind[good.ind %in% sire] #11
offspring.good <- good.ind[good.ind %in% offspring] #162


###Calculate mendelian incompatibility
##we have to translate the genotypes so they are all in alphab. order, regardless of ref/alt allele
##we do this because we predict mismatches based on 
#potential genotypes from each parent pairs' gamete combinations (whether offspring is included)
##pairings of gametes could generate all different orders of genotypes; 
#thus we translate them to a standard format (alphab.)

#transform, samples as rows, snps as columns
geno <- t(geno) 
dim(geno)
snps <- colnames(geno) 

#replace missing genotype 0 with 0000
sum(geno[,snps]=="0")

for(r in 1:nrow(geno)){
  geno[r,snps] = gsub("^0$","0000",geno[r,snps])
}
sum(geno[,snps]=="0")#should now be zero

##check classes of heterozygote genotypes and see if they are in alphab order, before translation (substitution)
genotypes<-list()
genotypes_all<-NULL
for (l in 1:length(snps)){
  df<-as.data.frame(table(geno[,snps[l]]))
  genotypes[[l]]<-as.character(df$Var1)
  genotypes_all<-c(genotypes_all,geno[,snps[l]])
}
#should now be many fewer classes than before; with the alleles in alphab. order
table(genotypes_all)

###If genotypes are not in alphab order, do the substituion as following
#check how many genotype classes are currently included
nucleotides<-rbind(c("A","C"),c("A","G"),c("A","T"),c("C","G"),c("C","T"),c("G","T"))
het.translation<-as.data.frame(matrix(nrow=0,ncol=2))
p=4
gamA<-round(p/2-0.001,0)
gamB<-round(p/2+0.001,0)
for (n in 1:nrow(nucleotides)){
    AA<-as.data.frame(permutations(n=2,r=gamA,v=nucleotides[n,],repeats.allowed=TRUE))
    two.nuc<-unlist(lapply(as.list(data.frame(t(AA))), function(x) paste(x, collapse="")))
    AB<-as.data.frame(permutations(n=2,r=gamB,v=nucleotides[n,],repeats.allowed=TRUE))
    three.nuc<-unlist(lapply(as.list(data.frame(t(AB))), function(x) paste(x, collapse="")))
    combo.nuc<-unique(c(two.nuc,three.nuc))
    B<-as.data.frame(permutations(n=length(combo.nuc),r=2,v=combo.nuc,repeats.allowed=TRUE))
    paste.nuc<-unique(unlist(lapply(as.list(data.frame(t(B))), function(x) paste(x, collapse=""))))
    keep<-as.character(lapply(as.list(paste.nuc), function(x) length(unlist(strsplit(as.character(x),"")))))==p
    paste.nuc<-paste.nuc[keep]
    paste.nuc<-paste.nuc[!paste.nuc %in% c(paste(rep(nucleotides[n,1],times=p),collapse=""),paste(rep(nucleotides[n,2],times=p),collapse=""))]
    het.translation.sub<-as.data.frame(matrix(nrow=0,ncol=2))
    for (t in 1:length(paste.nuc)){
      df<-as.data.frame(table(strsplit(as.character(paste.nuc[t]),"")))
      het.translation.sub<-rbind(het.translation.sub,cbind(paste.nuc[t],paste( c( as.character(rep(df$Var1[1],times=df$Freq[1])) , as.character(rep(df$Var1[2],times=df$Freq[2])) ),collapse="" ) ))
      
    }
    het.translation.sub<-het.translation.sub[order(het.translation.sub$V2),]
    het.translation<-rbind(het.translation,het.translation.sub)
    print(nrow(het.translation))
  }
colnames(het.translation)<-c("orig","named")
head(het.translation)

##replace alternative heterozygotes with alphabetical heterozygotes (e.g. AAAT for TAAA or AATA)
#het.translation[,c(1,2)] = apply(het.translation[,c(1,2)], 2, function(x) as.character(x));
#geno[,snps] = apply(geno[,snps], 2, function(x) as.character(x));

#for(g in 1:nrow(het.translation)){
#  print(paste(c("finding ",as.character(het.translation[g,"orig"])," and replacing with ",as.character(het.translation[g,"named"])),collapse=""))
#  geno[geno==as.character(het.translation[g,"orig"])]<-as.character(het.translation[g,"named"])
#}


##check classes of heterozygote genotypes are alphabetical
#genotypes<-list()
#genotypes_all<-NULL
#for (l in 1:length(snps)){
#  df<-as.data.frame(table(geno[,snps[l]]))
#  genotypes[[l]]<-as.character(df$Var1)
#  genotypes_all<-c(genotypes_all,geno[,snps[l]])
#}
#should now be many fewer classes than before; with the alleles in alphab. order
#table(genotypes_all)

#add lifestage information
dam.geno <- geno[dam.good,]
sire.geno <- geno[sire.good,]
offspring.geno <- geno[offspring.good,]

dim(dam.geno)
dim(sire.geno)
dim(offspring.geno)

trio.geno <- rbind(sire.geno, dam.geno, offspring.geno)
dim(trio.geno)

write.table(trio.geno, "trio.geno.txt", quote=F, row.names = T)

info<-as.data.frame(cbind(rownames(trio.geno),
                          c(rep("adult",times=(length(sire.good)+length(dam.good))),rep("offspring",times=length(offspring.good))),
                          c(rep("male",times=length(sire.good)),
                            rep("female",times=length(dam.good)),
                            rep("offspring",times=length(offspring.good)))))
colnames(info)<-c("Sample","lifestage","sex")
table(info$lifestage)
table(info$sex)


##compare parental combo genotypes with each offspring, make ploidy screen first step, make gametes half of ploidy
#double reduction not currently considered
missing.genos="0000"
locus_mismatches<-as.data.frame(matrix(nrow=0,ncol=(length(snps)+3)))
summed_mismatches<-as.data.frame(matrix(nrow=0,ncol=6))
for (m in 1:nrow(dam.geno)){
    malleles<-strsplit(as.character(t(dam.geno[m,snps])),"")#returns list with vector for each locus
    malleles.combos<-lapply(malleles, function(x) combn(x,2,simplify=FALSE))
    mom.gametes<-lapply(lapply(malleles.combos, 
                               function(x) lapply(x, function(x) paste(x, collapse=""))), 
                        function(x) unique(unlist(x)))

    for (p in 1:nrow(sire.geno)){
      palleles<-strsplit(as.character(t(sire.geno[p,snps])),"")#returns list with vector for each locus
      palleles.combos<-lapply(palleles, function(x) combn(x,2,simplify=FALSE))
      pop.gametes<-lapply(lapply(palleles.combos, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      
      pot.offspr.genos<-list()
      
      for (l in 1:length(pop.gametes)) {
        raw.pot.offspr.genos<-as.character(unlist(as.list(unite(expand.grid(mom.gametes[[l]],pop.gametes[[l]]),
                                                                "geno",sep="", remove=TRUE))))
        raw.pot.offspr.genos[raw.pot.offspr.genos %in% het.translation[het.translation$orig %in% raw.pot.offspr.genos,"orig"]] <- as.character(het.translation[het.translation$orig %in% raw.pot.offspr.genos,"named"][match(raw.pot.offspr.genos, het.translation[het.translation$orig %in% raw.pot.offspr.genos,"orig"], nomatch = 0)])
        par.alleles<-unique(c(unique(unlist(strsplit(mom.gametes[[l]],""))),unique(unlist(strsplit(pop.gametes[[l]],"")))))
        if(sum(par.alleles %in% "0")>0) {
          raw.pot.offspr.genos<-missing.genos
        }
        pot.offspr.genos[[l]]<-raw.pot.offspr.genos
      }
        names(pot.offspr.genos)<-snps
        pot.offspr.genos ###this all possible genotypes per locus based on gametes produced from all parents;

        ###collect mismatches from all; 
        poss.ploidy=ploidy
        for (o in 1:nrow(offspring.geno)){
        print(paste0("Now comparing dam #",m,", ",rownames(dam.geno)[m],
                      ", and sire #",p,", ",rownames(sire.geno)[p],
                      ", with offspring #",o,", ", rownames(offspring.geno)[o]))
        genos<-as.list(offspring.geno[o,snps])


        match<-mapply(function(x, y) ifelse(x %in% y,0,1), genos, pot.offspr.genos) #match 0, mismatch=1;
        #exclude combinations with 0
        match[as.logical(as.character(unlist(lapply(pot.offspr.genos, function(x) unique(ifelse(x %in% missing.genos,TRUE,FALSE))))))|genos %in% missing.genos]<-NA
        
        mismatches<-sum(na.omit(match))
        no.surveyed<-sum(!is.na(match))
        perc.mismatches<-mismatches/no.surveyed
        
        locus_mismatches<-rbind(as.matrix(locus_mismatches), 
                                c(rownames(offspring.geno)[o], 
                                  rownames(dam.geno)[m],
                                  rownames(sire.geno)[p],
                                  match))
        summed_mismatches<-rbind(summed_mismatches,
                                 cbind(rownames(offspring.geno)[o],
                                       rownames(dam.geno)[m],
                                       rownames(sire.geno)[p], 
                                       mismatches,no.surveyed,perc.mismatches))
        }
      }
  }
    
colnames(summed_mismatches)<-c("offspring","compared_female_parent","compared_male_parent","mismatches","loci.compared","percent.mismatch")
locus_mismatches<-as.data.frame(locus_mismatches)
colnames(locus_mismatches)<-c("offspring","compared_female_parent","compared_male_parent",snps)
head(summed_mismatches)
dim(summed_mismatches)
locus_mismatches[1:5,1:10]
dim(locus_mismatches)

write.table(summed_mismatches,file="Mendelian_incompatibility_summary.txt",row.names = FALSE,quote=FALSE)
write.table(locus_mismatches, "Mendelian_incompatibility_snps.txt", quote=F, row.names = F)



