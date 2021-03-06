library(tidyverse)
library(ggpubr)

rm(list=ls())
summed_mismatches <- read.table("Mendelian_incompatibility_summary.txt",header=T)


parentage <- read.csv("./data/gtseq/parentage.csv", header=T)
head(summed_mismatches)
dim(summed_mismatches)

parentage.mis <- inner_join(summed_mismatches, parentage, 
                            by=c("offspring"="Sample", "compared_female_parent"="Female_parent", "compared_male_parent"="Male_parent"))

#for those with uncertain male parents, male 1 or male 2, the one with smaller micmatches is the true male parent

both.true <- parentage.mis %>% 
  group_by(offspring, compared_female_parent) %>% 
  top_n(-1, percent.mismatch) %>% 
  select(Sample=offspring, Female_parent=compared_female_parent, Male_parent=compared_male_parent)

both.true <- as.data.frame(both.true)

write.table(both.true, "./true_parentage.csv", quote=F, row.names = F)


###MI for 2 true parents
both.mis <- inner_join(summed_mismatches, both.true, by=c("offspring"="Sample", 
                                                          "compared_female_parent"="Female_parent", 
                                                          "compared_male_parent"="Male_parent"))


###MI for 1 true parent
one.mis <- NULL
for (i in 1:nrow(both.true)){
  o=as.character(both.true[i,"Sample"])
  m=as.character(both.true[i, "Female_parent"])
  p=as.character(both.true[i, "Male_parent"])
  
  temp1 <- summed_mismatches %>% 
    filter(offspring==o & compared_female_parent==m & compared_male_parent!= p)
  
  temp2 <- summed_mismatches %>% 
    filter(offspring==o & compared_female_parent!=m & compared_male_parent==p )
  
  one.mis <- rbind(one.mis, temp1, temp2)
}


###MI for 0 true parent
zero.mis <- anti_join(summed_mismatches, both.mis) 
zero.mis <- anti_join(zero.mis, one.mis)

df <- rbind(both.mis, one.mis, zero.mis)
df$type <- c(rep("Two", nrow(both.mis)), rep("One", nrow(one.mis)), rep("None", nrow(zero.mis)))
df$type <- as.factor(df$type)
df$type <- factor(df$type, levels=c("Two", "One", "None"))


df %>% 
  group_by(type) %>% 
  summarise(mean=mean(percent.mismatch), sd=sd(percent.mismatch))

pdf(".MI.pdf", width = 12, height = 9)
df %>% 
  ggplot(aes(x=percent.mismatch, fill=type))+
  geom_density(alpha=0.7)+
  theme_bw(base_size = 20)+
  xlab("Percentage of Mendelian Incompatibilities")+
  ylab("Density")+
  theme(legend.position = "none")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~type, ncol=1)

df %>% 
  ggplot(aes(x=type, y=percent.mismatch, fill=type))+
  geom_boxplot(alpha=0.7)+
  theme_bw(base_size = 20)+
  ylab("Percentage of Mendelian Incompatibilities")+
  theme(legend.position = "none",
        axis.title.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  geom_signif(comparisons = list(c("Two", "One"), c("Two","None"), c("One","None")), 
              map_signif_level=TRUE)

dev.off()


