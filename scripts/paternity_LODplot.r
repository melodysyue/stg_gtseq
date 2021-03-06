library(tidyverse)

rm(list=ls())

#paternity 
p <- read.table("./data/gtseq/o_parentage_paternity.txt", header=T, sep="\t")
#true.paternity <- read.table("./data/true_parentage.csv", header=T)

p1 <- p %>% 
  select(Offspring.ID, Mother.ID, Candidate.father.ID, Pair.LOD.score.1) %>% 
  group_by(Offspring.ID, Mother.ID) %>% 
  summarize(true=max(Pair.LOD.score.1), fake=min(Pair.LOD.score.1)) %>% 
  gather(father, LOD, -Offspring.ID, -Mother.ID) 


p1$father <- factor(p1$father, levels=c("true", "fake"))

pdf("paternity.pdf", width = 12, height = 9)
p1 %>% 
  ggplot(aes(x=father, y=LOD, fill=father))+
  geom_boxplot(alpha=0.7)+
  theme_bw(base_size = 20)+
  ylab("LOD Score")+
  theme(legend.position = "none",
        axis.title.x = element_blank())+
  scale_fill_brewer(palette = "Set2")+
  scale_x_discrete(labels=c("Confirmed", "Not confirmed"))

dev.off()



