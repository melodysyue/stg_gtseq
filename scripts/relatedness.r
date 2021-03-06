library(tidyverse)

rm(list=ls())


r <- read.table("./data/gtseq/o_relatedness_edited.txt", check.names = F, header=T)
r[1:5,1:5]

r <- r %>% 
  gather(ID2, relatedness, -ID)


#self
self <- r %>% 
  filter(ID==ID2)

#parent-offspring
parentage <- read.table("./data/gtseq/true_parentage.txt", header=T)
parentage <- as.data.frame(parentage)


dam_offspring=left_join(parentage,r,by=c("Sample"="ID", "Female_parent"="ID2")) %>% 
  select(ID=Sample, ID2=Female_parent, relatedness)

sire_offspring=left_join(parentage,r,by=c("Sample"="ID", "Male_parent"="ID2")) %>% 
  select(ID=Sample, ID2=Male_parent, relatedness)

parent_offspring <- rbind(dam_offspring, sire_offspring)

#Full-sib 
relationship <- read.csv("./data/gtseq/parentage.csv", header=T)

full_group1 <- relationship %>% 
  group_by(Female_parent, Male_parent) %>% 
  summarize(n=n()) %>% 
  group_by(Female_parent) %>% 
  mutate(n_male=n()) %>% 
  filter(n_male==1) %>% 
  as.data.frame()

full1 <- NULL

for (i in 1:nrow(full_group1)){
  m=full_group1[i, "Female_parent"]
  p=full_group1[i, "Male_parent"]
  
  o <- relationship %>% 
    filter(Female_parent==m & Male_parent==p) %>% 
    pull(Sample)
  
  temp <- r %>% 
    filter(ID %in% o & ID2 %in% o & ID != ID2)
  
  full1 <- rbind(full1, temp)
}

full_group2 <- parentage %>% 
  group_by(Female_parent, Male_parent) %>% 
  summarize(n=n()) %>% 
  as.data.frame()

full2 <- NULL

for (i in 1:nrow(full_group2)){
  m=full_group2[i, "Female_parent"]
  p=full_group2[i, "Male_parent"]
  
  o <- parentage %>% 
    filter(Female_parent==m & Male_parent==p) %>% 
    pull(Sample)
  
  temp <- r %>% 
    filter(ID %in% o & ID2 %in% o & ID != ID2)
  
  full2 <- rbind(full2, temp)
}

full <- rbind(full1, full2)
full$ID <- as.character(full$ID)
full$ID2 <- as.character(full$ID2)

#remover edundant pairs;
full_nodup <- full %>% 
  mutate(normalized = purrr::map2_chr(ID, ID2, ~paste(sort(c(.x, .y)), collapse = " "))) %>%
  group_by(normalized) %>% 
  summarise(ID = dplyr::first(ID),
            ID2 = dplyr::first(ID2),
            relatedness=mean(relatedness)) %>%
  select(-normalized)


#half siblings, same moms, but different dads;
parentage %>% 
  group_by(Female_parent, Male_parent) %>% 
  summarize(n=n())

moms <-c("18BLA126", "18BLA163")

half_group <- parentage %>% 
    group_by(Female_parent, Male_parent) %>% 
    summarize(n=n()) %>% 
    filter(Female_parent %in% moms)

half <- NULL

for (i in 1:length(moms)){
  
  f=moms[i]
  p=half_group %>% 
    filter(Female_parent == f) %>% 
    pull(Male_parent)
  
  o1 <- parentage %>% 
    filter(Female_parent==f & Male_parent==p[1]) %>% 
    pull(Sample)
  o2 <- parentage %>% 
    filter(Female_parent==f & Male_parent==p[2]) %>% 
    pull(Sample)
  
  temp <- r %>% 
    filter(ID %in% o1 & ID2 %in% o2 & ID != ID2)
  
  half <- rbind(half, temp)
}

#

dim(self)
dim(parent_offspring)
dim(full_nodup)
dim(half)
all <- rbind(self, parent_offspring, full_nodup, half)


all$type <- c(rep("Self", nrow(self)),
              rep("Parent-Offspring", nrow(parent_offspring)),
              rep("Full-sibling", nrow(full_nodup)),
              rep("Half-sibling", nrow(half)))

all$type <- factor(all$type, levels=c("Self","Parent-Offspring", "Full-sibling", "Half-sibling"))

summary(full_nodup)


pdf("./relatedness.pdf", width = 12, height = 9)
all %>% 
  ggplot(aes(x=type, y=relatedness, fill=type))+
  geom_boxplot(alpha=0.7)+
  theme_bw(base_size = 20)+
  ylab("Relatedness")+
  theme(legend.position = "none",
        axis.title.x = element_blank())+
  scale_fill_brewer(palette = "Set2")

dev.off()

