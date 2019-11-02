########################################
#### trophic flow sankey ####
########################################

#### prelims ####
library(tidyverse)
library(networkD3)
library(bipartite)



########################################
#### data ####

## tri-trophic, make her long
tri.long<-read.csv('/Users/collnell/Dropbox/Projects/canada/clean data/canada_cleaned_wide.csv')%>%
  mutate(G_sp = case_when(G_sp == 'Sciaphila duplex' ~ 'Pseudosciaphila duplex', 
                          G_sp == "Anteraea polyphemus" ~ 'Antheraea polyphemus',
                          G_sp == "Campae perlata" ~ 'Campaea perlata',
                          G_sp == 'Ennomos subsignarius' ~ 'Ennomos subsignaria',
                          G_sp == 'Malacosoma americanum' ~ 'Malacosoma americana',
                          G_sp == 'Sparganothis pettitana' ~ 'Cenopis pettitana',
                          G_sp == 'Lambina fiscellaria' ~ 'Lambdina fiscellaria',
                          G_sp == 'Pococera asperatella' ~'Pococera aplastella',
                          TRUE ~ as.character(G_sp)),
         MOTH_family = case_when(G_sp == 'Dioryctria abietivorella' ~ 'Pyralidae', TRUE ~ as.character(MOTH_family)))%>%
  mutate(PTER_sp = case_when(PARA_family == 'Pteromalide' ~ PARA_sp), 
         TACH_sp = case_when(PARA_family == 'Tachinidae' ~ PARA_sp),
         PTER_n = case_when(PARA_family == 'Pteromalide'~ PARA_n),
         HYME_sp = case_when(PARA_order == 'Hymenoptera' ~ PARA_sp_only),
         DIPT_sp = case_when(PARA_order == 'Diptera' ~ PARA_sp_only),
         TREE_sp_only = case_when(!(word(TREE_G_sp, 2) == 'sp.')~TREE_G_sp))%>%
  filter(!(TREE_genus == 'Gleditsia' & TREE_family == 'Sapindaceae'), TREE_genus != 'Viburnum')%>%
  left_join(lep.info, by=c('G_sp'='HERB_sp'))
## PARA_n is the total number of parasitoids for each collection, TACH_n etc is taxa-specific



## add parasitoid count by taxa to wide df
para.nums<-tri.long%>%select(ID, TACH_n:PTER_n)%>%group_by(ID)%>%summarize_at(vars(TACH_n, DIPT_n, HYME_n, BRAC_n, ICH_n, PTER_n), funs(sum(., na.rm=TRUE)))

tri.wide<-read.csv('clean data/canada_cleaned_all.csv')%>%
  mutate(G_sp = case_when(G_sp == "Anteraea polyphemus" ~ 'Antheraea polyphemus',
                          G_sp == "Campae perlata" ~ 'Campaea perlata',
                          G_sp == 'Ennomos subsignarius' ~ 'Ennomos subsignaria',
                          G_sp == 'Malacosoma americanum' ~ 'Malacosoma americana',
                          G_sp == 'Sparganothis pettitana' ~ 'Cenopis pettitana',
                          G_sp == 'Sciaphila duplex' ~ 'Pseudosciaphila duplex', 
                          G_sp == 'Lambina fiscellaria' ~ 'Lambdina fiscellaria',
                          G_sp == 'Pococera asperatella' ~ 'Pococera aplastella',
                          TRUE ~ as.character(G_sp)))%>%
  left_join(para.nums, by='ID')%>%
  mutate(MOTH_host_n = replace_na(Adults,0)+replace_na(PARA_n,0), 
         MOTH_host_tach = replace_na(Adults, 0)+replace_na(TACH_n), 
         MOTH_host_para = replace_na(Adults, 0)+replace_na(PARA_n),MOTH_pre=replace_na(Larvae,0)+replace_na(Pupae,0))%>%
  mutate(TREE_sp = word(TREE_G_sp, 2))%>%
  mutate(TREE_sp_only = case_when(!(TREE_sp %in% c('sp.', 'sp')) ~ TREE_G_sp))%>%
  mutate(MOTH_family = case_when(G_sp == 'Dioryctria abietivorella' ~ 'Pyralidae', TRUE ~ as.character(MOTH_family)))%>%
  filter(!(TREE_genus == 'Gleditsia' & TREE_family == 'Sapindaceae'), TREE_genus != 'Viburnum')%>%
  left_join(lep.info, by=c('G_sp'='HERB_sp'))

tri.out<-tri.wide%>%
  mutate(PARA_n = ifelse(is.na(PARA_n), 0, PARA_n),
         MOTH_larvae = ifelse(is.na(Larvae), Adults + PARA_n - Pupae, Larvae),
         MOTH_pupae = ifelse(is.na(Pupae), 0, Pupae),
         MOTH_host = MOTH_larvae + MOTH_pupae,
         MOTH_host = ifelse(MOTH_host < Adults, Adults, MOTH_host),
         MOTH_host = ifelse(MOTH_host == 0, 1, MOTH_host),
         MOTH_adult = Adults,
         MOTH_out = MOTH_adult+PARA_n,
         MOTH_host= ifelse(PARA_n != 0 & MOTH_host == 0, 1, MOTH_host), 
         MOTH_para = ifelse(MOTH_host < MOTH_out, MOTH_host-MOTH_adult, PARA_n),
         MOTH_para_hosts = MOTH_para+MOTH_adult,
         PARA_rate = MOTH_para/MOTH_para_hosts,
         PARA_rate = ifelse(PARA_rate == 'NaN', NA, PARA_rate))

## PLANT -HERBIVORE MATRICES
tri.ph.wide<-tri.out%>%
  group_by(HERB_family_2019, HERB_ID_SP, TREE_clade, TREE_family, TREE_genus)%>%
  summarize(HERB_n = sum(MOTH_host))%>%
  dcast(TREE_clade+TREE_family+TREE_genus~HERB_ID_SP, value.var='HERB_n', fill=0)

tri.ph.mat<-tri.ph.wide%>%column_to_rownames('TREE_genus')%>%dplyr::select(-TREE_clade, -TREE_family)
tri.ph.prop<-tri.ph.mat%>%decostand('total')

## HOST PARASITOPID MATRICES
tri.hp.long<-tri.long%>%
  left_join(tri.out%>%dplyr::select(ID, HERB_ID_SP, MOTH_host, MOTH_adult, MOTH_para, MOTH_para_hosts, PARA_rate))%>%
  group_by(ID, HERB_ID_SP, TREE_genus, PARA_order, PARA_family, PARA_sp, PARA_sp_only, MOTH_host, MOTH_adult, MOTH_para, MOTH_para_hosts, PARA_rate)%>%
  summarize(PARA_n_sp = sum(PARA_n))%>%
  filter(!is.na(PARA_order))

tri.hp.wide<-tri.long%>%
  left_join(tri.out%>%dplyr::select(ID, HERB_ID_SP, MOTH_host, MOTH_adult, MOTH_para, MOTH_para_hosts, PARA_rate))%>%
  group_by(HERB_ID_SP, PARA_order, PARA_family, PARA_sp_only)%>%
  summarize(PARA_n = sum(PARA_n), HERB_n = sum(MOTH_para_hosts))%>%
  filter(!is.na(PARA_sp_only), !is.na(HERB_ID_SP))%>%
  dcast(HERB_ID_SP~PARA_sp_only, value.var='PARA_n', fill=0)

hyme.mat<-tri.long%>%
  left_join(tri.out%>%dplyr::select(ID, HERB_ID_SP, MOTH_host, MOTH_adult, MOTH_para, MOTH_para_hosts, PARA_rate))%>%
  group_by(HERB_ID_SP, PARA_order, PARA_family, HYME_sp)%>%
  summarize(PARA_n = sum(HYME_n), HERB_n = sum(MOTH_para_hosts))%>%
  filter(!is.na(HYME_sp), !is.na(HERB_ID_SP))%>%
  dcast(HERB_ID_SP~HYME_sp, value.var='PARA_n', fill=0)
hyme.mat


dipt.mat<-tri.long%>%
  left_join(tri.out%>%dplyr::select(ID, HERB_ID_SP, MOTH_host, MOTH_adult, MOTH_para, MOTH_para_hosts, PARA_rate))%>%
  group_by(HERB_ID_SP, PARA_order, PARA_family, DIPT_sp)%>%
  summarize(PARA_n = sum(DIPT_n), HERB_n = sum(MOTH_para_hosts))%>%
  filter(!is.na(DIPT_sp), !is.na(HERB_ID_SP))%>%
  dcast(HERB_ID_SP~DIPT_sp, value.var='PARA_n', fill=0)

## bi trophic get em up

bi.long<-read.csv('/Users/collnell/Dropbox/Projects/canada/MM270.csv')%>%filter(!(PLANT_genus %in% c('Caragana','Ribes','Rosa','Shepherdia')))

# matrix of host plant associations of focal herbivores - aggregate by new species ids
bi.ph.wide<-bi.long%>%
  group_by(HERB_family_2019, HERB_sp_only, HERB_ID_SP, PLANT_genus)%>%
  summarize(n=sum(Records, na.rm=TRUE))%>%
  dcast(HERB_family_2019+HERB_sp_only+HERB_ID_SP~PLANT_genus, value.var='n', fill=0)

bi.ph.mat<-bi.ph.wide%>%column_to_rownames('HERB_ID_SP')%>%dplyr::select(-HERB_family_2019, -HERB_sp_only)

### filter ot focal species
# keep leps with 50 records, tree genera with 100
bi.trees<-colnames(bi.ph.mat[colSums(bi.ph.mat)>=50])# 27 w/ 100 leps, 32 if filtered at 50

## community matrices, relative abundances
tri.hp.mat<-tri.hp.wide%>%column_to_rownames('HERB_ID_SP')
tri.hp.prop<-tri.hp.mat%>%decostand('total')
bi.mat<-bi.ph.mat%>%dplyr::select(bi.trees)%>%t()
bi.mat.prop<-bi.mat%>%decostand('total')

# store tree names
bi.trees<-rownames(bi.mat)
########################################


