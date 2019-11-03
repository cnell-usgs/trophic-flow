########################################
#### trophic flow sankey ####
########################################

#### prelims ####
library(tidyverse);library(reshape2)
library(networkD3);library(bipartite)

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

# Load energy projection data
URL <- "https://cdn.rawgit.com/christophergandrud/networkD3/master/JSONdata/energy.json"
Energy <- jsonlite::fromJSON(URL)

# 2 dataframe - links and nodes
Energy$links
Energy$nodes

p <- sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
                   Target = "target", Value = "value", NodeID = "name",
                   units = "TWh", fontSize = 12, nodeWidth = 30)
p

# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(
  source=c("group_A","group_A", "group_B", "group_C", "group_C", "group_E"), 
  target=c("group_C","group_D", "group_E", "group_F", "group_G", "group_H"), 
  value=c(2,3, 2, 3, 1, 3)
)
links

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)
nodes

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p

###########################################################
## generate edgelist

bi.feed.long<-bi.mat%>%as.data.frame()%>%
  rownames_to_column('TREE_genus')%>%
  melt(variable.name='HERB_ID_SP')%>%
  filter(value>0)%>%
  left_join(lep.info)%>%
  mutate(feed=str_trim(feed), TREE_class=case_when(TREE_genus %in% c('Abies','Picea','Pinus','Tsuga','Thuja','Pseudotsuga','Juniperus', 'Larix') ~ 'Gymnosperm',TRUE ~ 'Angiosperm'))
bi.feed.long

edge.ph<-bi.feed.long%>%
  group_by(TREE_class,TREE_genus, feed)%>%
  summarize(n=sum(value), sp=length(unique(HERB_ID_SP)), fam=length(unique(HERB_family_2019)))
edge.ph

node.ph<-data.frame(name=c(as.character(edge.ph$TREE_genus), as.character(edge.ph$feed))%>%unique())

edge.ph$source<-match(edge.ph$TREE_genus, node.ph$name)-1
edge.ph$target<-match(edge.ph$feed, node.ph$name)-1


p<-sankeyNetwork(Links = as.data.frame(edge.ph)%>%filter(n>8), Nodes = node.ph,
              Source = 'source', Target='target', Value = 'sp', NodeID = 'name', sinksRight=FALSE)

p


edge.ph%>%filter(n>8)%>%view
## look at species list for feeding guilds
bi.feed.long%>%
  group_by(feed, HERB_family_2019)%>%
  summarize(sp_n=length(unique(HERB_sp_2019)), sps=paste0(unique(HERB_sp_2019), collapse=', '))%>%View

cat_n<-sum(edge.ph$n)
cat_n

feed.sp<-bi.feed.long%>%group_by(feed, HERB_sp_2019, TREE_class)%>%
  summarize(n=sum(value), genera=length(unique(TREE_genus)))%>%
  filter(n>10)%>%dcast(feed+HERB_sp_2019~TREE_class, value.var='n')%>%
  mutate(diet = case_when(is.na(Gymnosperm) ~ 'Angiosperm', is.na(Angiosperm) ~ 'Gymnosperm', !is.na(Angiosperm) & !is.na(Gymnosperm) ~ 'both'))
feed.sp

feed.sp%>%
  group_by(diet)%>%
  summarize(sp=length(unique(HERB_sp_2019))) ## 42 crossover species

# how many species of each diet category are in each feeding guild?
str(feed.sp)
feed.ga<-feed.sp%>%
  replace_na(list(Angiosperm=0, Gymnosperm=0))%>%
  mutate(n=Angiosperm+Gymnosperm)%>%
  group_by(feed, diet)%>%
  summarize(n_sum=sum(n, na.rm=TRUE), sp=length(unique(HERB_sp_2019)))
feed.ga

plotweb()

bi.feed.cast<-edge.ph%>%
  mutate(n_prop=n/cat_n)%>%
  dcast(TREE_class~feed, value.var='sp', fill=0)%>%
  mutate()
bi.feed.cast

low.n<-rowSums(bi.feed.cast%>%column_to_rownames('TREE_class'))
#low.n[1:2]<-1
low.n
mat<-bi.feed.cast%>%column_to_rownames('TREE_class')
plotweb(bi.feed.cast%>%column_to_rownames('TREE_class'),
        low.abun=colSums(mat), high.abun=rowSums(mat))
colSums(mat)

library(igraph)

mat/rowSums(mat)
rowSums(mat)

# first make herb feeding guilds
100*(log(1+colSums(bi.feed.cast[,-1]))/24.23641)

sum(log(1+colSums(bi.feed.cast[,-1]))) #327 total
 

log(1+colSums(bi.feed.cast[,-1]))

## scale plants relative to diversity on them
## herb circle relative diversity
## edge width is diversity, color proportional to host class cat


# sam ebut proportion data?
View(lep.info)

## sclae tree trophic level to total number of caterpillars associated with each
# or total diversty


# herbivore aslo scaled to total number
# vary sizes based on abundance


## width of edge is the total number of caterpiilars
# color within by whether or not they feed on one or both
# so need new category within each ofproportion feeding one or either


#######################################

#### host-parasitoid ####
