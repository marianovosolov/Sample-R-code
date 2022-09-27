#### This code generates the backbone to the figure of the Danish Ancestry in the Meso-Neo project run by Morten Allentoft and Martin Sikora.
#### The figure can be found as figure 4 in the preprint https://www.biorxiv.org/content/10.1101/2022.05.04.490594v2.full


library(tidyverse)
library(ggvis)
library(ggthemes)
library(ggprism)
library(rcarbon)
# sample_list_full<- read_delim("data/intersect_ancestry_isotope.txt",delim = "\t",col_names = F)

imput_list<- read_delim("data/intersect_ancestry_isotope_imputed.txt",delim = "\t",col_names = F)
###### ancestry ########
anc_data_full<- read_delim("data/dk_ph_4way_yam.admixture_forPlot.csv",delim = ",")

sample_list_sr_diet<- read_csv("data/Meso_Neo_data_11_forPlot_Denmark.csv") %>% 
  janitor::clean_names()
unique(anc_data_full$flag)
#### update sample list ####

sample_list_full_t<- sample_list_sr_diet %>% 
  rename(sampleId = ris_eno) %>% 
  select(sampleId,bp_date) %>% 
  mutate(sampleId = fct_reorder(sampleId,bp_date,.desc = T)) %>% 
  write_csv(.,"data/sample_list_191021.csv")

#### Diet and strontium #########

ds_full<- sample_list_sr_diet %>% 
  mutate(ris_eno = fct_reorder(ris_eno,bp_date,.desc = T))

diet_full_p1<- ggplot(ds_full,aes(x=ris_eno),na.rm=F)+
  geom_point(aes(y=x13cdiet),na.rm = F)+
  geom_vline(xintercept = 33.25,lwd = 1.5)+
  geom_vline(xintercept = 68.25,lwd = 1.5)+
  labs(x=NULL, y=NULL, title=NULL)+
  coord_equal(ratio = 3.5/4,expand = T)+
  theme_minimal()+
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
  axis.ticks.y = element_blank(),
  legend.title = element_text(size=8))
diet_full_p1

diet_full_p2<- ggplot(ds_full,aes(x=ris_eno),na.rm=F)+
  geom_point(aes(y=x15ndiet),na.rm = F)+
  geom_vline(xintercept = 33.25,lwd = 1.5)+
  geom_vline(xintercept = 69.25,lwd = 1.5)+
  labs(x=NULL, y=NULL, title=NULL)+
  coord_equal(ratio = 4.5/4,expand = T)+
  theme_minimal()+
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size=8))
diet_full_p2

# plot strontium
stro_full_p<- ggplot(ds_full,aes(ris_eno,sr86_sr87))+
  geom_point(na.rm = F)+
  labs(x=NULL, y=NULL, title=NULL)+
  geom_vline(xintercept = 33.25,lwd = 1.5)+
  geom_vline(xintercept = 68.25,lwd = 1.5)+
  coord_equal(ratio = 5500/4,expand = T)+
  theme_minimal()+
  theme(axis.ticks.x=element_blank(),axis.text.x = element_blank())+
  theme(legend.title = element_text(size=8))
stro_full_p

cowplot::plot_grid(diet_full_p1,diet_full_p2,stro_full_p,ncol=1)

##### Plot Ancestry ######

anc_data_full_filt<- anc_data_full %>% 
  filter(value!=0) %>% 
  mutate(sampleId = case_when(
    nchar(str_extract(sampleId,"[0-9]+"))==1~paste0("NEO00",str_extract(sampleId,"[0-9]+")),
    nchar(str_extract(sampleId,"[0-9]+"))==2~paste0("NEO0",str_extract(sampleId,"[0-9]+")),
    TRUE~sampleId)) %>% 
  mutate(sampleId = case_when(sampleId=="NEO962"~"NEO962b",TRUE~sampleId)) %>% 
  filter(sampleId!="NEO957") %>% 
  mutate(grey = case_when(coverage < 0.05 ~"yes",TRUE ~"no")) %>% 
  left_join(sample_list_full_t,by = "sampleId") %>% 
  filter(!is.na(value)) %>% 
  rename(popId = sampleId) %>% 
  mutate(popId = fct_reorder(popId,bp_date,.desc = T)) %>% 
  rename(Var2 = name)
unique(anc_data_full_filt$flag)

anc_data_full_filt %>% 
  group_by(flag) %>% 
  tally()


unique(anc_data_full$group)

anc_data_full_filt<- anc_data_full_filt %>%
  mutate(Var2 = factor(Var2, levels = c("X1","X2","X3","X4")))


anc_cols_new<- c("#698B22", #army green
                 "#CD1076", #purple
                 "#8B7355", #brown
                 "#FF6347" #orange
                 )

anc_full_p<- ggplot(anc_data_full_filt,aes(x = popId,y = value,fill = Var2,alpha = grey))+
  geom_col(size = 0.1,width = 1,na.rm = F)+
  theme_minimal() +
  scale_fill_manual(values = anc_cols_new)+
  labs(x=NULL, y=NULL, title=NULL)+
  geom_vline(xintercept = 34.55,lwd = 1.5)+
  geom_vline(xintercept = 69.55,lwd = 1.5)+
  # labs(x = "Individuals", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_alpha_manual(values = c(1,0.8))+
  coord_equal(ratio = 80/6,expand = T)+
  scale_x_discrete(expand = expand_scale(add = 1),position = "top") +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x.top = element_text(angle = 90,vjust = 0.5),
    panel.grid = element_blank(),legend.position = "none")
  
anc_full_p
# length(unique(anc_data_full_filt$popId))
# cowplot::plot_grid(anc_full_p,eye.color.p,hair.color.p,height.p,ncol = 1,align = "v")

### Phenotype andres grayed out####
pheno.andres<- read_csv("data/phenotype_withLowCov_Andres_290921.csv") %>% 
  mutate(sampleId = case_when(
    nchar(str_extract(sampleId,"[0-9]+"))==1~paste0("NEO00",str_extract(sampleId,"[0-9]+")),
    nchar(str_extract(sampleId,"[0-9]+"))==2~paste0("NEO0",str_extract(sampleId,"[0-9]+")),
    TRUE~sampleId)) %>% 
  mutate(sampleId = case_when(sampleId=="NEO962"~"NEO962b",TRUE~sampleId)) %>% 
  mutate(cov_level = case_when(coverage < 0.1~"yes",TRUE~"no"))

pheno.andres.complete<- sample_list_full_t %>% 
  left_join(pheno.andres,by = "sampleId") %>% 
  # rename(sampleId=pop) %>% 
  mutate(sampleId = fct_reorder(sampleId,bp_date))

# length(unique(pheno.andres.complete$sampleId))

# setdiff(unique(anc_data_full_filt$popId),unique(pheno.andres.complete$sampleId))
eye.andres<- pheno.andres.complete %>% 
  select(sampleId,pigmQC,bp_date,cov_level,contains("Eye")) %>% 
  pivot_longer(cols = -c(sampleId,pigmQC,bp_date,cov_level),names_to = "eye_color",values_to = "values") %>% 
  mutate(sampleId = fct_reorder(sampleId,bp_date,.desc = T))

hair.andres<- pheno.andres.complete %>% 
  select(sampleId,pigmQC,bp_date,cov_level,contains("Hair")) %>% 
  pivot_longer(cols = -c(sampleId,pigmQC,bp_date,cov_level),names_to = "hair_color",values_to = "values") %>% 
  mutate(sampleId = fct_reorder(sampleId,bp_date,.desc = T))

height.andres<- pheno.andres.complete %>% 
  select(sampleId,pigmQC,bp_date,cov_level,pHeight) %>% 
  mutate(sampleId = fct_reorder(sampleId,bp_date,.desc = T))


#phenotype andres plots
height.p<- ggplot(height.andres,aes(x=sampleId,alpha=cov_level),na.rm=F)+
  geom_point(aes(y=pHeight),na.rm = F)+
  # geom_vline(xintercept = 25.25,lwd = 1.5)+
  # geom_vline(xintercept = 42.25,lwd = 1.5)+
  labs(x=NULL, y=NULL, title=NULL)+
  coord_equal(ratio = 4/4,expand = T)+
  theme_minimal()+
  theme(axis.ticks.x=element_blank(),axis.text.x = element_blank())+
  theme(axis.ticks.y = element_blank(),legend.position = "none")+
  theme(legend.title = element_text(size=8))+
  scale_alpha_manual(values = c(1,0.6))
height.p

eye.cols<- c("#a0a8a5","#5c350e","#2390de")
eye.color.p<- ggplot(eye.andres,aes(x =sampleId,y = values,fill = forcats::fct_rev(eye_color),alpha = cov_level))+
  geom_col(size = 0.1,width = 1,na.rm = F)+
  theme_minimal() +
  labs(x=NULL, y=NULL, title=NULL)+
  # geom_vline(xintercept = 25.25,lwd = 1.5)+
  # geom_vline(xintercept = 42.25,lwd = 1.5)+
  # labs(x = "Individuals", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  coord_equal(ratio = 40/4,expand = T)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),legend.position = "none")+
  scale_fill_manual(values = eye.cols)+
  scale_alpha_manual(values = c(1,0.7))
eye.color.p

hair.cols<- c("#f58700","#692c01","#f5df4e","#242120")
hair.color.p<- ggplot(hair.andres,aes(x = sampleId,y = values,fill = forcats::fct_rev(hair_color),alpha = cov_level))+
  geom_col(size = 0.1,width = 1,na.rm = F)+
  theme_minimal() +
  labs(x=NULL, y=NULL, title=NULL)+
  # geom_vline(xintercept = 25.25,lwd = 1.5)+
  # geom_vline(xintercept = 42.25,lwd = 1.5)+
  # labs(x = "Individuals", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  coord_equal(ratio = 40/4,expand = T)+
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),legend.position = "none")+
  scale_fill_manual(values = hair.cols)+
  scale_alpha_manual(values = c(1,0.7))
hair.color.p
cowplot::plot_grid(anc_full_p,eye.color.p,hair.color.p,height.p,ncol = 1,align = "v")

##### timeline #######
time_data_full<- sample_list_full_t %>% 
  # select(pop,age_time_bp) %>% 
  # distinct() %>% 
  mutate(y_axis = 0) %>% 
  mutate(pop = fct_reorder(sampleId,bp_date))
# mutate(y_axis2 = 0.1) %>% 
# pivot_longer(cols = c(y_axis,y_axis2),values_to = "value",names_to = "axis")

time_vec<- seq(from = min(time_data_full$bp_date),to = max(time_data_full$bp_date,by = 129))

time_full_p<-
  ggplot(time_data_full,aes(x=bp_date,y=y_axis))+
  geom_point()+
  geom_line()+
  scale_x_reverse(expand = c(0, 0),limits = c(10950,2450),guide = "prism_minor",minor_breaks = seq(2450,10950, 500))+
  scale_y_continuous(expand = c(0,0.5))+
  ggrepel::geom_text_repel(aes(label = sampleId),angle = 90,size = 3.5,nudge_y = 0.3,
                           force             = 0.5,
                           direction         = "x",
                           hjust             = 0,
                           segment.size      = 0.2)+
  #,limits = c(11000,500),guide = "prism_minor",minor_breaks = seq(500,11000, 500)
  # coord_equal(1/5000)+
  # ylim(0.1,-0.1)+
  # theme_classic()+
  labs(x=NULL, y=NULL, title=NULL)+
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_rect("white"),
    legend.position = "none")
time_full_p

###### pollen #########

pollen<- read_csv("data/pollen_veg_composition_new.csv") %>% 
  janitor::clean_names() %>% 
  # select(-age_bc) %>% 
  pivot_longer(cols = -c(age_bp,age_bc),names_to = "veg",values_to = "values") %>% 
  mutate(veg = factor(veg,levels = c("crops","grassland","secondary_forest","forest")))


# col vectors

pollen_cols<- c("#d8d678","#f4824e","#87c440","#138742")

pollen_p_v2<- ggplot(pollen,aes(age_bp,values,fill = veg))+
  geom_area(color="black")+
  scale_x_reverse(expand = c(0, 0),guide = "prism_minor",minor_breaks = seq(4493, 6905, 100))+
  scale_fill_manual(values = pollen_cols)+
  geom_vline(xintercept = 4750,lwd = 1.5)+
  geom_vline(xintercept = 5854,lwd = 1.5)+
  coord_equal(ratio = 10/4,expand = T)+
  theme(panel.background = element_rect("white"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = )+
  # scale_x_continuous(expand = c(7000, 0)) + 
  scale_y_continuous(expand = c(0, 0))
pollen_p_v2


#### Y major plot #####

ychr_data<- read_delim("data/neo.impute.1000g.prune.admixture.Q.pong_new.tsv",delim = "\t") %>% 
  filter(maxK==7) %>% 
  filter(country=="Denmark") %>%
  filter(value!=0) %>% 
  filter(str_detect(popId,"NEO")) %>% 
  mutate(popId = case_when(
    nchar(str_extract(popId,"[0-9]+"))==1~paste0("NEO00",str_extract(popId,"[0-9]+")),
    nchar(str_extract(popId,"[0-9]+"))==2~paste0("NEO0",str_extract(popId,"[0-9]+")),
    TRUE~popId)) %>% 
  mutate(popId = case_when(popId=="NEO962"~"NEO962b",TRUE~popId)) %>% 
  mutate(ageAverage = case_when(popId=="NEO962b"~5785,TRUE~ageAverage)) %>% 
  mutate(Var2 = as.character(Var2)) %>% 
  filter(!is.na(value)) %>% 
  mutate(impute_samp = case_when(popId%in%imput_list$X1~"impute",TRUE~"not_impute")) %>% 
  right_join(sample_list_full_t,by = c("popId" = "sampleId")) %>% 
  select(popId,bp_date,hgYMajor,hgMT,sex) %>% 
  mutate(popId = fct_reorder(popId,bp_date,.desc = T)) %>% 
  mutate(hgYMajor=factor(hgYMajor,levels=c("I2","R1b","I1","Q1","L1","I","R1a"))) %>% 
  distinct()

Ycolors<- read_csv("data/hginfo2_updated.csv") %>% 
  filter(hgId%in%ychr_data$hgYMajor) %>% 
  arrange(factor(hgId, levels = unique(ychr_data$hgYMajor)))

length(unique(ychr_data$popId))

ychr_plot<- ggplot(ychr_data, aes(popId,forcats::fct_rev(hgYMajor),fill=hgYMajor))+
  geom_tile(na.rm = F)+
  scale_fill_manual(values = Ycolors$hex_colors)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())+
  coord_equal(ratio = 4/4,expand = T)+
  theme(panel.background = element_rect("white"),legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )+
  labs(x=NULL, y=NULL, title=NULL)+
  scale_x_discrete(breaks = NULL)

ychr_plot
###### MT plot ######
MT_data<- ychr_data %>% 
  # filter(!is.na(hgMT)) %>%
  mutate(MT=str_extract(hgMT,"[A-Z]")) %>% 
  mutate(MT=factor(MT,levels=c("U","H","K","R","J","N","T","V","W"))) %>% 
  mutate(popId = fct_reorder(popId,bp_date,.desc = T))

MT_colors<- c("#cd2990",rep("#000000",8))

MT_plot<-  ggplot(MT_data,aes(popId,forcats::fct_rev(MT),fill=MT))+
  geom_tile(na.rm = F)+
  theme(axis.text.x = element_text(angle = 90),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())+
  coord_equal(ratio = 4/4,expand = T)+
  theme(panel.background = element_rect("white"),legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )+
  labs(x=NULL, y=NULL, title=NULL)+
  scale_fill_manual(values = MT_colors)

MT_plot


#### Sample sex ######


sex_data<- ychr_data %>% 
  select(popId,sex,bp_date) %>% 
  mutate(popId = fct_reorder(popId,bp_date,.desc = T))

sex_plot<-  ggplot(sex_data,aes(popId,sex))+
  geom_tile(na.rm = F)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())+
  coord_equal(ratio = 5/4,expand = T)+
  theme(panel.background = element_rect("white"),legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  )+
  labs(x=NULL, y=NULL, title=NULL)


sex_plot

######## cultural succession #########3

cs<- read_csv("data/cultural_sussession.csv") %>% 
  mutate(culture = factor(culture,levels = c("Maglemose",
                                             'Kongemose',
                                             'Ertebølle',
                                             "Funnel Beaker", 
                                             "Pitted Ware",
                                             "Single Grave",
                                             "Dagger",
                                             "Bronze Age"))) %>% 
  select(-c(start_BC,end_BC)) %>% 
  pivot_longer(cols=-c(culture,item),names_to = "state",values_to = "date")


cs_plot<- ggplot(cs,aes(date, forcats::fct_rev(culture), group=item)) +
  geom_line(size = 2) +
  coord_equal(ratio = 400/4,expand = T)+
  scale_x_reverse(expand = c(0, 0),
                  guide = "prism_minor",
                  minor_breaks = seq(0, 9000, 500),
                  breaks = seq(min(cs$date),max(cs$date),by=1000))+
  theme_classic()
cs_plot  

######## merge all together #########
cowplot::plot_grid(cs_plot,time_full_p,anc_full_p,sex_plot,ychr_plot,MT_plot,eye.color.p,hair.color.p,height.p,stro_full_p,diet_full_p1,diet_full_p2,pollen_p_v2, ncol = 1,align = "v")
