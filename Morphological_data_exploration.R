#### This script explores the patternes in lizaard morphological traits between the two islands Pod Kopiste and Pod Mrcaru.

library(tidyverse)
library("FactoMineR")
library("factoextra")

#Read the pop data
pop.data<- read_csv("data/pop_data_with_bamNames.csv") %>% 
  dplyr::select(-X1)
#read the morpho data and change the sex to be more informative, and filter to have only the ind that we have data for
morpho.data<- read_csv("data/morpho_data/morpho_data_podarcis_sic_MFA.csv") %>% 
  filter(ind_name%in%pop.data$user_id)


head(morpho.data)

no.tail<- morpho.data %>% 
  dplyr::select(-c(real_tail_length,total_tail_length,locality,ind_name))
## plot pairs to see autocorelation
only.tail<- morpho.data %>% 
  dplyr::select(c(real_tail_length,total_tail_length)) 
pairs(only.tail) # leave both of the tail traits for GWAS

only.body<- morpho.data %>% 
  dplyr::select(c(3:6))
pairs(only.body) # leave bodyh , bodyw

only.head<- morpho.data %>% 
  dplyr::select(c(7:15))

pairs(only.head)# leave open, bite force, head hight

only.limbs<- morpho.data %>% 
  dplyr::select(c(18:28))

pairs(only.limbs)# leave everyhitng except forell and hindll, tibia, lth
#### MFA for the morpho data ####

my.groups<- c(4,9,11)
my.type<- c("s","s","s")
my.names<- c("body","head","limbs")
res.mfa <- MFA(no.tail,
               group = my.groups,
               type = my.type,
               name.group = my.names,
               graph = FALSE)

# plots
eig.val <- get_eigenvalue(res.mfa)
head(eig.val)

fviz_screeplot(res.mfa)
group <- get_mfa_var(res.mfa, "group")
group
# plot the groups
fviz_mfa_var(res.mfa, "group")
# contribution plots
fviz_contrib(res.mfa, "group", axes = 1)
fviz_contrib(res.mfa, "group", axes = 2)

### Extract the data for GWAS ####

sub.morpho.data<- morpho.data %>% 
  dplyr::select(ind_name,locality,bodyh, bodyw,open,bite_force,head_hight,real_tail_length,
                total_tail_length,interl,femur,metatars,hindll,humerus,radius,metacar,ltf) %>% 
  left_join(pop.data,by = c("ind_name"= "user_id")) %>% 
  dplyr::select(-c(q30,mreads,ngi_id,island)) %>% 
  dplyr::select(ID_name,ind_name:ltf)



write_csv(sub.morpho.data,"data/morpho_data/GWAS_morpho_data040620.csv")
##### Morpho2 - all the ones that were excluded before #####
sub.morpho.data2<- morpho.data %>% 
  dplyr::select(-c(locality,bodyh, bodyw,open,bite_force,head_hight,real_tail_length,
                total_tail_length,interl,femur,metatars,hindll,humerus,radius,metacar,ltf)) %>% 
  left_join(pop.data,by = c("ind_name"= "user_id")) %>% 
  dplyr::select(-c(q30,mreads,ngi_id,island)) %>% 
  dplyr::select(-ID_name)
write_csv(sub.morpho.data2,"data/morpho_data/GWAS_morpho_data_round2_090620.csv")
#load the diversity data
alpha.div.df<- read.csv("data/alpha_div_per_ind.csv")

alpha.div.df$evenness<- alpha.div.df$alpha.div/12959808
##### ANOVA for the different traits that were chosen for GWAS #####

nest.data<- morpho.data %>% 
  left_join(alpha.div.df,by = c("ind_name" = "ind.id")) %>% 
  dplyr::select(-c(pop,evenness)) %>% 
  mutate(locality = gsub(" ","_",locality)) %>% 
  pivot_longer(cols = -c(ind_name,locality,alpha.div),names_to = "trait",values_to = "value") %>% 
  group_by(trait) %>% 
  nest() %>% 
  mutate(anova.res = map(data,function(df) aov(value ~ locality, df))) %>% 
  mutate(anova.sum = map(anova.res,function(res) broom::tidy(res))) %>% 
  mutate(lm.div = map(data,function(df) lm(value ~ alpha.div+locality, df))) %>% 
  mutate(lm.raw = map(lm.div,function(res) summary(res))) %>% 
  mutate(lm.sum = map(lm.div,function(res) broom::tidy(res))) %>% 
  unnest(anova.sum) %>% 
  unnest(lm.sum,names_repair = "universal") 

## anova results plot
nest.data %>% 
  select(trait,data,p.value...9,term...4) %>% 
  filter(term...4 == "locality") %>% 
  unnest(data) %>% 
  mutate(fill_col = ifelse(p.value...9 < 0.05,"sig","not_sig")) %>% 
  ggplot(.,aes(locality,value,fill = fill_col))+
  geom_violin()+
  geom_boxplot(width = 0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#D65108","#0075C4"))+
  theme_minimal()+
  facet_wrap(.~trait,scales = "free_y",ncol = 3)

#only significant
nest.data %>% 
  select(trait,data,p.value...9,term...4) %>% 
  filter(term...4 == "locality") %>% 
  unnest(data) %>% 
  mutate(fill_col = ifelse(p.value...9 < 0.05,"sig","not_sig")) %>% 
  ggplot(.,aes(locality,value,fill = fill_col))+
  geom_violin()+
  geom_boxplot(width = 0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#D65108","#0075C4"))+
  theme_minimal()+
  facet_wrap(.~trait,scales = "free_y",ncol = 3)
unique(nest.data$lm.raw)

data.gg.lm<- nest.data %>% 
  select(trait,data,p.value...15,term...11,estimate) %>% 
  filter(term...11 == "alpha.div") %>% 
  unnest(data) 

data.gg.lm %>% 
  filter(trait == "open") %>% 
  ggplot(.,aes(alpha.div,value,color = locality))+
  geom_point()

  mutate(fill_col = ifelse(p.value...15 < 0.05,"sig","not_sig")) %>% 
  ggplot(.,aes(alpha.div,value,fill = fill_col))+
  geom_point()+
  scale_fill_manual(values = c("#D65108","#0075C4"))+
  theme_minimal()+
  facet_wrap(.~trait,scales = "free_y",ncol = 3)
