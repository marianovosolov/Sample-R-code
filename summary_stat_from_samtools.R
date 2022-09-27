#### This script takes the output from samtools for all the samples and puts the needed values into tables

#load the library
library(tidyverse)
##### Reads mapped and duplicates summary ##### 
#create a data frame that will hold the summary of the reads number and the number of duplicates
read_summary<- data.frame(sample_name=rep(NA,94),real_sample=rep(NA,94),island=rep(NA,94),reads_mapped=rep(NA,94),total_reads=rep(NA,94),perc_mapped=rep(NA,94),total_dup=rep(NA,94),perc_dup=rep(NA,94),data_type=rep(NA,94))
seq_data<- as.data.frame(read_csv("~/sharedhome/L.Dalen_19_09_sample_info.csv"))
#list all the files with txt
all_files_flagstat<- list.files(path = "~/sharedhome/samtools_summary_merged",pattern = "*\\.txt",full.names = T)
#get all the basenames of the filess
all_files_names<- list.files(path = "~/sharedhome/samtools_summary_merged",pattern = "*\\.txt",full.names = F)

#run a loop that gets the name of the sample (=the first part in name of the file), the original tag (KXX or MXX), the island it is from, the amount of reads mapped and the percentage of duplicate reads in the mapping
for(i in 1:length(all_files_flagstat)){
  read_summary$sample_name[i]<- word(all_files_names[i],1,sep="\\.")
  a<- seq_data[which(seq_data[1]==paste(word(all_files_names[i],1,sep="\\_"),word(all_files_names[i],2,sep="\\_"),sep="_")),2]
  read_summary$real_sample[i]<- a
  read_summary$island[i]<- str_sub(a,1,1)
  flagstat<- strsplit(read_file(all_files_flagstat[i]),split = "\n")
  read_summary$reads_mapped[i]<- as.numeric(word(flagstat[[1]][5],1))
  read_summary$total_reads[i]<- as.numeric(word(flagstat[[1]][1],1))
  read_summary$perc_mapped[i]<- (read_summary$reads_mapped[i]/read_summary$total_reads[i])*100
  read_summary$total_dup[i]<- as.numeric(word(flagstat[[1]][4],1))
  read_summary$perc_dup[i]<- (as.numeric(word(flagstat[[1]][4],1))/read_summary$reads_mapped[i])*100
}
#save the table to RDS
write_csv(read_summary,"~/sharedhome/reads_mapped_all_samples130120.csv")

#### Reads coverage for each sample #####
#list all the files with csv
all_files_coverage<- list.files(path = "~/sharedhome/sharedhome/samtools_summary_merged",pattern = "*\\.csv",full.names = T)
#get all the basenames of the filess
all_files_names_cov<- list.files(path = "~/sharedhome/sharedhome/samtools_summary_merged",pattern = "*\\.csv",full.names = F)

#create the data frame that will hold the data
coverage_summary<- data.frame(sample_name=rep(NA,94),real_sample=rep(NA,94),island=rep(NA,94),mean_coverage=rep(NA,94),data_type=rep(NA,94))

for(j in 1:length(all_files_coverage)){
  coverage_summary$sample_name[j]<- word(all_files_names_cov[j],1,sep="\\.")
  a<- seq_data[which(seq_data[1]==paste(word(all_files_names_cov[j],1,sep="\\_"),word(all_files_names_cov[j],2,sep="\\_"),sep="_")),2]
  coverage_summary$real_sample[j]<- a
  coverage_summary$island[j]<- str_sub(a,1,1)
  cov_data<- read_csv(all_files_coverage[j],col_names = F) %>% 
    separate(col = X1,into = c("loc_total","read_num"),sep = " ") %>% 
    type_convert()
  coverage_summary$mean_coverage[j]<- weighted.mean(cov_data$read_num,cov_data$loc_total)
  coverage_summary$data_type[j]<- str_extract(all_files_names_cov[j],"[A-z]+.(?=_)")
}


write_csv(coverage_summary,"~/sharedhome/coverage_summary_130120.csv")




