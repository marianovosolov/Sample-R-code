source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
library("ShortRead")
library("R.utils")
setwd("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/")
#read the data
srx35_unzip<- R.utils::gunzip("SRX350671/SRX350671.fastq.gz")
srx35<- readLines(srx35_unzip)
srx35_no_space <- readFastq(srx35[nchar(srx35)!=0])
which(nchar(srx35_no_space)==0)
head(srx35_no_space)
class(srx35_no_space)
write(srx35_no_space,file = "SRX350671/SRX350671_new.fastq")
srx35_fq<- readFastq("SRX350671/SRX350671_new.fastq")
srx11<- readFastq("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/SRX1607911/SRX1607911.fastq")
class(srx11)
# srx12<- readFastq("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/SRX1607912/SRX1607912.fastq")
# srx13<- readFastq("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/SRX1607913/SRX1607913.fastq")
# srx14<- readFastq("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/SRX1607914/SRX1607914.fastq")
# srx15<- readFastq("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/SRX1607915/SRX1607915.fastq")
# srx16<- readFastq("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/SRX1607916/SRX1607916.fastq")
#save all the files into a list
file.list<- list(srx11,srx12,srx13,srx14,srx15,srx16)

names(file.list)<- c("srx11","srx12","srx13","srx14","srx15","srx16")
try.data<- srx12[1:6]
try.data@id
a<- gsub(".*H","aaaH",try.data@id)
a
a<- gsub(".*.1z","",a)
a
a<- gsub("z.*","",a)
a
#change the problematic part of the read name
#run a loop over all the files and change the names of the reads
i=2
for(i in 1:length(file.list)){
  file.list[[i]]@id<- BStringSet(gsub(".*H","aaaH",file.list[[i]]@id))
  file.list[[i]]@id<- BStringSet(gsub(" l","zzz",file.list[[i]]@id))
  file.list[[i]]@id<- BStringSet(gsub(".*aaa","",file.list[[i]]@id))
  file.list[[i]]@id<- BStringSet(gsub("zzz.*","",file.list[[i]]@id))
  print(head(file.list[[i]]@id))
}



file.list[[2]]@id<- BStringSet(paste0("H",file.list[[2]]@id))
file.list[[3]]@id<- BStringSet(paste(file.list[[3]]@id,"/1",sep=""))
head(file.list[[3]]@id)
writeFastq(object = file.list[[3]],file="SRX1607913/SRX1607913_1.fastq",mode="w",compress=F)
#save the files back to one big file

for (fl in file.list) {
  fq = fl
  writeFastq(fq, "All_SRX16079/all_read_OD_norway_2018.fastq", mode="a", compress=F)
}
f=1
#also save each file indipendently
for (f in 1:length(file.list)){
  path.to.file<- paste("D:/Maria Novosolov/Post-doc/Molecular Data/Reads/Bergen_Norway/","SRX160791",f,"/","SRX160791",f,"_1",".fastq",sep="")
  fq.file<- file.list[[f]]
  writeFastq(object = fq.file,file=path.to.file,mode="w",compress=F)
}


