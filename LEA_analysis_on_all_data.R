### This script is using the LEA package to analyse the data in a vcf format calculate the admixture proportions, and then imputing the data. 


lib_loc = "/science/willerslev/users-shared/science-snm-willerslev-rlw363/Rlibrary"
library(vcfR)
library(LEA,lib.loc = lib_loc)
setwd("~/sharedhome/LEA")
all.vcf<- read.vcfR("all_vcf_merged_PS_noM48.vcf.gz")
# Scanning file to determine attributes.
# File attributes:
#   meta lines: 209410
# header_line: 209411
# variant count: 20268640
# column count: 102
# Meta line 209410 read in.
# All meta lines processed.
# gt matrix initialized.
# Character matrix gt created.
# Character matrix gt rows: 20268640
# Character matrix gt cols: 102
# skip: 0
# nrows: 20268640
# row_num: 0
# Processed variant: 20268640
# All variants processed
gt <- extract.gt(all.vcf, element = 'GT', as.numeric = TRUE)
saveRDS(gt,"vcf_matrix_noM48.rds")
readRDS("vcf_matrix_noM48.rds")
#the matrix is in the wrong dimationality so I need to transpose it

#check how many are 3 or bigger
which(gt>=4, arr.ind=TRUE)
# no value of 4 or larger
gt2<- gt
#change all the NA to 9
gt2[is.na(gt2)]<-9
#change all the cells with 3 to NA
gt2[gt2==3]<- NA
#check the length of all the NA's
length(which(is.na(gt2)))
#only 3876 cells had 3 in them
#remove all the lines that have NA's
gt3<- gt2[which(complete.cases(gt2)==T),]
#20267012 sites remaned after removing the sites with 3 alleles
# Transpose the matrix
gt.t<- t(gt3)
#save to rds the new matrix and the transposed
saveRDS(gt3,"no_NA_sno_PS_matrix_noM48.rds")
saveRDS(gt.t,"trans_no_NA_sno_PS_matrix_noM48.rds")
#write a lfmm object
write.lfmm(gt.t,output.file = "PS_noM48_snp_M2.lfmm")
# 1628 sites were excluded because they had more than two ALT alleles
lfmm2geno(input.file = "PS_noM48_snp_M2.lfmm",output.file = "PS_noM48_snp_M2.geno")

#start an empty project
project.PS = NULL
# calculate the admixture clusters
project.PS = snmf("PS_noM48_snp_M2.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new",
               CPU = 44)
# plot cross-entropy criterion for all runs in the snmf project
ent.plot<- plot(project.PS, col = "blue", pch = 19, cex = 1.2)
ent.plot
pdf("nuc_plots/ent_plot.pdf")
ent.plot
dev.off()

#Best K is 2
# select the best run for K = 4
best = which.min(cross.entropy(project.PS, K = 2))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
bp<- barchart(project.PS, K = 2, run = best,border = NA, space = 0,col = my.colors,xlab = "Individuals", ylab = "Ancestry proportions",main = "Ancestry matrix")
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .3)


# Population differentiation tests
p = snmf.pvalues(project.PS,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 2)

#Plotting the p-values
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)



#### Imputate the data ######
project.missing = snmf("PS_noM48_snp_M2.geno", K = 2, entropy = TRUE, repetitions = 10,project = "new",CPU = 44)

## If the R was stoped for some reason I can reload the project with 
project = load.snmfProject("project_name.snmfProject")
#calculate the best run with the lowest cross-entropy
best = which.min(cross.entropy(project.missing, K = 2))

# Impute the missing genotypes
impute(project.PS, "PS_noM48_snp_M2_imputed.lfmm",method = 'mode', K = 2, run = best)

