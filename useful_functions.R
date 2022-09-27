#### Function to plot results from PCAngsd ######
library(tidyverse)

plot_PCAngsd<- function(data_file,var_data,var_names,plot_cols = NULL,pc_to_plot=c("PC1","PC2")){
  data<- PCAngsd_dataLoad(data_path = data_file)
  data_df<- as.data.frame(data$vectors[,1:5])
  data_df<- cbind(var_data[,var_names[1]],var_data[,var_names[2]],data_df)
  names(data_df)<- c(var_names[1],var_names[2],"PC1","PC2","PC3","PC4","PC5")
  pca_plot<-ggplot(data_df,aes(get(pc_to_plot[1]),get(pc_to_plot[2]),col=get(var_names[2]),label = get(var_names[1])))+
    geom_point()+
    scale_color_manual(name = var_names[2],values = plot_cols)+
    theme_minimal()+
    theme(legend.position = "bottom")
  return(pca_plot)
  
}

#load the covariance matrix and convet it to PCA values
PCAngsd_dataLoad<- function(data_path){
  cov_mat<- as.matrix(read.table(here::here(data_path)))
  cov_mat_eg<- eigen(cov_mat)
  return(cov_mat_eg)
}

# calculate and print the contribution table for each PC
PCAngsd_prop_table<- function(data_path){
  data<- PCAngsd_dataLoad(data_path = data_path)
  eg_values<- data$values
  knitr::kable(col.names = c("PC1","PC2","PC3","PC4","PC5"),rbind(
    SD = sqrt(eg_values),
    Proportion = eg_values/sum(eg_values),
    Cumulative = cumsum(eg_values)/sum(eg_values))[,1:5])
}



#### Function to plot sfs data ####
library(tidyverse)
nnorm <- function(x) x/sum(x)

p_realSFS<- function(data_path,pop_name,main=NULL,fold = TRUE){
  data1<- rbind(
    pop1=scan(here::here(data_path))[-1]
  )
  data1 <- t(apply(data1,1,nnorm))
  if(fold==TRUE){
    dim_plot = round(length(data1)/2)
  }else{
    dim_plot = round(length(data1))
  }
  barplot(data1[1:dim_plot],beside=T,legend=c(pop_name),names=1:dim_plot,main=main)
}


##### Functions to read and plot neighbor joining tree ####
library(tidyverse)
library(ape)

read_ibsMat<- function(data_path,col.names = col.names){
  data<- readr::read_delim(here::here(data_path),delim = "\t",col_names = F) %>% 
    select(-ncol(.))
  colnames(data)<- col.names
  data_mat<- as.matrix(data)
  return(data_mat)
}

plot_ibsTree<- function(data_file,col.nam,color,pop.data,root=F,cex=NULL){
  data_mat<- read_ibsMat(data_file,col.nam)
  tree<- fastme.ols(data_mat)
  pop.col<- pop.data %>% 
    mutate(col_island = case_when(island=="K"~color[1],TRUE~color[2])) %>% 
    select(user_id,island,col_island)
  
  if(root==T){
    tree<- phytools::midpoint.root(tree)
    pop.col.phy<- pop.col[match(tree$tip.label,pop.col$user_id),]
    plot(tree,tip.color = pop.col.phy$col_island,label.offset = 0.0001,cex=cex)
  }else{
    pop.col.phy<- pop.col[match(tree$tip.label,pop.col$user_id),]
    plot(tree,tip.color = pop.col.phy$col_island,label.offset = 0.0001)
  }
  
}


##### Calculate fst peaks #####
fst_peaks<- function(vec,limits = c()){
  #find the highest Fst value
  max_fst<- max(vec)
  # Compare it to the a critical value
  # c1 <- 0.26
  c1 <- limits[1]
  # c2 <- 0.21
  c2<- limits[2]
  peak_fst<- data.frame(loc = NA,fst_val = NA)
  
  # if it's higher save it to a new list and remove it
  for(i in 1:length(fst_high)){
    max_fst<- max(fst_high,na.rm = T)
    if(max_fst > c1){
      peak_fst[i,1]<- which(fst_high==max_fst)
      peak_fst[i,2] <- max_fst
      loc_r<- which(fst_high==max_fst)
      loc_l<- which(fst_high==max_fst)
      fst_high[which(fst_high==max_fst)]<- NA
      fst_next_r<- max_fst
      fst_next_l<- max_fst
      while(fst_next_r >= c2){
        fst_next_r<- fst_high[loc_r+1]
        if(is.na(fst_next_r)){
          break
        }
        fst_high[which(fst_high==fst_next_r)]<- NA
        loc_r <- loc_r+1
      }
      while(fst_next_l >= c2){
        fst_next_l<- fst_high[loc_l-1]
        if(is.na(fst_next_l)){
          break
        }
        fst_high[which(fst_high==fst_next_l)]<- NA
        loc_l <- loc_l-1
      }
    }else{
      break
    }
  }
  return(peak_fst)
}
