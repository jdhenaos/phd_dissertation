library(tidyverse)
library(data.table)
library(reshape)
library(ggforestplot)

## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

#general setup
dir <- "/workspaces/phd_dissertation/kimono/analysis/data/raw_data/BRCA/breast_cancer/networks/"

#load library
setwd(dir)

#get all files for one missingness and one noise lvl
get_files_missingness <- function(files_dir,missingness,noise){
  files <- list.files(files_dir)
  files[files %like% paste0("noise_",noise,"_missing_",missingness,"_")]
}

#load -> scale -> extract betweeness
get_scaled_betweeness_DT <- function(file_iteration){
  dat <- fread(file_iteration)
  #scale betweeness by (N-1)(N-2) |N = number of nodes
  N <- length(dat$V1)
  dat$Betweenness <- dat$Betweenness / ((N-1)*(N-2))
  dat <- dat[,c(2,4)] 
  dat
}

#load -> scale -> extract betweeness
get_scaled_betweeness_DT <- function(file_iteration){
  dat <- fread(file_iteration)
  #scale betweeness by (N-1)(N-2) |N = number of nodes
  N <- length(dat$V1)
  dat$Betweenness <- dat$Betweenness / ((N-1)*(N-2))
  dat <- dat[,c(2,4)] 
  dat
}

#combine multiple seed interations
get_betweeness_missingness <- function(files_dir,iterations){
  sub_result <- c()
  sub_result <- paste(files_dir,iterations[1],sep='/') %>% get_scaled_betweeness_DT()
  for (i in 2:length(iterations)) {
    tmp_DT <- paste(files_dir,iterations[i],sep='/') %>% get_scaled_betweeness_DT()
    colnames(tmp_DT) <- c('Node',i)
    sub_result <- merge(sub_result,tmp_DT, by = 'Node',all = TRUE)
  }
  data.table(sub_result[,1] , rowMeans(sub_result[,-1]))
}

#############################################3333

#benchmarks
noise_lvls <- c(0)
missingness_lvls <- c(0,0.1,0.2,0.3,0.4,0.5)

#method

#experiment_list <- c("multi_omics_random","multi_omics_sample","single_omics")
experiment_list <- c("moran","mos","so")
method_list <- c("galasso_false","galasso_true","knn_kimono","salasso_false","salasso_true","cocolasso","hmlasso")

overall <- c()

for (experiment in experiment_list) {
  for (method in method_list) {
    files_dir <- paste(dir,method,experiment,sep="/")
    tryCatch(
      expr = { 
        for (noise in noise_lvls) {
          
          result <- c()
          tmp_missingness <- get_files_missingness(files_dir,missingness_lvls[1],noise)
          result <- get_betweeness_missingness(files_dir,tmp_missingness)
          colnames(result) <- c('Node',missingness_lvls[1])
          
          for (i in 2:length(missingness_lvls)) {
            tryCatch(
              expr = {
                
                tmp_missingness <- get_files_missingness(files_dir,missingness_lvls[i],noise)
                tmp_DT <- get_betweeness_missingness(files_dir,tmp_missingness)
                colnames(tmp_DT) <- c('Node',missingness_lvls[i])
                result <- merge(result,tmp_DT, by = 'Node',all = TRUE)
                
              },
              error = function(e){
                
              },
              warning = function(w){
                
              },
              finally = {
                
              }
            )    
          }
          
          
          r <- result[order(result[,2],decreasing = TRUE),][1:200]
          r[is.na(r)] <- 0
          
          tmp_dat <- data.frame(melt(r), method, noise,experiment)
          
          overall <- rbind(tmp_dat,overall)
          
        }
      },
      error = function(e){
        
      },
      warning = function(w){
        
      },
      finally = {
        
      }
    )   
  }
}

##plotting
overall$method %>% unique
exclude <- overall$noise %in% c(0)   

gg_dat <- overall[ exclude,]
#ggplot(gg_dat, aes(x=variable, y=value,fill=method))+
#  scale_y_continuous(trans='sqrt') + geom_boxplot() + facet_grid(experiment~noise)

gg_dat$method <- as.factor(gg_dat$method)
gg_dat <- gg_dat %>% mutate(method = fct_relevel(method,  "knn_kimono"  ,   "galasso_false", "galasso_true" ,"salasso_false", "salasso_true" , "cocolasso" ,  "hmlasso" ))
gg_dat$stripe <- as.numeric(gg_dat$method)

my_grey <- t_col(color = "grey",percent = 90)

palette_OkabeIto <- c(  "#CC79A7", "#D55E00",  "#E69F00","#0072B2", "#56B4E9", "#009E73", "#F0E442")

levels(gg_dat$method) <- c( "knn_KiMONo" , "GALasso" ,"GALasso\n (adaptive weights)", "SALasso", "SALasso \n(adaptive weights)", "CoCoLasso",  "HMLasso")

gg_dat <- gg_dat[gg_dat$experiment=='so',]