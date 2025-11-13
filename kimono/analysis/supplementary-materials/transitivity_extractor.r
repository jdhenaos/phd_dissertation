library(tidyverse)
library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)
library(viridis)

#general setup
dir <- "/workspaces/phd_dissertation/kimono/analysis/data/derived_data/"
file_dir <- "/workspaces/phd_dissertation/kimono/analysis/figures/BRCA/"
#file_dir <- "/workspaces/phd_dissertation/kimono/analysis/figures/MIBC/"

#load library
setwd(dir)

dat_raw <- fread("220411_all_project_info.csv")
#dat_raw <- fread("221116_all_project_info.csv")

#dat_raw <- cbind(dat_raw[,c(1,2,3,16,18)],nodes=dat_raw$V6+dat_raw$V8+dat_raw$V10)
#dat_raw <- dat_raw[dat_raw$iteration == 1,]

dat_raw$method <- gsub(pattern = "_moran|_mos|_mor|_so|_sor",replacement = "",x = dat_raw$method)

dat <- dat_raw[tolower(dat_raw$method) %in% c("galasso_false","galasso_true","knn_kimono","salasso_false","salasso_true","cocolasso","hmlasso"),]


#ggplot(dat, aes(x=clustering, y=n_genes+n_methylation+n_cnv, color=as.factor(method)) )+
#  geom_point(alpha=0.5) + facet_grid(experiment~noise)

#ggplot(dat, aes(y=clustering, x=missingness, color=as.factor(method)) )+
#  geom_point(alpha=0.5) + facet_grid(experiment~noise)

#ggplot(dat, aes(x=missingness, y=clustering,fill=as.factor(method)))+
#  scale_y_continuous(trans='sqrt') + geom_boxplot() + facet_grid(experiment~noise)


dat$method <- as.factor(dat$method) 
dat <- dat %>% mutate(method = fct_relevel(method,  "knn_kimono",    "galasso_false" ,"galasso_true" ,  "salasso_false" , "salasso_true" ,"cocolasso",  "hmlasso"  ))
dat$experiment <- as.factor(dat$experiment)
dat <- dat %>% mutate(experiment = fct_relevel(experiment,  "so",    "moran" ,"mos" ,  "sor" , "mor"))

palette_OkabeIto <- c(  "#CC79A7", "#D55E00",  "#E69F00","#0072B2", "#56B4E9", "#009E73", "#F0E442")

dat <- aggregate(dat[, -c(1:8)], list('noise'=dat$noise,
                                      'missingness'=dat$missingness,
                                      'method'=dat$method,
                                      'experiment'=dat$experiment), function(x){mean(x,na.rm=TRUE)})

dat_raw$method <- as.factor(dat_raw$method) 
dat_raw <- dat_raw %>% mutate(method = fct_relevel(method,  "knn_kimono",    "galasso_false" ,"galasso_true" ,  "salasso_false" , "salasso_true" ,"cocolasso",  "hmlasso"  ))
dat_raw$experiment <- as.factor(dat_raw$experiment)
dat_raw <- dat_raw %>% mutate(experiment = fct_relevel(experiment,  "so",    "moran" ,"mos" ,  "sor" , "mor"))

dat_raw <- dat_raw[ dat_raw$noise== 0 & dat_raw$experiment %in% 'so' & !dat_raw$method %in% "bdcoco_lasso",]
dat <- dat[ dat$noise== 0 & dat$experiment %in% 'so' & !dat$method %in% "bdcoco_lasso",]

dat_raw$adaptive <- grepl(pattern = "true",x = dat_raw$method)
dat$adaptive <- grepl(pattern = "true",x = dat$method)

################################################

method <- c("knn_kimono",    "galasso_false" ,"galasso_true" ,  "salasso_false" , "salasso_true" ,"cocolasso",  "hmlasso")
missingness <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
variable <- c("transitivity")

printer <- function(which_method,which_missingness,which_variable){

    sub_dat <- dat[which(dat$method == which_method & dat$missingness == which_missingness),]
    sub_raw_dat <- dat_raw[which(dat_raw$method == which_method & dat_raw$missingness == which_missingness),]

    m <- sub_dat$clustering
    d <- sd(sub_raw_dat$clustering,na.rm=TRUE)
    
    return(paste0(round(m,3),"+-",round(d,2)))
}

for(j in variable){
    cat(j,"\n")
    for(m in missingness){
        cat("missingness: ",m,"\n") 
        for(n in method){
            label <- printer(which_method=n,which_missingness=m,which_variable=j)
            cat(n,": ",label,"\n")
        }
    }
    cat("\n")
}



