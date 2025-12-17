library(dplyr)
library(DirichletReg)

prefix.plots <- "C:/Users/delza/Documents/juan/240821_paper_V3/figures/"
prefix.tables <- "C:/Users/delza/Documents/juan/240821_paper_V3/tables/"

meta <- read.csv("data/raw_data/meta.csv", sep=";", header=TRUE,stringsAsFactors = FALSE)
meta2 <- meta[meta["type"] == "measured",]
rownames(meta2)<-meta2$ID_pseudo

df_port <- read.table("data/derivative_data/20220803_plate_normalized20220803_data_donwshift_imputed.csv", sep=" ", header=TRUE, stringsAsFactors = FALSE)

to.remove <- colnames(meta[,9:27])[!colnames(meta[,9:27]) %in% colnames(df_port)]
meta <- meta[,which(!colnames(meta) %in% to.remove)]

meta2$D_ASTHMA_SEVERITYGRADE_SCREEN[meta2$D_ASTHMA_SEVERITYGRADE_SCREEN == 9]<-0
meta2[rownames(meta2),colnames(df_port)] <- df_port[rownames(meta2),colnames(df_port)]

meta2$year <- as.factor(meta2$year)
meta2$plate_new <- as.factor(meta2$plate_new)
meta2$D_SEX<-as.factor(meta2$D_SEX)
meta2$D_ASTHMA_SEVERITYGRADE_SCREEN <-as.factor(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN)
meta2$D_SmokingStatus<-as.factor(meta2$D_SmokingStatus)

#########################################
############### COMP ######################
sputum_cell_prop_col<-c("MS_DIFF_MAKROS", "MS_DIFF_NEUT", "MS_DIFF_EOS", "MS_DIFF_LYM",
                    "MS_DIFF_FLIMMEREPITHEL", "MS_DIFF_MONOS")


dirig.meta <- subset(meta2, D_ASTHMA_SEVERITYGRADE_SCREEN!=0)
dirig.meta$Y <- DR_data(mutate_all(dirig.meta[sputum_cell_prop_col], function(x) as.numeric(as.character(x))))

res2<-DirichReg(Y ~ Eotaxin.3 + G.CSF + IFN.g + IL.10+ IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24
                + IL.33 + IL.4 + IL.5+ IL.8 + Periostin + SCGB1A.1 + TNF.alpha +
                D_SEX+D_AGE + D_SmokingStatus+year + ICS_YN,
                data=dirig.meta, control=list(iterlim=5000))

summary(res2)
confint(res2,level = 0.95)
