library(corrplot)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(reshape2)
library(limma)

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

meta2[, colnames(df_port)] <- removeBatchEffect(t(meta2[, colnames(df_port)]),batch = meta2$year) %>% t()

to.cor.plot <- meta2[,which(colnames(meta2) %in% colnames(df_port))] %>%
  cbind(meta2$MS_DIFF_EOS,
        meta2$MS_DIFF_NEUT,
        meta2$MB_DIFF_EOS_abs_corrected.FZB,
        meta2$MB_DIFF_NEUT_abs
        ) %>% na.omit()

colnames(to.cor.plot) <- c("CCL-26","G-CSF","IFN-g","IL-10",
                           "IL-13","IL-17","IL-1a","IL-37",
                           "IL-24","IL-33","IL-4","IL-5",
                           "IL-8","POSTN","SCGB1A1","TNF-a",
                           "Sputum Eosinophils","Sputum Neutrophils","Blood Eosinophils","Blood Neutrophils")


corr_matrix <- cor(to.cor.plot)
sparse_corr <- corr_matrix
sparse_corr[abs(sparse_corr) < 0.8] <- 0
sparse_corr <- sparse_corr[which(rowSums(sparse_corr) != 1.0),colSums(sparse_corr) != 1.0]

pdf(file = "results/kruskal/correlation.pdf",width = 7/2.54,height = 7/2.54)
corrplot(sparse_corr,method="num",number.cex = 0.6,cl.cex=0.6,tl.cex=0.6,tl.col = 'black',type = 'upper',col=colorRampPalette(c("#313695", "white", "#A50026"))(200))
dev.off()

