library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)
library(limma)
library(ggpubr)
library(reshape2)

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
meta2$OCS_mg <- ifelse(is.na(meta2$OCS_mg),0,meta2$OCS_mg)

meta2[, colnames(df_port)] <- removeBatchEffect(t(meta2[, colnames(df_port)]),batch = meta2$year) %>% t()

to.box.plot <- data.frame()
for (i in colnames(df_port)) {
  to.bind <-  data.frame(rep(i,nrow(meta2)),
                         meta2[,i],
                         as.factor(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN),row.names = NULL)
  colnames(to.bind) <- c('Protein','Expression','Condition')
  to.box.plot <- rbind(to.box.plot,to.bind)
}

cyt <- "TNF.alpha"

pair <- pairwise.wilcox.test(meta2[,which(colnames(meta2) == cyt)],meta2$D_ASTHMA_SEVERITYGRADE_SCREEN,p.adjust.method = "BH")
pair <- melt(pair$p.value)
pair <- pair[which(!is.na(pair$value)),]

q <- ggplot(data = subset(to.box.plot, (Protein == cyt)), aes(x = Condition, y = Expression)) + 
    geom_boxplot(aes(fill = Condition), width = 0.4,outlier.shape = NA) + theme_classic() +
      geom_jitter(size=0.7) +
        scale_fill_manual(values = c("#ffeda0","#feb24c","#fc4e2a","#bd0026"),
      labels=c("Healthy", "Mild", "Moderate", "Severe"),name="") +
        scale_x_discrete(labels=c("Healthy", "Mild", "Moderate", "Severe")) +
        xlab("Asthma severity") + ylab("Normalized abundances") +
        labs(title = "TNF-a - Kruskal-Wallis (p.adj= 0.039)",
      subtitle = "Wilcoxon rank-sum (post-hoc) - p.adj:") +
  theme(legend.position = "none", text = element_text(size=8), plot.title = element_text(size=8)) +
  stat_compare_means(method = "kruskal",label.y=5) +
  stat_compare_means(method = "wilcox",comparisons = list(
    c("1","3"),
    c("2","3")
  ))

boxp2 <- ggplot_build(q)
boxp2$data[[4]]$annotation <- formatC(as.vector(sapply(pair[pair$value < 0.05,"value"], function(x){rep(x,3)})),digits = 3)
boxp2$data[[3]]$label <- ""
boxp2$data[[4]]$textsize <- boxp2$data[[3]]$size <- 2.5

pdf(file = "results/kruskal/tnfa_kruskal.pdf",width = 7/2.54,height = 7/2.54)
plot(ggplot_gtable(boxp2))
dev.off()