library(DirichletReg)
library(ggplot2)

prefix.plots <- "C:/Users/delza/Documents/juan/240821_paper_V3/figures/tesis/composition/"

meta <- read.csv("C:/Users/delza/Documents/juan/230328_input_data/meta.csv", sep=";", header=TRUE,stringsAsFactors = FALSE)
meta2 <- meta[meta["type"] == "measured",]
rownames(meta2)<-meta2$ID_pseudo

df_port <- read.table("C:/Users/delza/Documents/juan/230328_input_data/20220803_plate_normalized20220803_data_donwshift_imputed.csv", sep=" ", header=TRUE, stringsAsFactors = FALSE)

to.remove <- colnames(meta[,9:27])[!colnames(meta[,9:27]) %in% colnames(df_port)]
meta <- meta[,which(!colnames(meta) %in% to.remove)]

meta2$D_ASTHMA_SEVERITYGRADE_SCREEN[meta2$D_ASTHMA_SEVERITYGRADE_SCREEN == 9]<-0
meta2[rownames(meta2),colnames(df_port)] <- df_port[rownames(meta2),colnames(df_port)]

meta2$year <- as.factor(meta2$year)
meta2$plate_new <- as.factor(meta2$plate_new)
meta2$D_SEX<-as.factor(meta2$D_SEX)
meta2$D_ASTHMA_SEVERITYGRADE_SCREEN <-as.factor(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN)
meta2$D_SmokingStatus<-as.factor(meta2$D_SmokingStatus)

meta2$ICS_YN <- ifelse(meta2$QMA_YN_ICS_Einzelinhalator == 1|meta2$QMA_YN_SystemicSteroids == 1,1,0) %>% as.factor()
meta2$OCS_YN <- meta2$QMA_YN_SystemicSteroids %>% as.factor()

##################################################
################# COMPOSITION ####################
##################################################

sputum_cell_prop_col<-c("MS_DIFF_MAKROS", "MS_DIFF_NEUT", "MS_DIFF_EOS", "MS_DIFF_LYM",
                        "MS_DIFF_FLIMMEREPITHEL", "MS_DIFF_MONOS")

dirig.meta <- subset(meta2, D_ASTHMA_SEVERITYGRADE_SCREEN!=0)
dirig.meta$Y <- DR_data(mutate_all(dirig.meta[sputum_cell_prop_col], function(x) as.numeric(as.character(x))))

res2<-DirichReg(Y ~ Eotaxin.3 + G.CSF + IFN.g + IL.10+ IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24
                + IL.33 + IL.4 + IL.5+ IL.8 + Periostin + SCGB1A.1 + TNF.alpha +
                  D_SEX+D_AGE + D_SmokingStatus+year + ICS_YN,
                data=dirig.meta, control=list(iterlim = 50000))

summary(res2)
confint(res2,level = 0.95)

#######################################

cyt <- list(Eotaxin.3 = "CCL26",
            G.CSF = "G-CSF",
            IFN.g = "IFN-g",
            IL.10 = "IL-10",
            IL.13 = "IL-13",
            IL.17 = "IL-17",
            IL.1alpha = "IL-1a",
            IL.1F7 = "IL-37",
            IL.24 = "IL-24",
            IL.33 = "IL-33",
            IL.4 = "IL-4",
            IL.5 = "IL-5",
            IL.8 = "IL-8",
            Periostin = "POSTN",
            SCGB1A.1 = "SCGB1A1",
            TNF.alpha = "TNF-a")

cells <- list(MS_DIFF_MAKROS = "Macrophages",
              MS_DIFF_NEUT = "Neutrophils",
              MS_DIFF_EOS = "Eosinophils",
              MS_DIFF_LYM = "Lymphocytes",
              MS_DIFF_FLIMMEREPITHEL = "Ciliated epithelial cells",
              MS_DIFF_MONOS = "Monocytes")


i <- 1

results <- confint(res2,level = 0.95)
sums <- summary(res2)
to.plot <- data.frame()

for (i in seq(6)) {
  to.plot <- rbind(
    to.plot,
    data_frame(
      cell = cells[[i]],
      cyts = unlist(cyt),
      coefs = results$coefficients[[i]][-1][1:16],
      ci.1 = results$ci[[1]][[i]][,1][-1][1:16],
      ci.2 = results$ci[[1]][[i]][,2][-1][1:16]
    )
  )
}

to.plot$p <- NA

for (j in seq(16)) {
  to.plot[which(to.plot$cyts == cyt[[j]]),"p"] <- sums$coef.mat[which(rownames(sums$coef.mat) == names(cyt[j])),4]
}


to.plot <- to.plot[which(to.plot$p < 0.05),]

to.plot.sig <- to.plot[which(to.plot$cyts %in% c("IL-14", "IL-10", "G-CSF", "POSTN", "IL-5", "CCL26","TNF-a")),]
to.plot.non.sig <- to.plot[which(!to.plot$cyts %in% c("IL-14", "IL-10", "G-CSF", "POSTN", "IL-5", "CCL26","TNF-a")),]


p <- to.plot.sig %>%
  ggplot(aes(x = coefs, y = cyts)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci.1, xmax = ci.2), height = 0.05) +
  facet_wrap(~ cell, scales = "fixed") +
  theme_bw() +
  labs(
    x = "Coefficient (95% CI)",
    y = "Cytokine"
  ) + theme(strip.background =element_rect(fill="#4292C6")) + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"sig.pdf"),plot = p, width = 14,height = 8,units = "cm")

p <- to.plot.non.sig %>%
  ggplot(aes(x = coefs, y = cyts)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci.1, xmax = ci.2), height = 0.05) +
  facet_wrap(~ cell, scales = "fixed") +
  theme_bw() +
  labs(
    x = "Coefficient (95% CI)",
    y = "Cytokine"
  ) + theme(strip.background =element_rect(fill="#4292C6")) + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"non_sig.pdf"),plot = p, width = 14,height = 14,units = "cm")





