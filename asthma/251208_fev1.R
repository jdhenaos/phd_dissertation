library(MASS)
library(dplyr)
library(effects)
library(ggplot2)

prefix.plots <- "C:/Users/delza/Documents/juan/240821_paper_V3/figures/tesis/fev1/"

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


##################################
############## MODERATE ##########
##################################

FEV1 <-meta2[!is.na(meta2$L_SPIRO_FEV1perc), ]

fit <- lm(L_SPIRO_FEV1perc ~ Eotaxin.3 + G.CSF + IFN.g + IL.10+ IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5+ IL.8 + Periostin + SCGB1A.1 + TNF.alpha
          + year+D_SEX+D_AGE + D_SmokingStatus + OCS_YN + ICS_YN,
          data = subset(FEV1, D_ASTHMA_SEVERITYGRADE_SCREEN==2))
summary(fit)

res<-summary(fit)$coefficients
sig <- res[res[,"Pr(>|t|)"]<=0.05,]

eff_plot <- Effect(focal.predictors = "TNF.alpha" , fit, xlevels=list(TNF.alpha=seq(min(meta2$TNF.alpha) - 1,max(meta2$TNF.alpha) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = TNF.alpha, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = TNF.alpha, y = L_SPIRO_FEV1perc), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = TNF.alpha, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "TNF-a (Moderate)\n(Beta=-4.68; p=0.03)",
       x = "Normalized abundances",
       y = "FEV1 [%pred.]") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"TNFa_moderate.pdf"),plot = p, width = 7,height = 7,units = "cm")

##################################
############## SEVERE ############
##################################

fit <- lm(L_SPIRO_FEV1perc ~ Eotaxin.3 + G.CSF + IFN.g + IL.10+ IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5+ IL.8 + Periostin + SCGB1A.1 + TNF.alpha
          + year+D_SEX+D_AGE + D_SmokingStatus + OCS_YN + ICS_YN,
          data = subset(FEV1, D_ASTHMA_SEVERITYGRADE_SCREEN==3))
summary(fit)

res<-summary(fit)$coefficients
sig <- res[res[,"Pr(>|t|)"]<=0.05,]

eff_plot <- Effect(focal.predictors = "IL.10" , fit, xlevels=list(IL.10=seq(min(meta2$IL.10) - 1,max(meta2$IL.10) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.10, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.10, y = L_SPIRO_FEV1perc), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.10, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-10 (Severe)\n(Beta=4.52; p=0.02)",
       x = "Normalized abundances",
       y = "FEV1 [%pred.]") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL10_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")







































