library(MASS)
library(dplyr)
library(effects)
library(ggplot2)

prefix.plots <- "C:/Users/delza/Documents/juan/240821_paper_V3/figures/tesis/feno/"

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
############## MILD ##############
##################################

fit.overdisp<- glm.nb( L_FeNO_Value ~ Eotaxin.3 + G.CSF + IFN.g + IL.10+ IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5+ IL.8 + Periostin + SCGB1A.1 + TNF.alpha
                       + year+D_SEX+D_AGE + D_SmokingStatus + ICS_YN, data=subset(meta2, D_ASTHMA_SEVERITYGRADE_SCREEN==1), link = log)
summary(fit.overdisp )

res<-summary(fit.overdisp)$coefficients
sig <- res[res[,"Pr(>|z|)"]<=0.05,][-1,]

eff_plot <- Effect(focal.predictors = rownames(sig)[1] , fit.overdisp, xlevels=list(IL.1F7=seq(min(meta2$IL.1F7) - 1,max(meta2$IL.1F7) + 1,0.1)))
eff_data <- as.data.frame(eff_plot)

p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.1F7, y = fit), color = "blue") +
  geom_point(data = meta2, aes(x = IL.1F7, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.1F7, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-37 (Mild)\n(Beta=0.27; p=0.05)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL37_mild.pdf"),plot = p, width = 7,height = 7,units = "cm")

eff_plot <- Effect(focal.predictors = rownames(sig)[2] , fit.overdisp, xlevels=list(IL.33=seq(min(meta2$IL.33) - 1,max(meta2$IL.33) + 1,0.1)))
eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.33, y = fit), color = "blue") +
  geom_point(data = meta2, aes(x = IL.33, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.33, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-33 (Mild)\n(Beta=-0.11; p=0.01)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL33_mild.pdf"),plot = p, width = 7,height = 7,units = "cm")

##################################
############## MODERATE ##########
##################################

fit.overdisp<- glm.nb( L_FeNO_Value ~ Eotaxin.3 + G.CSF + IFN.g + IL.10+ IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5+ IL.8 + Periostin + SCGB1A.1 + TNF.alpha
                       + year+D_SEX+D_AGE + D_SmokingStatus + OCS_YN + ICS_YN, data=subset(meta2, D_ASTHMA_SEVERITYGRADE_SCREEN==2), link = log)
summary(fit.overdisp )

res<-summary(fit.overdisp)$coefficients
sig <- res[res[,"Pr(>|z|)"]<=0.05,][-1,]

eff_plot <- Effect(focal.predictors = "IL.5" , fit.overdisp, xlevels=list(IL.5=seq(min(meta2$IL.5) - 1,max(meta2$IL.5) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.5, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.5, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.5, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-5 (Moderate)\n(Beta=0.14; p=0.003)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL5_moderate.pdf"),plot = p, width = 7,height = 7,units = "cm")

##################################
############## SEVERE ############
##################################

fit.overdisp<- glm.nb( L_FeNO_Value ~ Eotaxin.3 + G.CSF + IFN.g + IL.10+ IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5+ IL.8 + Periostin + SCGB1A.1 + TNF.alpha + year+D_SEX+D_AGE + D_SmokingStatus + OCS_YN + ICS_YN, data=subset(meta2, D_ASTHMA_SEVERITYGRADE_SCREEN==3), link = log)

res<-summary(fit.overdisp)$coefficients
sig <- res[res[,"Pr(>|z|)"]<=0.05,][-1,]

eff_plot <- Effect(focal.predictors = "IFN.g" , fit.overdisp, xlevels=list(IFN.g=seq(min(meta2$IFN.g) - 1,max(meta2$IFN.g) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IFN.g, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IFN.g, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IFN.g, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IFN-g (Severe)\n(Beta=-0.14; p=0.02)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IFNg_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")

####################################

eff_plot <- Effect(focal.predictors = "IL.13" , fit.overdisp, xlevels=list(IL.13=seq(min(meta2$IL.13) - 1,max(meta2$IL.13) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.13, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.13, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.13, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-13 (Severe)\n(Beta=-0.13; p=0.01)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL13_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")

######################################

eff_plot <- Effect(focal.predictors = "IL.17" , fit.overdisp, xlevels=list(IL.17=seq(min(meta2$IL.17) - 1,max(meta2$IL.17) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.17, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.17, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.17, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-17 (Severe)\n(Beta=-0.17; p=0.02)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL17_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")

########################################

eff_plot <- Effect(focal.predictors = "IL.1F7" , fit.overdisp, xlevels=list(IL.1F7=seq(min(meta2$IL.1F7) - 1,max(meta2$IL.1F7) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.1F7, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.1F7, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.1F7, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-37 (Severe)\n(Beta=0.16; p=0.04)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL37_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")

##########################################

eff_plot <- Effect(focal.predictors = "IL.24" , fit.overdisp, xlevels=list(IL.24=seq(min(meta2$IL.24) - 1,max(meta2$IL.24) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.24, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.24, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.24, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-24 (Severe)\n(Beta=0.06; p=0.04)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL24_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")

###########################################

eff_plot <- Effect(focal.predictors = "IL.4" , fit.overdisp, xlevels=list(IL.4=seq(min(meta2$IL.4) - 1,max(meta2$IL.4) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.4, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.4, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.4, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-4 (Severe)\n(Beta=0.15; p=0.04)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL4_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")

################################################

eff_plot <- Effect(focal.predictors = "IL.5" , fit.overdisp, xlevels=list(IL.5=seq(min(meta2$IL.5) - 1,max(meta2$IL.5) + 1,0.1)))

eff_data <- as.data.frame(eff_plot)

# Create a ggplot with effects plot and data points
p <- ggplot() +
  geom_line(data = eff_data, aes(x = IL.5, y = fit), color = "blue") + 
  geom_point(data = meta2, aes(x = IL.5, y = L_FeNO_Value), color = "black", shape = 16) +
  geom_ribbon(data = eff_data, aes(x = IL.5, ymin = lower, ymax = upper), fill = "blue", alpha = 0.08) +
  labs(title = "IL-5 (Severe)\n(Beta=0.13; p=0.01)",
       x = "Normalized abundances",
       y = "FeNO (ppb)") +
  theme_bw() + theme(text = element_text(size = 8))

ggsave(filename = paste0(prefix.plots,"IL5_severe.pdf"),plot = p, width = 7,height = 7,units = "cm")
































































































