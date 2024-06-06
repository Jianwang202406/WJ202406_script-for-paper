rm(list = ls())
library(ggplot2)
library(RColorBrewer)
############ α diversity
data<- read.csv(file = 'MF_percent.csv',header=TRUE)
colors<-c("darkolivegreen3","gold","dodgerblue4","darkseagreen",
          
          "chartreuse4","darkorange","burlywood2","brown3","#984EA3","cyan3","grey50")

data$sample <- factor(data$sample, levels=c("MH","DXAL","WY", "CBS", "DLS","BTM","JGS","WYS","DHS","JFL","XSBN"), ordered=TRUE)
data$MF <- factor(data$MF, levels=c("Lignin","Tannin","Aromatic", "Protein", "Lipid","Carbohydrate","Others","UHCs"), ordered=TRUE)



p1<-ggplot(data, aes( x = sample, y = percent, fill = MF))+
  
  geom_bar(position = "fill", stat = "identity")+
  
  theme_bw()+
  
  scale_fill_manual(values=colors)+ 
  
  scale_y_continuous(expand = c(0,0))+
  
  labs(x="",y="Relative Abundance (%)",fill="")+
  
  theme(text=element_text(size = 15),
        
        axis.text.y=element_text(size = 15,color = "black"),
        
        axis.text.x=element_text( size = 15,color = "black",angle = 45, hjust =1, vjust = 1),
        
        legend.title=element_text( size = 15), 
        
        legend.text=element_text( size = 12))+ 
  
  theme(panel.grid = element_blank(),
        
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  
  guides(fill=guide_legend(keywidth = 0.8, keyheight = 1.5))

p1
###################################    NMDS

library(vegan)
DOM <- read.delim(file = "otu_like_fill_intensity_NMDS.txt",header = T, row.names = 1)
DOM <- data.frame(t(DOM))
hel <- decostand(DOM, method = 'hellinger')

group <- read.delim('design.txt', stringsAsFactors = FALSE)



bray_dis <- vegdist(hel , method = 'bray')


set.seed(123)

nmds_dis <- metaMDS(bray_dis, k = 2)

stress <-nmds_dis$stress

nmds_dis_site <- data.frame(nmds_dis$points)



stress_text <- paste("Stress  =", round(stress, 4))
adonis_text <- paste(paste("Adonis  =", round(adonis$R2, 2)), "**")[1]
anosim_text <- paste(paste("Anosim  =", round(anosim$statistic, 2)), "**")

library(ggplot2)
library(ggrepel)




nmds_dis_site$ID<- rownames(nmds_dis_site)


merged<-merge(nmds_dis_site,group,by="ID",all.x=TRUE)
library(ggplot2)
library(ggrepel)
library(ggthemes)

library(ggalt)
col1 <- colorRampPalette(c( "#FF8C00","#B0E2FF", "#63B8FF",  "#FFAEB9","#FF3030" ))




p2 <- ggplot(data = merged, aes(x=MDS1, y=MDS2,color=MAT)) +
  geom_point(size = 3.5,alpha=0.9)+

  scale_color_gradientn(name="",colours=col1(21),
                        limits=c(-5,25),
                        
                        guide = guide_legend(key.width=2,key.height=1))+
  geom_hline(yintercept=0,linetype=3,size=1)+ 
  geom_vline(xintercept=0,linetype=3,size=1)+

  geom_encircle(aes(group = group,fill=group),expand=0,spread=0.5,s_shape=0.9,alpha = 0.2,linetype="dashed",size=2,colour="black")+
  geom_text_repel(data = merged, aes(label = sample), color = 'black', size = 2.5,box.padding = unit(0.3, 'lines'), show.legend = FALSE)+     
  scale_fill_manual(values=c( "#FF8C00","#B0E2FF", "#FFAEB9"))+
  labs(x = 'NMDS1', y = 'NMDS2') +
  scale_y_continuous(limit = c(-0.16, 0.12),breaks=c(-0.15,-0.1,-0.05,0.0,0.05,0.1))+
  annotate('text', label = paste('Stress =', round(nmds_dis$stress, 4)), x = 0.137, y = -0.16, size = 5, colour = 'black')+ 
  annotate('text', label = paste(paste("Adonis  =", round(adonis$R2, 2)), "**")[1], x = 0.14, y = -0.145, size = 5, colour = 'black')+
  annotate('text', label = paste(paste("Anosim  =", round(anosim$statistic, 2)), "**"), x = 0.14, y = -0.13, size = 5, colour = 'black')+
  theme_bw()+
  theme(axis.title = element_text(  size = 15,colour = "black"))+
  theme(axis.text = element_text(  size =15,color="black"),axis.text.x=element_text( hjust=1, vjust=1, size = 15,color="black"))+
  
  theme(legend.text =  element_blank(),legend.position = "none",legend.title = element_text( size =15,color="black"))+
  

  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))
p2

#################################shannon
dom<-read.csv(file = "heatmap_dom _all1.csv",header = T)

dom$sample <- factor(dom$sample, levels=c("MH","DXAL","WY", "CBS", "DLS","BTM","JGS","WYS","DHS","JFL","XSBN"), ordered=TRUE)
cols <-colorRampPalette(c("#8B8989", "#FF8C00", "#008B45", "#5D478B", "#436EEE", "#CDAD00"))


p3<-ggplot(dom,mapping = aes(x =sample, y =Shannon_dom,fill=sample,color=sample))+
  
  geom_boxplot(linetype="solid",size=1, width=0.5,alpha=0.6,outlier.shape = NA)+
  scale_fill_manual(name=" ",values = cols(11))+ 
  scale_color_manual(name=" ",values = cols(11))+          
  ylim(7.2,8.2)+
  labs(y=paste("Shannon index", sep=""))+
  
  geom_jitter(position=position_jitter(0.1),shape=21,alpha=0.8,size=3)+
  
  theme_bw() +
  theme(plot.title = element_text(size =15, colour = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = 15,color="black"),
        axis.title = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(),
        legend.text =  element_blank(),legend.position = "none",
        legend.title = element_text( size = 15,color="black"))+
  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))
p3

###########################richness

p4<-ggplot(dom,mapping = aes(x =sample, y =Richness_dom,fill=sample,color=sample))+
  
  geom_boxplot(linetype="solid",size=1, width=0.5,alpha=0.6,outlier.shape = NA)+
  scale_fill_manual(name=" ",values = cols(11))+ 
  scale_color_manual(name=" ",values = cols(11))+         
  ylim(5000,7400)+
  labs(y=paste("Observed richness", sep=""))+
  
  geom_jitter(position=position_jitter(0.1),shape=21,alpha=0.8,size=3)+
  
  theme_bw() +
  theme(plot.title = element_text(size =15, colour = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = 15,color="black"),
        axis.title = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(),
        legend.text =  element_blank(),legend.position = "none",
        legend.title = element_text( size = 15,color="black"))+
  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))
p4

############################# Random Forest-SHAP

library(tidyverse)

DOM<- read.csv("DOM2.csv", header = T)
head(DOM)
names(DOM)

DOM_data  <- DOM %>%dplyr::select(MAT,pH,TC,TN,MAP,
                                  Bacterial.Shannon.index,
                                  Fungal.Shannon.index,
                                  Bacterial.PCoA1,
                                  Fungal.PCoA1)



#### Bacterial diversity ####
library(randomForest)
library(rfPermute)
library(ggplot2)
library(ggpubr)
names(DOM_data)

DOM_Shannon_B  <- DOM %>%dplyr::select(MAT,pH,TC,TN,MAP)


set.seed(123)
Rf_bacteria <- randomForest(Bacterial.Shannon.index~MAT+pH+TC+TN+MAP
                            
                            ,
                            
                            data = DOM_data ,
                            ntree=500,
                            importance=TRUE,
                            proximity=TRUE )
summary(Rf_bacteria)


importance_otu.scale <- data.frame(importance(Rf_bacteria, scale = TRUE), check.names = FALSE)
importance_otu.scale


######  SHAP value  ###########
library(fastshap)
library(ggplot2)

shap_bac <- explain(Rf_bacteria, X = subset(DOM_Shannon_B ), nsim = 10, 
                    pred_wrapper = predict)
shap_bac
summary(abs(shap_bac))


p1<-autoplot(shap_bac,fill="#104E8B", width = 0.8)+labs(title = NULL, x = 'Bacierial diversity ', y = 'Mean absolute SHAP value', fill = NULL)+
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.12),breaks = c(0, 0.04,0.08,0.12))+
  scale_x_continuous(expand = c(0, 0.12), limit = c(0, 0.02,0.04,0.06,0.08,0.10,0.12))+
  annotate('text', label = sprintf('italic(R^2) == %.2f', 72.5), 
  x = 1.5, y = 0.07, size = 5, family = "serif",parse = TRUE)+
  annotate('text', label = sprintf('italic(p) < %.2f', 0.05), 
  x = 1.1, y =0.07, size = 5, family = "serif",parse = TRUE)+
  
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  theme(plot.margin = unit(c(0.6,0.5,0.5,0.5),"cm"))+
  theme(legend.position="none")
p1

library(ggpmisc)

formula <- y ~ poly(x, 2, raw = TRUE)

p2<-autoplot(shap_bac,type = "dependence",feature = "Tannin.diversity", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
 
  stat_fit_glance(method = 'lm',
                  method.args = list(formula =y ~ poly(x, 2, raw = TRUE)),
                  
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=4)+
  scale_y_continuous(limits=c(-0.4,0.4))+
  
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
 
  theme(legend.position="none")
p2

formula <- y ~ x
formula <- y ~ poly(x, 2, raw = TRUE)
p4<-autoplot(shap_bac,type = "dependence",feature = "Carbohydrate.diversity", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
 
  
  stat_fit_glance(method = 'lm',
                  
                
                  method.args = list(formula =y ~ poly(x, 2, raw = TRUE)),
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=5)+
  
  scale_y_continuous(limits=c(-0.06,0.06))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  
  theme(legend.position="none")
p4

p3<-autoplot(shap_bac,type = "dependence",feature = "Lipid.diversity", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
  
  stat_fit_glance(method = 'lm',
                  
                  method.args = list(formula =y ~ x),
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=4)+
  scale_y_continuous(limits=c(-0.2,0.2))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  
  theme(legend.position="none")


p3

p5<-autoplot(shap_bac,type = "dependence",feature = "pH", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
 
  stat_fit_glance(method = 'lm',
                  
                  method.args = list(formula =y ~ x),
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=4)+
  scale_y_continuous(limits=c(-0.4,0.4))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
 
  theme(legend.position="none")


p5

####  fungal diversity ####
names(DOM_data)

set.seed(345)
Rf_fungi <- randomForest(Fungal.diversity~
                           MAT+pH+
                           Aromatic.diversity+
                           
                           Carbohydrate.diversity+
                           
                           Lignin.diversity+
                           
                           Lipid.diversity+
                           
                           Protein.diversity+
                           
                           Tannin.diversity,
                         data = DOM_data,
                         ntree=500,
                         importance=TRUE,
                         proximity=TRUE)
summary(Rf_fungi)


importance_otu.scale <- data.frame(importance(Rf_fungi, scale = TRUE), check.names = FALSE)
importance_otu.scale



######  SHAP value  ###########
library(fastshap)
library(ggplot2)

shap_fungi <- explain(Rf_fungi, X = subset(DOM_Shannon_B), nsim = 10, 
                      pred_wrapper = predict)
shap_fungi
summary(abs(shap_fungi))


p6<-autoplot(shap_fungi,fill="#104E8B", width = 0.8)+labs(title = NULL, x = 'Fungal diversity ', y = 'Mean absolute SHAP value', fill = NULL)+
  annotate('text', label = 'Shannon_Bacieria', x = 2, y = 0.095, size = 6,family = "serif") +
  annotate('text', label = sprintf('italic(R^2) == %.2f', 28.5), 
  x = 1.6, y = 0.055, size = 6, family = "serif",parse = TRUE)+
  annotate('text', label = sprintf('italic(p) < %.2f', 0.05), 
  x = 1.2, y = 0.055, size = 6, family = "serif",parse = TRUE)+
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.06),breaks = c(0,0.02,0.04,0.06))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  theme(plot.margin = unit(c(0.6,0.5,0.5,0.5),"cm"))+
  theme(legend.position="none")
p6

library(ggpmisc)

formula <- y ~ poly(x, 2, raw = TRUE)

p7<-autoplot(shap_fungi,type = "dependence",feature = "Lignin.diversity", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed") + 
  stat_fit_glance(method = 'lm',
                  
                  method.args = list(formula =y ~ x),
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=4)+
  scale_y_continuous(limits=c(-0.2,0.2))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))+
  theme(legend.position="none")
p7

formula <- y ~ x
p8<-autoplot(shap_fungi,type = "dependence",feature = "Carbohydrate.diversity", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed") +
  
  stat_fit_glance(method = 'lm',
                  
                  method.args =  list(formula =y ~ x),
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=4)+
  
  scale_y_continuous(limits=c(-0.2,0.2))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))+
  theme(legend.position="none")
p8

p10<-autoplot(shap_fungi,type = "dependence",feature = "Protein.diversity", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed") + 
  stat_fit_glance(method = 'lm',
                  
                  method.args = list(formula =y ~ poly(x, 2, raw = TRUE)),
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=4)+
  scale_y_continuous(limits=c(-0.1,0.3))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))+
  theme(legend.position="none")


p10
formula <- y ~ x


formula <- y ~ poly(x, 2, raw = TRUE)
p9<-autoplot(shap_fungi,type = "dependence",feature = "MAT", X = DOM,shape=21)+
  geom_point(size=1.5,color="#104E8B")+
  stat_smooth(lty=1,se=TRUE, method="lm",formula=formula,color="black",show.legend = F)+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed") + 
  stat_fit_glance(method = 'lm',
                  
                  method.args = list(formula =y ~ poly(x, 2, raw = TRUE)),
                  
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  
                  parse = TRUE,label.x = "left", label.y = "top",size=4)+
  
  
  scale_y_continuous(limits=c(-0.2,0.3))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12))+
  theme(plot.margin = unit(c(0.9,0.8,0.8,0.8),"cm"))+
  theme(legend.position="none")


p9





