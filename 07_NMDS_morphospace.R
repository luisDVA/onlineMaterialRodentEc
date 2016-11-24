# nMDS
library(dplyr)

# Reading in the morphology data and assigning rownames
morphCorrected <- read.csv("sizeCorrected.csv",stringsAsFactors = FALSE) %>% dplyr::select(-X)
rownames(morphCorrected) <- morphCorrected$sp
morphCorrected$sp <- NULL

library(cluster)
# dissimilarity matrix
distGmorph <- daisy(morphCorrected,metric="gower",stand=F)

library(MASS)
# run nMDS
# the magic value is one of many that minimizes stress (derived by testing from 0 to 0.5 in 
### increments of 0.02)
nmdsMorph <- sammon(distGmorph,magic=0.44,niter = 200)

# put the two dimensions into a DF
NDMSaxesRod <- as.data.frame(nmdsMorph$points)
NDMSaxesRod$sp <- rownames(NDMSaxesRod)

# Read FDA predictions
FDAdietPred <- read.csv("dietFDAvariatesPreds.csv", stringsAsFactors = FALSE) %>% 
  rename(sp=Row.names)
FDAlocPred <- read.csv("locFDAvariatesPreds.csv", stringsAsFactors = FALSE) %>% 
  rename(sp=Row.names)

# merge
morphospaceRoddiet <- left_join(NDMSaxesRod,FDAdietPred) %>% 
  rename(nMDS1=V1,nMDS2=V2)
morphospaceRodloc <- left_join(NDMSaxesRod,FDAlocPred) %>% 
  rename(nMDS1=V1,nMDS2=V2)


library(ggplot2)
library(ggrepel)

dietplot <- 
  ggplot(morphospaceRoddiet,aes(x=nMDS1,y=nMDS2,fill=dietPred,shape=dietPred,label=miscl))+
  geom_vline(xintercept = 0,color="light grey")+
  geom_hline(yintercept = 0,color="light grey")+
  geom_point(size=5,aes(color=dietPred))+
  geom_text_repel(fontface = 'bold',
                  point.padding = unit(0.8, 'lines'),
                  segment.color = '#555555',
                  segment.size = 0.5,
                  force = 0.5)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),legend.position="none")+
  scale_shape_manual(values = c(25,21,21,23,22))+
  scale_color_manual(values=c("white","white","black","white","white"))+
  scale_fill_manual(values=c("#97ce9c","#cff09e","#cff09e","#5a948b","#0b486b"))+
  ylab("nMDS axis 2")+xlab("nMDS axis 1")

locplot <- 
  ggplot(morphospaceRodloc,aes(x=nMDS1,y=nMDS2,shape=locPred,fill=locPred,label=miscl))+
  geom_point(size=5,aes(color=locPred))+
  geom_vline(xintercept = 0,color="light grey")+
  geom_hline(yintercept = 0,color="light grey")+
  geom_text_repel(fontface = 'bold',
                  point.padding = unit(0.8, 'lines'),
                  segment.color = '#555555',
                  segment.size = 0.5,
                  force = 0.4)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),legend.position="none")+
  scale_color_manual(values=c("white","black","white","white","white","white","white","black","white","black"))+
  scale_shape_manual(values = c(21,21,25,22,22,21,23,23,24,24))+
  scale_fill_manual(values=c("#CFF09E","#CFF09E","blue","#FA7921", "#79BD9A", 
                             "#0E590A", "#FDE74C","#FDE74C","#5BC0EB","#5BC0EB" ))

library(gridExtra)
grid.arrange(dietplot,locplot,ncol=2)
