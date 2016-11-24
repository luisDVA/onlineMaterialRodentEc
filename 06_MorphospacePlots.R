######### FDA and canonical plots

# required packages
library(dplyr)
library(mda)
library(ggplot2)
library(ggrepel)


# Read ecology data
ecolData <- read.csv("DietLocomotion.csv",stringsAsFactors = FALSE)%>% 
  select(sp,DietCond,Locomotion)
ecolData$sp <- gsub(" ","_",ecolData$sp)

# Reading in the morphology data and assigning rownames
morphCorrected <- read.csv("sizeCorrected.csv",stringsAsFactors = FALSE) %>% select(-X)
rownames(morphCorrected) <- morphCorrected$sp
morphCorrected[,2:13] <- signif(morphCorrected[,2:13],3)

#merge all (sort=F keeps the table in the right order)
ddA <- merge(morphCorrected,ecolData,sort=FALSE)

####### training Data
# DIET
ddDiet <- ddA %>% filter(DietCond!="U") %>% select(-Locomotion)
row.names(ddDiet) <- ddDiet$sp
ddDiet$sp <- NULL
# Locomotion
ddLoc <- ddA %>% filter(Locomotion!="U") %>% select(-DietCond)
row.names(ddLoc) <- ddLoc$sp
ddLoc$sp <- NULL

######### unknown ecology taxa
#diet
dduDiet <- ddA %>% filter(DietCond=="U") %>% select(-Locomotion)
row.names(dduDiet) <- dduDiet$sp
dduDiet$sp <- NULL
#locomotion
dduLoc <- ddA %>% filter(Locomotion=="U") %>% select(-DietCond)
row.names(dduLoc) <- dduLoc$sp
ddLoc$sp <- NULL

######### run fda on training data, all Chars, DIET
# fit FDA
dietFDA <- fda(factor(DietCond)~.,data=ddDiet)

# predict variates for training data
dietFDAvariatesTR <- predict(dietFDA,type="variates",newdata = ddDiet)
# predict classes for training data
dietFDAclassesTR <- predict(dietFDA,type="class")
# variates for testing data (unknown diet)
dietFDAvariatesU <- predict(dietFDA,type="variates",newdata = dduDiet)
dietFDAclassesU <- predict(dietFDA,type="class",newdata = dduDiet)

# bind/fortify
dietMorphospace <- data.frame(rbind(dietFDAvariatesTR,dietFDAvariatesU),
                   c(as.character(dietFDAclassesTR),as.character(dietFDAclassesU))) %>% 
          select(DF1=1,DF2=2,DF3=3,dietPred=4)

# to label misclassified sp
actualDiets <- ddA %>% select(sp,DietCond)
#bind true onto df
dietMorphospace <-  merge(dietMorphospace,actualDiets,by.x="row.names",by.y="sp")
dietMorphospace$dietPred <- as.character(dietMorphospace$dietPred)
# identify misclassified sp
dietMorphospace$miscl <- dietMorphospace$dietPred==dietMorphospace$DietCond

# to label sp with unknown ecology 
dietMorphospace$dietPred[which(dietMorphospace$DietCond=="U")] <- paste(as.character(dietMorphospace$dietPred[which(dietMorphospace$DietCond=="U")]),"p",sep = "")
dietMorphospace$dietPred <- factor(dietMorphospace$dietPred)

# tidy up label column
dietMorphospace$miscl[which(dietMorphospace$miscl==T)] <- NA
dietMorphospace$miscl[which(dietMorphospace$DietCond=="U")] <- NA
dietMorphospace$miscl[which(dietMorphospace$miscl==F)] <- dietMorphospace$DietCond[which(dietMorphospace$miscl==F)]

# write to disk
write.csv(dietMorphospace,file="dietFDAvariatesPreds.csv")

# plot figure
 
ggplot(dietMorphospace,aes(x=DF1,y=DF2,fill=dietPred,shape=dietPred,label=miscl))+
  geom_vline(xintercept = 0,color="light grey")+
  geom_hline(yintercept = 0,color="light grey")+
  geom_point(size=5,aes(color=dietPred))+
  geom_text_repel(fontface = 'bold',
                  point.padding = unit(0.8, 'lines'),
                  segment.color = '#555555',
                  segment.size = 0.5,
                  force = 0.5)+
  theme_bw()+
  theme(panel.grid.major=element_blank())+
  scale_shape_manual(values = c(25,21,21,23,22))+
  scale_color_manual(values=c("white","white","black","white","white"))+
  scale_fill_manual(values=c("#97ce9c","#cff09e","#cff09e","#5a948b","#0b486b"))

# to export as vector graphics
postscript("fig2.eps", height = 5.5, width = 7.5,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)


dev.off()


###########################################################

#### fda on training data, all characters, Locomotion
# fit the FDA
locFDA <- fda(factor(Locomotion)~.,data=ddLoc)

# variates for training data
locFDAvariates <- predict(locFDA,type="variates",newdata = ddLoc)
# classes for training data
locFDAclasses <- predict(locFDA,type="class")
# variates for testing data (unknown loc)
locFDAvariatesU <- predict(locFDA,type="variates",newdata = dduLoc)
#classes for testing data (unknown loc)
locFDAclassesU <- predict(locFDA,type="class",newdata = dduLoc)

# bind into DF
locMorphospace <- data.frame(rbind(locFDAvariates,locFDAvariatesU),
                              c(as.character(locFDAclasses),as.character(locFDAclassesU))) %>% 
                  select(DF1=1,DF2=2,DF3=3,DF4=4,DF5=5,DF6=6,locPred=7)

# to label misclassified sp
actualLocs <- ddA %>% select(sp,Locomotion)
locMorphospace <-  merge(locMorphospace,actualLocs,by.x="row.names",by.y="sp")
locMorphospace$locPred <- as.character(locMorphospace$locPred)
# to label sp with unknown ecology
locMorphospace$miscl <- locMorphospace$locPred==locMorphospace$Locomotion
locMorphospace$locPred[which(locMorphospace$Locomotion=="U")] <- paste(as.character(locMorphospace$locPred[which(locMorphospace$Locomotion=="U")]),"p",sep = "")
locMorphospace$locPred <- factor(locMorphospace$locPred)
# tidy up labels column
locMorphospace$miscl[which(locMorphospace$miscl==T)] <- NA
locMorphospace$miscl[which(locMorphospace$Locomotion=="U")] <- NA
locMorphospace$miscl[which(locMorphospace$miscl==F)] <- locMorphospace$Locomotion[which(locMorphospace$miscl==F)]

# write to disk
write.csv(locMorphospace,file="locFDAvariatesPreds.csv")

# plotting

ggplot(locMorphospace,aes(x=DF1,y=DF2,shape=locPred,fill=locPred,label=miscl))+
  geom_point(size=5,aes(color=locPred))+
  geom_vline(xintercept = 0,color="light grey")+
  geom_hline(yintercept = 0,color="light grey")+
  geom_text_repel()+
  theme_bw()+
  theme(panel.grid.major=element_blank())+
  scale_color_manual(values=c("white","black","white","white","white","white","white","black","white","black"))+
  scale_shape_manual(values = c(21,21,25,22,22,21,23,23,24,24))+
  scale_fill_manual(values=c("#CFF09E","#CFF09E","blue","#FA7921", "#79BD9A", 
                              "#0E590A", "#FDE74C","#FDE74C","#5BC0EB","#5BC0EB" ))
