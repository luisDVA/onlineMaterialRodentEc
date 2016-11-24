# body size ANOVA
# load libraries
library(phytools)
library(dplyr)

# Read ecology data
ecolData <- read.csv("DietLocomotion.csv",stringsAsFactors = FALSE)%>% 
  select(sp,DietCond,Locomotion)
ecolData$sp <- gsub(" ","_",ecolData$sp)

# read morphology data
morphData <- read.csv("meansMasses.csv",stringsAsFactors = FALSE)
morphData$sp <- gsub(" ","_",morphData$sp)

#merge all (sort=F keeps the table in the right order)
ddA <- merge(ecolData,morphData,by="sp",sort=F)

####### training Data
# DIET
ddDiet <- ddA %>% filter(DietCond!="U") %>% filter(DietCond!="I") %>% select(-Locomotion)
row.names(ddDiet) <- ddDiet$sp

# Locomotion
ddLoc <- ddA %>% filter(Locomotion!="U") %>% select(-DietCond)
row.names(ddLoc) <- ddLoc$sp

# Reading in the trees
rodTrees <- read.nexus("FStreeBlock.nex")

######## trim trees
# define keep.tip function by L. Revell
keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

# drop extra tips from trees
# tips to keep
### for diet
tipsWantDiet <- rownames(ddDiet)
rodTreesDiet <-lapply(rodTrees,keep.tip,tip=tipsWantDiet)
class(rodTreesDiet) <- "multiPhylo"
### for locomotion
tipsWantLoc <- rownames(ddLoc)
rodTreesLoc <-lapply(rodTrees,keep.tip,tip=tipsWantLoc)
class(rodTreesLoc) <- "multiPhylo"

# fit the PGLS for diet and body size

# sort the DF with tip labels
ddDiet <- ddDiet[match(rodTreesDiet[[1]]$tip.label,ddDiet$sp),]

# Phyl ANOVA for body mass (diet) looped across the tree set
phylANOVAdietBodyMfit <- list()

for (i in 1:100) {
  phylANOVAdietBodyMfit[[i]]<-phylANOVA(rodTreesDiet[[i]],ddDiet$DietCond,log(ddDiet$massGrams),posthoc = F)
}
# get the significance testing values from the lists and count the outputs
phylANOVAdietBMpvals <- phylANOVAdietBodyMfit %>% lapply("[[",2) %>% simplify2array()
length(which(phylANOVAdietBMpvals>0.05))

# Phyl ANOVA for CBL (diet)
phylANOVAdietCBLfit <- list()

for (i in 1:100) {
  phylANOVAdietCBLfit[[i]]<-phylANOVA(rodTreesDiet[[i]],ddDiet$DietCond,log(ddDiet$CBL),posthoc = F)
}
# get the significance testing values from the lists and count the outputs
phylANOVAdietCBLpvals <- phylANOVAdietCBLfit %>% lapply("[[",2) %>% simplify2array()
length(which(phylANOVAdietCBLpvals>0.05))

# Phyl ANOVA for HB (diet)
phylANOVAdietHBfit <- list()

for (i in 1:100) {
  phylANOVAdietHBfit[[i]]<-phylANOVA(rodTreesDiet[[i]],ddDiet$DietCond,log(ddDiet$HB),posthoc = F)
}
# get the significance testing values from the lists and count the outputs
phylANOVAdietHBpvals <- phylANOVAdietHBfit %>% lapply("[[",2) %>% simplify2array()
length(which(phylANOVAdietHBpvals>0.05))

#### Locomotion
# sort the DF with tip labels
ddLoc <- ddLoc[match(rodTreesLoc[[1]]$tip.label,ddLoc$sp),]

# Phyl ANOVA for body mass (locomotion)
phylANOVAlocBodyMfit <- list()

for (i in 1:100) {
  phylANOVAlocBodyMfit[[i]]<-phylANOVA(rodTreesLoc[[i]],ddLoc$Locomotion,log(ddLoc$massGrams),posthoc = F)
}
# get the significance testing values from the lists and count the outputs
phylANOVAlocBMpvals <- phylANOVAlocBodyMfit %>% lapply("[[",2) %>% simplify2array()
length(which(phylANOVAlocBMpvals>0.05))

# Phyl ANOVA for CBL (locomotion)
phylANOVAlocCBLfit <- list()

for (i in 1:100) {
  phylANOVAlocCBLfit[[i]]<-phylANOVA(rodTreesLoc[[i]],ddLoc$Locomotion,log(ddLoc$CBL),posthoc = F)
}
# get the significance testing values from the lists and count the outputs
phylANOVAlocCBLpvals <- phylANOVAlocCBLfit %>% lapply("[[",2) %>% simplify2array()
length(which(phylANOVAlocCBLpvals>0.05))

# Phyl ANOVA for HB (diet)
phylANOVAdietHBfit <- list()

for (i in 1:100) {
  phylANOVAdietHBfit[[i]]<-phylANOVA(rodTreesDiet[[i]],ddDiet$DietCond,log(ddDiet$HB),posthoc = F)
}
# get the significance testing values from the lists and count the outputs
phylANOVAdietHBpvals <- phylANOVAdietHBfit %>% lapply("[[",2) %>% simplify2array()
length(which(phylANOVAdietHBpvals>0.05))


