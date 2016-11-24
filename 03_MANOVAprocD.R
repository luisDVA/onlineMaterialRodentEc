# Procrustes Distance MANOVA (prior to Canonical Variates Analysis)

library(dplyr)
library(geomorph)

# Read diet data
ecolData <- read.csv("DietLocomotion.csv",stringsAsFactors = FALSE)%>% 
  select(sp,DietCond,Locomotion)
ecolData$sp <- gsub(" ","_",ecolData$sp)

# Reading in the morphology data and assigning rownames
morphCorrected <- read.csv("sizeCorrected.csv",stringsAsFactors = FALSE) %>% select(-X)
rownames(morphCorrected) <- morphCorrected$sp

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

# fit the PGLS for diet
dietProcDfit <- list()
summDietProcDfits <- list()

for (i in 1:100) {
dietProcDfit[[i]] <- procD.pgls(ddDiet[,1:12]~factor(ddDiet$DietCond),phy=rodTreesDiet[[i]])
summDietProcDfits[[i]] <- summary(dietProcDfit[[i]])
}

# p value for the diet anovas
dietMANOVApvals <- summDietProcDfits %>% lapply("[[",1) %>% lapply("[[",7) %>% lapply("[[",1) %>% simplify2array()
length(which(dietMANOVApvals<0.05)) ## for reporting

# fit the PGLS for locomotion
locProcDfit <- list()
summLocProcDfits <- list()

for (i in 1:100) {
  locProcDfit[[i]] <- procD.pgls(ddLoc[,1:12]~factor(ddLoc$Locomotion),phy=rodTreesLoc[[i]])
  summLocProcDfits[[i]] <- summary(locProcDfit[[i]])
}

# p value for the diet anovas
locMANOVApvals <- summLocProcDfits %>% lapply("[[",1) %>% lapply("[[",7) %>% lapply("[[",1) %>% simplify2array()
length(which(locMANOVApvals<0.01)) # for reporting

