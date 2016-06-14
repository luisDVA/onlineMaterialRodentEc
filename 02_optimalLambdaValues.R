# IDENTIFYING OPTIMAL LAMBDA VALUES
# packages and functions (install if needed)
library(dplyr)
library(geiger)
library(ape)
library(class)
library(geiger)
library(lattice)
library(mda)
library(nnet)
library(ggplot2)
source("phylo.fda.v0.2noPlotting.R") # modified to remove plotting


# Reading in the tree block
rodTrees <- read.nexus("FStreeBlock.nex")

# Read diet and locomotion data
ecolData <- read.csv("DietLocomotion.csv",stringsAsFactors = FALSE) %>% 
  select(sp,DietCond,Locomotion)
ecolData$sp <- gsub(" ","_",ecolData$sp)

# Reading in the morphology data and assigning rownames
morphCorrected <- read.csv("sizeCorrected.csv",stringsAsFactors = FALSE) %>% select(-X)
rownames(morphCorrected) <- morphCorrected$sp

#merge all (sort=F keeps the table in the right order)
ddA <- merge(morphCorrected,ecolData,sort=FALSE)
row.names(ddA) <- ddA$sp
ddA <- ddA %>% tbl_df()

# define keep.tip function by L. Revell
keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

## dropping extra tips from trees

# tips to keep
# for DIET
DietTipsWant <- ddA %>% filter(DietCond!="U") %>% select(sp) %>% collect %>% .[["sp"]]
LocTipsWant <- ddA %>% filter(Locomotion!="U") %>% select(sp) %>% collect %>% .[["sp"]]

# keep only the wanted tips in the tree block
# for diet
RodTreesFDIET <-lapply(rodTrees,keep.tip,tip=DietTipsWant)
class(RodTreesFDIET) <- "multiPhylo"
# for locomotion
RodTreesFLOC <-lapply(rodTrees,keep.tip,tip=LocTipsWant)
class(RodTreesFLOC) <- "multiPhylo"

##### Data frames with only complete ecol. data

# for diet
ddDiet <- ddA %>%  filter(DietCond!="U") %>% 
  select(-Locomotion,-DietCond)  
row.names(ddDiet) <- ddDiet$sp
ddDiet$sp <- NULL
ddDiet <- as.matrix(ddDiet,rownames.force = TRUE)
ddDiet <- signif(ddDiet,4)

# for locomotion
ddLoc <- ddA %>%  filter(Locomotion!="U") %>% 
  select(-Locomotion,-DietCond)  
row.names(ddLoc) <- ddLoc$sp
ddLoc$sp <- NULL
ddLoc <- as.matrix(ddLoc,rownames.force = TRUE)
ddLoc <- signif(ddLoc,4)


# Defining groups
grDiet <- ddA %>%  filter(DietCond!="U") %>% select(DietCond) %>% 
  collect %>% .[["DietCond"]] %>% factor()  

grLoc <- ddA %>%  filter(Locomotion!="U") %>% select(Locomotion) %>% 
  collect %>% .[["Locomotion"]] %>% factor()  

#### Loop through tree blocks

# first for diet

optLamDiet <- list()
pb   <- txtProgressBar(1, 100, style=3)
X <- ddDiet

for (i in 1:100)
{
  
  optLamDiet[[i]] <- optLambda(X,grps = grDiet,RodTreesFDIET[[i]])
  setTxtProgressBar(pb, i)
}

# bind the lambda vs RSS tables
RSStableDiet <- data.frame(apply(simplify2array(lapply(optLamDiet, "[[", 2)),2,rbind))
# add index of trees
RSStableDiet$tree <- rep(1:100,each=101)

#optional, visualize Lamda values
# ggplot(RSStableDiet,aes(x = Lambda,y=RSS,color=tree))+geom_point()

# put all the optimal lambda values in a vector
DietLambdaVals <- simplify2array(lapply(optLamDiet, "[[", 1))



############ REPEAT FOR LOCOMOTION #######
optLamLoc <- list()
pb   <- txtProgressBar(1, 100, style=3)
X <- ddLoc
for (i in 1:100)
{
  
  optLamLoc[[i]] <- optLambda(X,grLoc,RodTreesFLOC[[i]])
  setTxtProgressBar(pb, i)
}

# bind the lambda vs RSS tables
RSStableLoc <- data.frame(apply(simplify2array(lapply(optLamLoc, "[[", 2)),2,rbind))
# add index of trees
RSStableLoc$tree <- factor(rep(1:100,each=101))

# for rounding
roundRSSTabLoc <- RSStableLoc
roundRSSTabLoc$RSS <- round(roundRSSTabLoc$RSS,3)

# visualize
ggplot(RSStableLoc,aes(x = Lambda,y=RSS,color=tree))+geom_point()

# optimal values
RSStableLoc %>% group_by(tree) %>% slice(Lambda)
roundRSSTabLoc %>% group_by(tree) %>% slice(which.min(RSS))


# put all the optimal lambda values in a vector
LocLambdaVals <- simplify2array(lapply(optLamLoc, "[[", 1))

OptimalLambdaVals <- cbind.data.frame(tree=seq(1:100),DietLambdaVals,LocLambdaVals)
summary(OptimalLambdaVals)

# Confidence intervals (for reporting)
qnorm(0.975)*sd(OptimalLambdaVals$DietLambdaVals)/sqrt(100)  
qnorm(0.975)*sd(OptimalLambdaVals$LocLambdaVals)/sqrt(100)  

# Write table of optimal Values
write.csv(OptimalLambdaVals,"optimalLambdaVals.csv")
