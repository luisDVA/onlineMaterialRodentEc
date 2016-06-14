# phylogenetic FDA for Locomotion Mode
# modified from the tutorial by Lars Schmitz, 2014
################################################################################################

# Loading libraries (install first if needed)
library(dplyr)
library(ape)
library(class)
library(geiger)
library(lattice)
library(mda)
library(nnet)
library(ggplot2)
source("phylo.fda.v0.2noPlotting.R")

# Reading in the trees
rodTrees <- read.nexus("FStreeBlock.nex")

# Reading in the data and assigning row names
morphCorrected <- read.csv("sizeCorrected.csv",stringsAsFactors = FALSE) %>% select(-X)
rownames(morphCorrected) <- morphCorrected$sp

# adding ecological data
# read diet data
ecolData <- read.csv("DietLocomotion.csv",stringsAsFactors = FALSE) %>% select(sp,DietCond,Locomotion)
ecolData$sp <- gsub(" ","_",ecolData$sp)

#merge all (sort=F keeps the table in the right order)
ddA <- merge(morphCorrected,ecolData,sort=FALSE)
row.names(ddA) <- ddA$sp


# define keep.tip function by L. Revell
keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

# drop extra tips from trees

# tips to keep
# for LOCOMOTION
LocTipsWant <- ddA %>% filter(Locomotion!="U") %>% select(sp) %>% collect %>% .[["sp"]]

# keep only the wanted tips in the tree block
# for diet
# for locomotion
RodTreesFLOC <-lapply(rodTrees,keep.tip,tip=LocTipsWant)
class(RodTreesFLOC) <- "multiPhylo"


##### Data frames with only complete ecol. data

ddLoc <- ddA %>%  filter(Locomotion!="U") 
ddLoc2 <- ddLoc
# Defining groups
gALoc <- ddA %>%  filter(Locomotion!="U") %>% select(Locomotion) %>% 
  collect %>% .[["Locomotion"]] %>% factor()  

ddLoc <- ddLoc %>% select(-DietCond,-sp,-Locomotion)  
ddLoc <- as.matrix(ddLoc,rownames.force = TRUE)
ddLoc <- signif(ddLoc,4)


# creating a dataframe that only contains taxa with known group affiliation
X <- ddLoc
# Performing the discriminant analysis


# read lambda values
lambdaOpts <- read.csv("optimalLambdaVals.csv",stringsAsFactors = F)

# to loop across the tree block
pFDALocOut <- list()
pb   <- txtProgressBar(1, 100, style=3)

for (i in 1:100) {
  pFDALocOut[[i]] <- phylo.fda(X,gALoc,RodTreesFLOC[[i]],val=lambdaOpts$LocLambdaVals[i],keep.fitted = T)  # Warning message re: priors can be ignored.
  
  setTxtProgressBar(pb, i)
}



# get the % explained by each dimension
pFDALocOut %>% lapply("[[",1) %>% simplify2array() %>% apply(1,summary)

# confusion rates for all trees
round(simplify2array(lapply(lapply(pFDALocOut, "[[", 10),attr,which="error")),4)


# get the discriminant scores
LocVariates <- list()
for (i in 1:100) {
  LocVariates[[i]] <- predict(pFDALocOut[[i]],type="variates")
}

# unknown ecologies
dduLoc <- ddA %>% filter(Locomotion=="U") %>% select(-DietCond)
row.names(dduLoc) <- dduLoc$sp
dduLoc$sp <- NULL


# get the predicted classes

locFDAclassesALL <- list()
for (i in 1:100) {
  locFDAclassesALL[[i]] <- predict(pFDALocOut[[i]],type="class")
}

# how many different sets of predicted classes
length(unique(locFDAclassesALL))
# get the counts for different sets of predictions
setsuniqueLocAll <- c()
for (i in 1:length(unique(locFDAclassesALL))){
  setsuniqueLocAll[i] <- which(sapply(locFDAclassesALL, function(z) all(z == unique(locFDAclassesALL)[[i]]))) %>% length()
}
#print
setsuniqueLocAll

# bind the unique class predictions together
PredLocclasses <- as.data.frame(simplify2array(unique(locFDAclassesALL)))
# check prediction differences
PredLocclasses$dif <- PredLocclasses$V1==PredLocclasses$V2

# predictions discrimination
#### using the set of predictions that was most frequent

# new copy of original data
ddLoc3 <- ddLoc2 %>% filter(Locomotion !="U")

# bind true and predicted values into DF
predictionsLoc <-  cbind.data.frame(true=factor(ddLoc3$Locomotion),predicted=locFDAclassesALL[[which.is.max(setsuniqueLocAll)]])
predictionsCharLoc <-  cbind.data.frame(true=ddLoc3$Locomotion,predicted=as.character(locFDAclassesALL[[which.is.max(setsuniqueLocAll)]]))


# compare true and predicted
predictionsCharLoc$same <- predictionsCharLoc$true == predictionsCharLoc$predicted
predictionsCharLoc %>% group_by(true) %>% summary()

predTFLoc <- predictionsCharLoc %>% group_by(true,same) %>% 
  summarize(nmany=n()) %>%  
  mutate(freq = nmany / sum(nmany))

# percent correctly classified by Diet Category
predTFLoc %>% filter(same==T) %>% mutate(frounded=round(freq,2))

