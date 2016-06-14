# pFDA for DIET TYPE

# phylogenetic FDA 
# modified from the tutorial by Lars Schmitz, 2014
################################################################################################

# Loading libraries (install first if needed)
library(dplyr)
library(geiger)
library(mda)
library(nnet)
library(ggplot2)
library(ggrepel)
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

# prepare data
ddDiet <- ddA 

# define vector of taxon names
taxaA <- ddA$sp
ddDiet2 <- ddDiet

# Defining groups
gADiet <- ddA  %>% select(DietCond) %>% 
  collect %>% .[["DietCond"]] %>% factor()  

# prepare matrix
ddDiet <- ddDiet %>% select(-DietCond,-sp,-Locomotion)  
ddDiet <- as.matrix(ddDiet,rownames.force = TRUE)
ddDiet <- signif(ddDiet,4)

#which species are uknown
testtaxan <- which(ddA$DietCond=="U")


# creating a dataframe that only contains taxa with known group affiliation
X <- ddDiet
# Performing the discriminant analysis

# read lambda values
lambdaOpts <- read.csv("optimalLambdaVals.csv",stringsAsFactors = F)

# to loop through tree block
pFDAdietOut <- list()
pb   <- txtProgressBar(1, 100, style=3)

for (i in 1:100) {
  pFDAdietOut[[i]] <- phylo.fda.pred(X,gADiet,taxaA,rodTrees[[i]],testtaxan,val=lambdaOpts$DietLambdaVals[i])  # Warning message re: priors can be ignored.
  
  setTxtProgressBar(pb, i)
}


#inputs here are:

# XA all measurements including unknown
# gA groups including unknown
# taxaA vector of species names incl. unknown
# tree (trimmed?)
# testtaxan vector with the index of the uknown species
# val is lambdas

# get the % explained by each dimension *see the range of values min-max
pFDAdietOut %>% lapply("[[",1) %>% simplify2array() %>% apply(1,summary)

# confusion rates for all trees
round(simplify2array(lapply(lapply(pFDAdietOut, "[[", 10),attr,which="error")),4)


# get the discriminant scores
dietVariates <- list()
for (i in 1:100) {
  dietVariates[[i]] <- predict(pFDAdietOut[[i]],type="variates")
}


# get the predicted classes

dietFDAclassesALL <- list()
for (i in 1:100) {
  dietFDAclassesALL[[i]] <- predict(pFDAdietOut[[i]],type="class")
}

# how many different sets of predicted classes
length(unique(dietFDAclassesALL))

# get the counts for different sets of predictions
setsuniquedietAll <- c()
for (i in 1:length(unique(dietFDAclassesALL))){
  setsuniquedietAll[i] <- which(sapply(dietFDAclassesALL, function(z) all(z == unique(dietFDAclassesALL)[[i]]))) %>% length()
}

# print
setsuniquedietAll

# bind the two unique class predictions together
PredDietclasses <- as.data.frame(simplify2array(unique(dietFDAclassesALL)))
# check prediction differences
PredDietclasses$dif <- PredDietclasses$V1==PredDietclasses$V2 #essentially identical

# predictions discrimination
# using the set of predictions that was most frequent

# new copy of original data
ddDiet3 <- ddDiet2 %>% filter(DietCond!="U")

# bind true and predicted values into DF
predictionsDiet <-  cbind.data.frame(true=factor(ddDiet3$DietCond),predicted=dietFDAclassesALL[[which.is.max(setsuniquedietAll)]])

# compare true and predicted
predictionsDiet$same <- predictionsDiet$true == predictionsDiet$predicted
predictionsDiet %>% group_by(true) %>% summary()

predTFDiet <- predictionsDiet %>% group_by(true,same) %>% 
  summarize(nmany=n()) %>%  
  mutate(freq = nmany / sum(nmany))

# percent correctly classified by Diet Category
predTFDiet %>% filter(same==T) %>% mutate(frounded=round(freq,2))


