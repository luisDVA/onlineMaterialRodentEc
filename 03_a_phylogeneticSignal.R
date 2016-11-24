# Phylogenetic signal tests

# load libraries
library(geomorph)
library(dplyr)
library(phytools)

# load measurements
morphData <- read.csv("meansMasses.csv",stringsAsFactors = FALSE) %>%  
                    dplyr::select(-X,-dataProv)

# Reading in the trees
rodTrees <- read.nexus("FStreeBlock.nex")

# reorder the DF to match the tip labels
morphData <- morphData[match(rodTrees[[1]]$tip.label, gsub(" ","_",morphData$sp)),]
# set row names
row.names(morphData) <- gsub(" ","_",morphData$sp)

# split into craniodental/mandibular and external
externalMeas <- morphData %>% dplyr::select(sp,HB,T:UM)
crMeas <- morphData %>% dplyr::select(sp,CBL:ACP)

# assign row names 
row.names(externalMeas) <- gsub(" ","_",externalMeas$sp)
row.names(crMeas) <- gsub(" ","_",crMeas$sp)

# Kmult across tree block

# external measurements #################
KmultExternal <- list()
pb   <- txtProgressBar(1, 100, style=3)

for (i in 1:100) {
  KmultExternal[[i]] <- physignal(log(externalMeas[,2:8]),rodTrees[[i]])
  
  setTxtProgressBar(pb, i)
}
# value of K
extKvals <- KmultExternal %>% lapply("[[",1) %>% simplify2array() 
round(summary(extKvals),2)

# significance
extKpvals <- KmultExternal %>% lapply("[[",2) %>% simplify2array() 
length(which(extKpvals<0.05))

#### craniodental measurements ##############
Kmultcr <- list()
pb   <- txtProgressBar(1, 100, style=3)

for (i in 1:100) {
  Kmultcr[[i]] <- physignal(log(crMeas[,2:8]),rodTrees[[i]])
  
  setTxtProgressBar(pb, i)
}
# value of K
crKvals <- Kmultcr %>% lapply("[[",1) %>% simplify2array()
round(summary(crKvals),2)
# significance
crKpvals <- Kmultcr %>% lapply("[[",2) %>% simplify2array() 
length(which(crKpvals<0.05))
#######################

# optional plotting
KmultVals <- data.frame(crKvals,extKvals)
ggplot(KmultVals,aes(x=crKvals))+geom_histogram()+geom_density()

####################################
# Phylogenetic signal in body size

KbodyMass <- list()

for (i in 1:100) {
  KbodyMass[[i]] <- phylosig(rodTrees[[i]],log(morphData$massGrams),method="K",test=TRUE)
}
# K
bodySizeKvals <- KbodyMass %>% lapply("[[",1) %>% simplify2array() 
round(summary(bodySizeKvals),2)
# significance
bodySizeKpvals <- KbodyMass %>% lapply("[[",2) %>% simplify2array() 
length(which(bodySizeKpvals<0.05))

####################
# Phyl. signal for categorical traits
# CODE BELOW IS FROM Debastiani & Duarte (2016)
## Required packages

require(vegan)
require(FD)

## Description
# Calculate the phylogenetic signal based on a Mantel test, incorporating the Brownian motion evolutionary model.

## Arguments

# tree = Phylogenetic tree, as phylo object.
# traits = Matrix or data frame containing the traits data, with traits as columns and species as rows. Traits can be numeric, ordered, or factor, as the required by the gowdis function. (Symmetric or asymmetric binary variables should be numeric and only contain 0 and 1. character variables will be converted to factor).
# runs = Number of permutations in assessing significance.
# euclidean = Logical argument (TRUE or FALSE) to specify if use transformation of  euclidean properties in pairwise trait dissimilarities (Default euclidean = TRUE).
# sqrtPhylo = Logical argument (TRUE or FALSE) to specify if use square root transformation of phylogenetic distance (Default sqrtPhylo = FALSE).
# checkdata = Logical argument (TRUE or FALSE) to check if species sequence in the trait data follows the same order as in phylogenetic tree (Default checkdata = TRUE).


## Value

# perm.NULL = A vector of permuted Mantel statistic under the null hypothesis that the distances between both matrices are not related.
# perm.BM = A vector of permuted Mantel statistic under the null hypothesis that the traits evolve under a Brownian motion evolutionary model.
# r.Mantel = The Mantel observed statistic.
# p.NULL = The p value from no phylogenetic structure (standard p value in Mantel test).
# p.BM = The p value under simulation of traits from Brownian phylogenetic structure.

# source the evolutionary-model Mantel function
EM.mantel<-function(tree, traits, runs = 999, euclidean= TRUE, sqrtPhylo=FALSE, checkdata = TRUE, ...){
  phylo.dist<-cophenetic(tree)
  if(sqrtPhylo){
    phylo.dist<-sqrt(phylo.dist)
  }
  if(checkdata){
    if(is.null(tree$tip.label)){
      stop("\n Error in tip labels of tree\n")
    }
    if(is.null(rownames(traits))){
      stop("\n Error in row names of traits\n")
    }
    match.names <- match(rownames(traits),rownames(phylo.dist))
    if(sum(is.na(match.names)) > 0){
      stop("\n There are species from traits data that are not on phylogenetic tree\n")
    }
    phylo.dist <- phylo.dist[match.names, match.names]
  }
  if(length(tree$tip.label) > dim(traits)[1]){
    warning("Tree have more species that species in traits data")
  }
  if(dim(phylo.dist)[1] != dim(traits)[1] & checkdata == FALSE){
    stop("\n Different number of species in tree and in traits data, use checkdata = TRUE\n")
  }
  gow.dist<-gowdis(traits, ...)
  if(euclidean){
    gow.sim<-1-gow.dist
    gow.dist<-sqrt(1-gow.sim)
  }
  traits.attr<-attr(gow.dist, "Types", exact = TRUE)
  res.mantel<-mantel(phylo.dist,gow.dist,permutations=runs)
  res.BM<-matrix(NA,runs,1)
  for(k in 1:runs){
    traits_sim<-matrix(NA,length(tree$tip.label),dim(traits)[2])
    rownames(traits_sim)<-tree$tip.label
    for(i in 1:dim(traits)[2]){
      traits_sim[,i]<-rTraitCont(tree,model="BM")
    }
    traits_sim<-decostand(traits_sim,method="standardize",MARGIN=2)
    traits_sim<-as.data.frame(traits_sim)	
    for(i in 1:dim(traits)[2]){
      if(traits.attr[i] == "B" | traits.attr[i] == "A"){
        probs<-sum(traits[,i])/dim(traits)[1]
        threshold<-quantile(traits_sim[,i],probs=1-probs)
        traits_sim[,i]<-ifelse(traits_sim[,i]>=threshold,1,0)
      }
      if(traits.attr[i] == "N" | traits.attr[i] == "O"){
        n.levels<-length(levels(traits[,i]))
        traits.levels<-levels(traits[,i])
        probs<-cumsum(table(traits[,i]))/sum(table(traits[,i]))
        probs<-probs[1:(n.levels-1)]
        threshold<-quantile(traits_sim[,i],probs=probs)
        threshold<-c(min(traits_sim[,i]),threshold,max(traits_sim[,i]))
        temp<-matrix(NA,length(traits_sim[,i]),1)
        for(j in 1:n.levels){
          if(j < n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<threshold[j+1], traits.levels[j],temp)
          }
          if(j == n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<=threshold[j+1], traits.levels[j],temp)
          }
        }
        traits_sim[,i]<-as.factor(temp)
        if(traits.attr[i] == "O"){
          traits_sim[,i]<-ordered(temp,levels=levels(traits[,i]))
        }
      }
    }
    if(checkdata == TRUE){
      match.names <- match(rownames(traits),rownames(traits_sim))
      traits_sim<-traits_sim[match.names,,drop=FALSE]
    }
    gow.dist.BM<-gowdis(traits_sim, ...)
    if(euclidean){
      gow.sim.BM<-1-gow.dist.BM
      gow.dist.BM<-sqrt(1-gow.sim.BM)
    }
    res.mantel.BM<-mantel(phylo.dist,gow.dist.BM,permutations=0)
    res.BM[k,1]<-res.mantel.BM$statistic
  }
  p.BM<-(sum(ifelse(res.BM[,1]>=res.mantel$statistic,1,0))+1)/(runs+1)
  p.NULL<-res.mantel$signif
  r.Mantel<-res.mantel$statistic
  RES<-list(perm.NULL=res.mantel$perm,perm.BM=res.BM[,1],r.Mantel=r.Mantel,p.NULL=p.NULL,p.BM=p.BM)
  return(RES)
}

## load ecological traits
ecolData <- read.csv("DietLocomotion.csv",stringsAsFactors = FALSE) %>% dplyr::select(sp,DietCond,Locomotion)
ecolData$sp <- gsub(" ","_",ecolData$sp)

# subset for diet type
diet <- ecolData %>% dplyr::select(sp,DietCond)
# remove species with Unknown diet types
diet <- diet %>% dplyr::filter(DietCond!="U")

######## trim trees
# define keep.tip function by L. Revell
keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

# drop extra tips from trees
# tips to keep
tipsWant <-diet$sp
rodTreesDiet <-lapply(rodTrees,keep.tip,tip=tipsWant)
class(rodTreesDiet) <- "multiPhylo"

# reorder Diet DF to match tip labels
diet <- diet[match(rodTreesDiet[[1]]$tip.label, diet$sp),]
rownames(diet) <- diet$sp

# keep only the categories
diet <- diet[,2,drop=F]
diet$DietCond <- factor(diet$DietCond)

# run EMmantel
emMantelDiet <- list()
pb   <- txtProgressBar(1, 100, style=3)

for (i in 1:100) {
  emMantelDiet[[i]] <- EM.mantel(rodTreesDiet[[i]],diet)
  
  setTxtProgressBar(pb, i)
}

# emMantel correlations, diet
emMantelrsDiet <- emMantelDiet %>% lapply("[[",3) %>% simplify2array() 
round(summary(emMantelrsDiet),2)
# significance
emMantelpvalsBMdiet <- emMantelDiet %>% lapply("[[",5) %>% simplify2array() 
length(which(emMantelpvalsBMdiet<0.05))


########## EM Mantel for locomotion mode #######

# subset for locomotion
locomotion <- ecolData %>% dplyr::select(sp,Locomotion)
# remove species with Unknown diet types
locomotion <- locomotion %>% dplyr::filter(Locomotion!="U")

# trim trees
# drop extra tips from trees
# tips to keep
tipsWantLoc <-locomotion$sp
rodTreesLoc <-lapply(rodTrees,keep.tip,tip=tipsWantLoc)
class(rodTreesLoc) <- "multiPhylo"

# reorder Diet DF to match tip labels
locomotion <- locomotion[match(rodTreesLoc[[1]]$tip.label, locomotion$sp),]
rownames(locomotion) <- locomotion$sp

# keep only the categories
locomotion <- locomotion[,2,drop=F]
locomotion$Locomotion <- factor(locomotion$Locomotion)

# run EMmantel

# run EMmantel
emMantelLoc <- list()
pb   <- txtProgressBar(1, 100, style=3)

for (i in 1:100) {
  emMantelLoc[[i]] <- EM.mantel(rodTreesLoc[[i]],locomotion)
  
  setTxtProgressBar(pb, i)
}

# emMantel correlations, locomotion
emMantelrsLoc <- emMantelLoc %>% lapply("[[",3) %>% simplify2array() 
round(summary(emMantelrsLoc),2)
# significance
emMantelpvalsBMloc <- emMantelLoc %>% lapply("[[",5) %>% simplify2array() 
length(which(emMantelpvalsBMloc<0.05))


# plot diet and locomotion next to the trees for figure S1
library(ggtree)

# locomotor modes
locomotion$taxa <- row.names(locomotion)
locannotate <- locomotion %>% dplyr::select(taxa,Locomotion)
locplot <- ggtree(rodTreesLoc[[15]],layout = "fan") %<+% locannotate + geom_text(aes(label=Locomotion))
ggsave("loco2.eps",locplot)

#diet types
diet$taxa <- row.names(diet)
dietannotate <- diet %>% dplyr::select(taxa,DietCond)
dietannotate$DietCond <- gsub("H","",dietannotate$DietCond)
dietplot <- ggtree(rodTreesDiet[[15]],layout = "fan") %<+% dietannotate + geom_text(aes(label=DietCond))
ggsave("diet2.eps",dietplot)


