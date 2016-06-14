# To read the mammal phylogenies from the Faurby and Svenning MPE publication
# http://bios.au.dk/en/about-bioscience/organisation/ecoinformatics-and-biodiversity/data/
library(ape)

# read tree blocks FS2015 (download from website and put in working directory)
resolvedPh1 <- read.nexus("Fully_resolved_phylogeny_1.nex")
resolvedPh2 <- read.nexus("Fully_resolved_phylogeny_2.nex")
resolvedPh3 <- read.nexus("Fully_resolved_phylogeny_3.nex")

# pick at random 
randomtrees <- sort(sample(1:1000,100,replace = FALSE))
from1stBlock <- randomtrees[which(randomtrees<=333)]
from2ndBlock <- randomtrees[which(randomtrees>333&randomtrees<666)] -333 
from3rdBlock <- randomtrees[which(randomtrees>666)]-666 

# extract from blocks
block1 <- resolvedPh1[from1stBlock]
block2 <- resolvedPh2[from2ndBlock]
block3 <- resolvedPh3[from3rdBlock]

# bind into tree block
rodTrees <- c(block1,block2,block3)

# read spmeans data
morphData <- read.csv("meansMasses.csv",stringsAsFactors = FALSE)


# define keep.tip function by L. Revell
keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

# drop extra tips from trees

# tips to keep
# check if all are present
tipsWant <-gsub(" ","_",morphData$sp)

RodTreesF <-lapply(rodTrees,keep.tip,tip=tipsWant)
class(RodTreesF) <- "multiPhylo"
# write as Nexus file
write.nexus(RodTreesF,file="FStreeBlock.nex")



