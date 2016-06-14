# library preload
library(phytools)
library(dplyr)
options(digits=3)
# read treeblock
treeBlock <- read.nexus("FStreeBlock.nex")
# read morphology data
morphData <- read.csv("meansMasses.csv",stringsAsFactors = FALSE) %>% select(-X,-CBL,-HB,-dataProv)
morphData$sp <- gsub(" ","_",morphData$sp)

# vector for x, the independent size variable
sizeVec <- morphData$massGrams
names(sizeVec) <- morphData$sp

# matrix for Y, all other dependent variables
row.names(morphData) <- morphData$sp
Y <- morphData %>% select(-sp,-ACP,-massGrams) %>% as.matrix()


# loop through tree block
resBlock <- list()
# progress bar
pb   <- txtProgressBar(1, 100, style=3)

for (i in 1:100)
{
  
  resBlock[[i]] = phyl.resid(treeBlock[[i]], log(sizeVec), log(Y), method="lambda")
  setTxtProgressBar(pb, i)
}
# OPTIONAL: save result because for loop takes a long time
#saveRDS(object = resBlock,file="loopResult.rds") 
#resBlock <- readRDS("loopResult.rds") #in case we need to load it later

# model estimates 
betasSC <- apply(simplify2array(lapply(resBlock, "[[", 1)), c(1,2), mean) %>% round(.,3) %>% t() 
# Confidence Intervals
CIsbetas <- apply(simplify2array(lapply(resBlock, "[[", 1)), c(1,2), function(x)qnorm(0.975)*sd(x)/sqrt(100)) %>% round(.,3) %>% t() 

# bind estimates and CIs
betaCIs <- cbind(betasSC,CIsbetas)
colnames(betaCIs) <- c("int","slope","int 95% CI","slope 95% CI")
betaCIs <- as.data.frame(betaCIs) %>% select(1,3,2,4)
write.csv(betaCIs,quote = F)

# summarize residuals matrices across the entire tree block
# fourth element of the list 
residAllTrees <- apply(simplify2array(lapply(resBlock, "[[", 4)), c(1,2), mean)  %>% signif(3) %>% as.data.frame()
residAllTrees$sp <- row.names(residAllTrees)

# add ACP
acpDF <- morphData %>% select(sp,ACP)
morphCorrected <- left_join(residAllTrees,acpDF) %>% select(sp,everything())

#scale for better numerical fitting
morphCorrected[,2:13] <- scale(morphCorrected[,2:13])

# save table
write.csv(morphCorrected,file="sizeCorrected.csv")
