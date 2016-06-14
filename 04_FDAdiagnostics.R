######### FDA diagnostics
# flexible discriminant analysis

# packages (install if necessary)
library(dplyr)
library(mda)
library(ggplot2)
library(DiscriMiner)

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
# variance explained
round(dietFDA$percent.explained,2)
# training misclassification error
round(1-attributes(dietFDA$confusion)$error,2)

#use linDA to get crossvalidated error

linDA(ddDiet[,1:12],ddDiet$DietCond,validation = "crossval")$error_rate

# predict classes for training data
dietFDAclassesTR <- predict(dietFDA,type="class")

# bind true and predicted values into DF
predictionsDiet <-  cbind.data.frame(true=factor(ddDiet$DietCond),predicted=dietFDAclassesTR)

# compare true and predicted
predictionsDiet$same <- predictionsDiet$true == predictionsDiet$predicted
predictionsDiet %>% group_by(true) %>% summary()

predTFDiet <- predictionsDiet %>% group_by(true,same) %>% 
  summarize(nmany=n()) %>%  
  mutate(freq = nmany / sum(nmany))

# pecent correctly classified by category
predTFDiet %>% filter(same==T) %>% mutate(frounded=round(freq,2))

######### run fda on training data, all Chars, Locomotion
# fit FDA
locFDA <- fda(factor(Locomotion)~.,data=ddLoc)
# variance explained
round(locFDA$percent.explained,2)
# training misclassification error
round(1-attributes(locFDA$confusion)$error,2)

#use linDA to get crossvalidated error

linDA(ddLoc[,1:12],ddLoc$Locomotion,validation = "crossval")$error


# predict classes for training data
locFDAclassesTR <- predict(locFDA,type="class")


# bind true and predicted values into DF
predictionsLoc <-  cbind.data.frame(true=factor(ddLoc$Locomotion),predicted=locFDAclassesTR)
rownames(predictionsLoc) <- rownames(ddLoc)
# compare true and predicted
predictionsLoc$same <- predictionsLoc$true == predictionsLoc$predicted
predictionsLoc %>% group_by(true) %>% summary()

predTFLoc <- predictionsLoc %>% group_by(true,same) %>% 
  summarize(nmany=n()) %>%  
  mutate(freq = nmany / sum(nmany))

# pecent correctly classified by category
predTFLoc %>% filter(same==T) %>% mutate(frounded=round(freq,2))


