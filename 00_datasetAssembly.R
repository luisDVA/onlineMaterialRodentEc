##### Assembling the measurements dataset

library(dplyr)
library(tidyr)


# final measurements database
rodVals <- read.csv("measurementsFinal.csv",stringsAsFactors = F) %>% 
  select(-ID,-specimenID,-W..g.)

# summarize
rodSummary <- rodVals %>% group_by(sp,type) %>% tally()

# drop species with incomplete data (only crandiodental or external but not both)
# using the previously summarized data
incompleteSpecies <- rodSummary %>% count(sp) %>% filter(n<2)
rodValsC<- anti_join(rodVals,incompleteSpecies)

# split into craniodental and external measurements

RodCR <- rodValsC %>% filter(type=="s") %>% select(sp,CBL:ACP) 
RodEX <- rodValsC %>% filter(type=="e") %>% select(sp,HB:FF)

#calculate HB when it's missing
RodEX$HB <- ifelse (is.na(RodEX$HB),RodEX$TL-RodEX$T, RodEX$HB)
#calculate total length when it's missing
RodEX$TL <- ifelse (is.na(RodEX$TL), RodEX$HB+RodEX$T, RodEX$TL)

# get species averages
spMeansCR <- RodCR %>% group_by(sp) %>% summarise_each(funs(mean(.,na.rm=TRUE)))
spMeansEX <- RodEX %>% group_by(sp) %>% summarise_each(funs(mean(.,na.rm=TRUE))) %>% 
  select(-TL)

# put back together
ldSpecimens <- full_join(spMeansCR,spMeansEX)

# read Andrei M. data

amDATA <- read.csv("AMdata.csv",stringsAsFactors = FALSE) %>% 
          mutate(sp=gsub("_"," ",scname))%>% 
          select(-Wgrams,-sdWgrams,-nSpecimens,-extraSource,-scname) %>% 
          select(sp,everything())

# unratio craniodental measurements
for(i in 4:8){
  amDATA[,i] <- (amDATA[,i]/100)*amDATA$CBL
}
# unratio external measurements
for(i in 10:15){
  amDATA[,i] <- (amDATA[,i]/100)*amDATA$HB
}


# merge all the means
speciesMeans<- bind_rows(amDATA,ldSpecimens)

# read mass table
allMasses <- read.csv("bodySizes.csv",stringsAsFactors = F)
# merge mass and measurements
spMeansMasses <- inner_join(speciesMeans,allMasses)

# clean up
spMeansMasses <- spMeansMasses %>% select(-X) 
# rounding
spMeansMasses[,2:16] <- round(spMeansMasses[,2:16],2)

#write spMeans file

write.csv(spMeansMasses,file="meansMasses.csv")
