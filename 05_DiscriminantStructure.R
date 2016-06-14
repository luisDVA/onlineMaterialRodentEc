######### Discriminant structure

# load packages
library(dplyr)
library(mda)


# Read diet data
ecolData <- read.csv("DietLocomotion.csv",stringsAsFactors = FALSE)%>% 
  select(sp,DietCond,Locomotion)
ecolData$sp <- gsub(" ","_",ecolData$sp)

# Reading in the morphology data and assigning rownames
morphCorrected <- read.csv("sizeCorrected.csv",stringsAsFactors = FALSE) %>% select(-X)
rownames(morphCorrected) <- morphCorrected$sp

#merge all (sort=F keeps the table in the right order)
ddA <- merge(morphCorrected,ecolData,sort=FALSE)

####### training Data for corr structure
# DIET
ddDiet <- ddA %>% filter(DietCond!="U") %>% select(-Locomotion)
row.names(ddDiet) <- ddDiet$sp
ddDiet$sp <- NULL
# Locomotion
ddLoc <- ddA %>% filter(Locomotion!="U") %>% select(-DietCond)
row.names(ddLoc) <- ddLoc$sp
ddLoc$sp <- NULL
#measurements only
ddDietMeas <- ddDiet[,1:12]
ddLocMeas <- ddLoc %>% select(-Locomotion)


######### run fda on training data, all Chars, DIET
# fit FDA
dietFDA <- fda(DietCond~.,data=ddDiet)
# predict variates for training data
dietFDAvariates <- predict(dietFDA,type="variates")
# get rid of rownames that come from the original diet factor
rownames(dietFDAvariates) <- NULL


# for correlation tests
library(psych)


#bind and run test
corsDiet <- cbind(dietFDAvariates,ddDietMeas) %>% corr.test(adjust="holm")

# print values with the holm adjustment (values above the diagonal)
write.csv(round(corsDiet$r,2)[4:15,1:3],quote = F)
write.csv(t(round(corsDiet$p,2)[1:3,4:15]),quote=F)

# discriminant power and verify correlations 
library(DiscriMiner)


# Discriminant power
discPowerDiet <- round(discPower(ddDiet[,1:12],ddDiet$DietCond),3)
discPowerDiet$char <- row.names(discPowerDiet)
# get Discriminant power, sorted by Wilks lambda and print
discPowerDiet %>% select(char,everything()) %>% arrange(wilks_lambda) %>% 
  write.csv(quote=F)

####################LOCOMOTION########################
#### fda on training data, all characters, Locomotion

locFDA <- fda(factor(Locomotion)~.,data=ddLoc)
locFDAvariates <- predict(locFDA,type="variates")
row.names(locFDAvariates) <- NULL

## correlation test

corsLoc <- cbind(locFDAvariates,ddLocMeas) %>% corr.test(adjust="holm")

# for pasiting into word, values above diagonal
write.csv(round(corsLoc$r,2)[7:18,1:6],quote = F)
write.csv(t(round(corsLoc$p,2)[1:6,7:18]),quote=F)


# discriminant power, calculated without fossorial species because n=1
ddLocgop <- ddLoc
ddLocgop <-  rbind(ddLoc,ddLocgop[which(rownames(ddLocgop)=="Cratogeomys_fumosus"),])
# calculate discPower
dpowerLoc <- round(discPower(ddLocgop[,1:12],ddLocgop$Locomotion),3)
# sort and print
dpowerLoc %>% mutate(char=rownames(dpowerLoc)) %>%
  select(char,everything()) %>%   arrange(wilks_lambda) %>% write.csv(quote=F)



