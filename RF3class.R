####Build the three-class (CC, no till, and conventional till) RF model. 


library(dplyr)
library(caTools)
library(caret)
library(doParallel)
library(multcompView)

#set temp dir to slate folder
#install.packages('unixtools',repos='http://www.rforge.net/')
unixtools::set.tempdir("/N/slate/farellam")
tempdir()

#"Clean" Start to random forest model runs
#pacman::p_load(raster, gdalUtils, terra, dplyr, rgdal, sp, rhdf5, rlist, randomForest, caTools, caret, maptools, e1071, caretEnsemble, lattice, gridExtra,scales)

#load RFdata
RFdata <- read.csv("/N/project/CoverCrops/CoverCrops/data/RFdata/cloudmsk_RFdata_25mbuff.csv")
summary(as.factor(RFdata$site))
names(RFdata)

#RFdata <- RFdata[,-c(19:20)]

#get rid of NA and Inf vals
RFdata[sapply(RFdata, is.infinite)] <- NA
RFdata <- RFdata[rowSums(is.na(RFdata[,17:51]))!=35,]
  

#make new variables for the 2 and 3 class models -- these need to be in ABC order with Cover Crop being 1st
#two-class model: cover crop presence/absence 

RFdata$RF2class <- "CC"
RFdata$RF2class[RFdata$Cover_Crop == "N"] <- "NoCC" 
RFdata$RF2class <- as.factor(RFdata$RF2class)

#three-class model: cover crop specific category 0 = nothing (presumably conventional tillage), 1= 'No-Till', 2 = "Cover crop".  
#RFdata$RF3class[RFdata$Cover_Crop != "N" | RFdata$Tillage != "n" | RFdata$Tillage != "C" | RFdata$Tillage == "c"] <- "ZConv"
#RFdata$RF3class[RFdata$Tillage == "N" | RFdata$Tillage == "n"] <- "NoTill"
#RFdata$RF3class[RFdata$Cover_Crop == "A" | RFdata$Cover_Crop=="B" | RFdata$Cover_Crop=="BC"| RFdata$Cover_Crop=="C"| RFdata$Cover_Crop=="W" | RFdata$Cover_Crop=="G"] <- "CC"
#RFdata$RF3class <- as.factor(RFdata$RF3class)
#levels(as.factor(RFdata$RF3class))

#make empty column to hold 3class classes
RFdata$RF3class <- NA
#determine if fields had No Till 
RFdata$RF3class[RFdata$Tillage == 'N' | RFdata$Tillage == 'n'] <- 'NoTill'
#determine if fields had Conventional Tillage
RFdata$RF3class[RFdata$Tillage == 'C' | RFdata$Tillage == 'c'] <- 'ZConv'
#if a field had CC in the 2 class model it should also have CC in the 3 class model
RFdata$RF3class[RFdata$RF2class == 'CC'] <- 'CC'
RFdata$RF3class <- as.factor(RFdata$RF3class)
levels(as.factor(RFdata$RF3class))


#check the classes. 'NoTill' and 'ZConv' in the 3-class model should add up to 'NoCC' in the 2-class model  
plyr::count(RFdata$RF2class)
plyr::count(RFdata$RF3class)


#look at the counts x county
co_counts <- (RFdata) %>%
  group_by(site) %>%
  summarize(total = length(RF2class),
            YesCC_2class = sum(RF2class=="CC"),
            NoCC_2class = sum(RF2class=="NoCC"),
            Per_CC = ((sum(RF2class=="CC")) / (length(RF2class))) * 100,
            Conventional_3class = sum(RF3class == "ZConv", na.rm=TRUE),
            NoTill_3class = sum(RF3class == "NoTill", na.rm = TRUE),
            CC_3class = sum(RF3class == "CC", na.rm = TRUE))
co_counts
#2class total CoverCrops x county for 2015:
#Benton 9/254
#Gibson 64/272
#Posey 56/197
#Warren 16/227
#White 10/249
#if wanting to only analyze sites Mallory used


RFdata$sitedate <- paste(RFdata$site, RFdata$date, sep="")                   

#order data columns for RF
names(RFdata)
RFdata <- RFdata[,-c(1,3,4,5:8)]

#######**OMIT NA CLASS VALUES**#######
RFdata <- RFdata[!is.na(RFdata$RF3class),]#only do this for the 3class model!
######################################
####RANDOM FOREST ML####
set.seed(500)

#to randomly select samples from the whole data set
sample <- sample.split(RFdata$RF3class, SplitRatio = 0.8)
RFtrain <- subset(RFdata, sample==TRUE)
RFtest <- subset(RFdata, sample==FALSE)
print(paste(nrow(RFtrain), " samples in the training data", sep=""))
print(paste(nrow(RFtest), " samples in the testing data data", sep=""))


#upsample the training data to have equal classes
#RFtrain <- upSample(RFtrain, RFtrain$RF3class)

# Set up a resampling method in the model training process
tc <- caret::trainControl(method = "repeatedcv", # repeated cross-validation of the training data
                          number = 10, # number of folds
                          repeats = 10, # number of repeats
                          allowParallel = TRUE, # allow use of multiple cores if specified in training
                          verboseIter = TRUE,# view the training iterations
                          classProbs=TRUE, # calculate class probabilities
                          sampling='up') #to upsample the training data to have equal number of observations in each class

# Generate a grid search of candidate hyper-parameter values for inclusion into the models training
# These hyper-parameter values are examples. You will need a more complex tuning process to achieve high accuracies
# For example, you can play around with the parameters to see which combinations gives you the highest accuracy. 
rf.grid <- expand.grid(mtry=1:35) # number of variables available for splitting at each tree node

## Begin training the models. Took ~1hr for windshield all years on Carbonate using 5 cores (Fall 2015 ~15 min using 1 core)
# Train the random forest model
names(RFtrain)
cl <- makePSOCKcluster(12)
registerDoParallel(cl)

Sys.time()
rfmodel <- caret::train(x = RFtrain[,(10:44)], y = as.factor(RFtrain$RF3class),
                            method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
Sys.time()
rfmodel$finalModel
varImp(rfmodel)
stopCluster(cl)


saveRDS(rf_modelLST, "/N/project/CoverCrops/CoverCrops/RFmodels/RF3class_HLS_cloudmsk_25mbuff.rds")

rfmodel <- readRDS("/N/project/CoverCrops/CoverCrops/RFmodels/RF3class_HLS_cloudmsk_25mbuff.rds")



# compare predicted outcome and true outcome
RF3class_pred <- predict(rfmodel, RFtest[,10:44])
confusionMatrix(RF3class_pred, as.factor(RFtest$RF3class))


#looking at NDVI of misclassified obs
RFtest$pred <- predict(rfmodel, RFtest[,10:44])
cc <- subset(RFtest, RF3class=="CC")
cc_wrong <- subset(cc, pred %in% c("NoTill", "ZConv"))
cc_right <- subset(cc, pred=="CC")

summary(cc_wrong$medNDVI)
summary(cc_right$medNDVI)
hist(cc_wrong$medNDVI)
hist(cc_right$medNDVI)

###############################Testing on INDV Counties##############
#benton <- subset(all, site=="Benton")
#decatur <- subset(all, site=="Decatur")
#gibson <- subset(all, site=="Gibson")
#hamilton <- subset(all, site=="Hamilton")
#posey <- subset(all, site=="Posey")
#stjoe <- subset(all, site=="StJoe")
#warren <- subset(all, site=="Warren")
#washington <- subset(all, site=="Washington")
#white <- subset(all, site=="White")
#whitley <- subset(all, site=="Whitley")

#counties <- list(benton, decatur, gibson, hamilton, posey, stjoe, warren, washington, white, whitley)

#results <- matrix(NA, nrow=10, ncol=3)

#j = 0
#for(i in counties){
#  RF3class_pred <- predict(rf_modelLST, i[,10:44])
#  conf <- confusionMatrix(RF3class_pred, as.factor(i$RF3class))
#  kappa <- conf$overall[[2]]
#  ccacc <- conf$byClass[[1,1]]
#  results[j+1,1] <- unique(i$site)
#  results[j+1,2] <- round(kappa, digits=3)
#  results[j+1,3] <- round(ccacc, digits=3)
#  j <- j+1
#}
 
#colnames(results) <- c("county", "kappa", "CCacc")
#view(results)
#print(results)
