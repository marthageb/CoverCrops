####Build the two-class (conventional tillage and cover crop) RF model. Apply 2-class model to No Till and other tillage observations (fig. 3) and calculate misclassification rastes for the different CC types (Table XXXX)####


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
#assign the other tillage practices as othertill
RFdata$RF3class[RFdata$Tillage == 'E' | RFdata$Tillage == 'e' | RFdata$Tillage == 'm' | RFdata$Tillage == 'M' |
                  RFdata$Tillage == 'R' | RFdata$Tillage == 'S' | RFdata$Tillage == 'U' | RFdata$Tillage == "u"] <- 'otherTill'
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
            CC_3class = sum(RF3class == "CC", na.rm = TRUE),
            OtherTill = sum(RF3class == "otherTill", na.rm=TRUE))
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
#RFdata <- RFdata[!is.na(RFdata$RF3class),]#only do this for the 3class model!
######################################

all <- RFdata
RFdata <- subset(all, RF3class %in% c("ZConv", "CC"))
levels(RFdata$RF3class)
RFdata$RF3class <- (droplevels(RFdata$RF3class))



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
RF2class <- caret::train(x = RFtrain[,(10:44)], y = as.factor(RFtrain$RF3class),
                            method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
Sys.time()
RF2class$finalModel
varImp(RF2class)
stopCluster(cl)

saveRDS(RF2class, "/N/project/CoverCrops/CoverCrops/RFmodels/RF2class_Conv_CC.rds")
      

# compare predicted outcome and true outcome
RF3class_pred <- predict(RF2class, RFtest[,10:44])
confusionMatrix(RF3class_pred, as.factor(RFtest$RF3class))

notill <- subset(all, RF3class == "NoTill")
other <- subset(all, RF3class == "otherTill")

test <- rbind.data.frame(notill, other, RFtest)
test$predicted <- predict(RF2class, test[,10:44])


########################Determining CC type Misclassification###########
#look at CCtype
cc <- subset(RFtest, RF3class=="CC")
summary(as.factor(cc$Cover_Crop))

#define CC type
cc$cc.agg <- cc$Cover_Crop
cc$cc.agg[cc$Cover_Crop %in% c("0", "0B", "5N", "AB", "AL", "ann", "AP", "ba", "BC", "BG", "BLL", "BLS",
                               "bo", "BO", "BP", "BW", "cb", "CL", "CP", "CR", "GP", "H", "l", "L", "LO",
                               "O", "OB", "ob", "OB", "OBS", "OR", "Q", "R", "U", "X", "Y", "P")] <- 'Mix/Other'
cc$cc.agg[cc$Cover_Crop %in% c("a", "A")] <- "RyeGrass"
cc$cc.agg[cc$Cover_Crop == "B"] <- "Brassica"
cc$cc.agg[cc$Cover_Crop %in% c("c", "C")] <- "CerealRye"
cc$cc.agg[cc$Cover_Crop == "G"] <- "WinterGrain"
#cc$cc.agg[cc$Cover_Crop == "P"] <- "Mix"
cc$cc.agg[cc$Cover_Crop %in% c("w", "W")] <- "Wheat"
summary(as.factor(cc$cc.agg))

# compare predicted outcome and true outcome
RF2class_pred <- predict(RF2class, cc[,10:44])
cc$pred <- RF2class_pred 
confusionMatrix(RF2class_pred, as.factor(cc$RF3class))

wrong <- subset(cc, pred == "ZConv")
wrong.cc <- as.data.frame(summary(as.factor(wrong$cc.agg)))
names(wrong.cc) <- 'incorrect'
wrong.cc$type <- rownames(wrong.cc)
summary(as.factor(cc$cc.agg))
total.cc <- as.data.frame(summary(as.factor(cc$cc.agg)))
names(total.cc) <- 'total'
total.cc$type <- rownames(total.cc)

final <- merge(wrong.cc, total.cc, by="type")
final$misclass <- final$incorrect / final$total * 100

