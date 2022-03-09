####Create boxplots of various spectral indices for the different cover crop types. perform ANOVA and calculate TukeyHSD values


library(multcompView)


#set temp dir to slate folder
#install.packages('unixtools',repos='http://www.rforge.net/')
unixtools::set.tempdir("/N/slate/farellam")
tempdir()

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

#look at CCtype
cc <- subset(RFdata, RF2class=="CC")
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

co_counts <- (cc) %>%
        group_by(site) %>%
        summarize(Brassica = sum(cc.agg == "Brassica"),
                  CerealRye = sum(cc.agg == "CerealRye"),
                  RyeGrass = sum(cc.agg == "RyeGrass"),
                  Wheat = sum(cc.agg == "Wheat"),
                  WinterGrain = sum(cc.agg == "WinterGrain"),
                  Mix = sum(cc.agg == "Mix/Other"))
co_counts


#select only the most popularr CC types and then combine with NoTill data for plotting
#plotdata <- subset(cc, cc.agg %in% c("CerealRye", "RyeGrass", "Wheat", "winterGrain"))
notill <- subset(RFdata, RF3class =="NoTill")
notill$cc.agg <- "NoTill"
plotdata <- rbind.data.frame(cc, notill)

#order factor levels
plotdata$cc.agg <- ordered(plotdata$cc.agg, levels = c("Brassica", "CerealRye", "RyeGrass", "Wheat", "WinterGrain", "Mix/Other", "NoTill"))


#####BOX PLOTS####

########med NDVI###########
aov1 <- aov(medNDTI~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_medNDV.tif", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches

boxplot(plotdata$medNDVI ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(0,1.1),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="median NDVI Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")
text(1, 0.88, labels="bc")
text(2, 0.87, labels="d")
text(3, 0.95, labels="cd")
text(4, 0.93, labels="a")
text(5, 0.90, labels="a")
text(6, 0.85, labels="b")
text(7, 0.95, labels="e")

text(1.25, 1.07, labels=expression('F = 220.4; p = 2 x 10'^-16))
dev.off()


########med LST###########
aov1 <- aov(medLST~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_medLST", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches
boxplot(plotdata$medLST ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(-15, 20),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="median LST (°C) Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")
text(1, 18, labels="a")
text(2, 19, labels="b")
text(3, 14.5, labels="b")
text(4, 15, labels="a")
text(5, 18, labels="a")
text(6, 15, labels="b")
text(7, 20, labels="b")

text(1.25, -14.5, labels=expression('F = 32.56; p = 2 x 10'^-16))
dev.off()


########sd LST###########
aov1 <- aov(sdLST~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_sdLST", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches
boxplot(plotdata$sdLST ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(0, 20),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="LST sd from Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")
text(1, 15.75, labels="a")
text(2, 15, labels="c")
text(3, 15.75, labels="c")
text(4, 14.25, labels="a")
text(5, 16.25, labels="ab")
text(6, 16.25, labels="a")
text(7, 17.25, labels="b")

text(1.25, 20, labels=expression('F = 15.09; p = 2 x 10'^-16))
dev.off()

########med NDTI###########
aov1 <- aov(medNDTI~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_medNDTI", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches
boxplot(plotdata$medNDTI ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(-0.1, 0.4),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="median NDTI Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")
#text(1, 0.36, labels="ab")
#text(2, 0.33, labels="b")
#text(3, 0.34, labels="ab")
#text(4, 0.37, labels="a")
#text(5, 0.38, labels="ab")
#text(6, 0.35, labels="ab")
#text(7, 0.35, labels="b")

text(1.25, 0.4, labels=expression('F = 0.90; p = 0.496'))
dev.off()

########med STI###########
aov1 <- aov(medSTI~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_medSTI", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches
boxplot(plotdata$medSTI ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(0, 2.25),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="median STI Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")

text(1, 2.1, labels="a")
text(2, 2.05, labels="a")
text(3, 2.05, labels="b")
text(4, 2.2, labels="a")
text(5, 2.2, labels="a")
text(6, 2.1, labels="a")
text(7, 2.1, labels="a")

text(1.3, 0, labels=expression('F = 11.06; p = 2.52 x 10'^-12))
dev.off()

########med ExG###########
aov1 <- aov(medExG~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_medExG", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches
boxplot(plotdata$medExG ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(-0.2,0.5),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="median ExG Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")
text(1, 0.44, labels="ab")
text(2, 0.42, labels="ab")
text(3, 0.47, labels="ab")
text(4, 0.44, labels="a")
text(5, 0.4, labels="ab")
text(6, 0.43, labels="ab")
text(7, 0.42, labels="b")

text(1.25, 0.5, labels=expression('F = 3.13; p = 4.62 x 10'^-3))
dev.off()

########med ExGR###########
aov1 <- aov(medExGR~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_medExGR", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches
boxplot(plotdata$medExGR ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(-1.18,-0.4),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="median ExGR Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")
text(1, -0.45, labels="bcd")
text(2, -0.5, labels="d")
text(3, -0.43, labels="cd")
text(4, -0.48, labels="ab")
text(5, -0.5, labels="a")
text(6, -0.50, labels="c")
text(7, -0.40, labels="e")

text(1.25, -1.18, labels=expression('F = 59.83; p = <2 x 10'^-16))
dev.off()

########med STI###########
aov1 <- aov(medSTI~cc.agg, plotdata)
summary(aov1)
print(model.tables(aov1, "means", se=T))
tuk1 <- TukeyHSD(aov1)
multcompLetters4(aov1, tuk1)

tiff("/N/project/CoverCrops/CoverCrops/plots/CC_medSTI.tif", units='in', width=9, height=5, res=300)
par(mar=c(2.2,3.4,0.1,0.1),  family="Calibri Light") #2x2 plot 'mai' set the plot margins in inches, 'oma' sets the outer margins in inches

boxplot(plotdata$medSTI ~ plotdata$cc.agg, ylab="", xlab="", ylim=c(0.6,2.25),
        col=(c("#AFC2D5", "#CCDDD3", "#DFEFCA", "#FFF9A5", "#B48B7D", "#EF959D", "#69585F")))
title(ylab="median STI Nov 1 - Mar 31", line=2.3, cex.lab=1, family="Calibri Light")
text(1, 2.05, labels="bc")
text(2, 0.87, labels="d")
text(3, 0.95, labels="cd")
text(4, 0.93, labels="a")
text(5, 0.90, labels="a")
text(6, 0.85, labels="b")
text(7, 0.95, labels="e")

text(1.25, 1.07, labels=expression('F = 220.4; p = 2 x 10'^-16))
dev.off()
