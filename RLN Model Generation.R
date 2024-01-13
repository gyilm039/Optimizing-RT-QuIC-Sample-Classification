#for questions, please contact gyilm039@uottawa.ca

#to download packages if not already installed
install.packages(setdiff(c("readxl", "xlsx", "ggplot2", "ROCR"), rownames(installed.packages())))  

#to load necessary packages
library("readxl")
library("xlsx")
library("ggplot2")
library("ROCR")

#to get the directory of your R code and input data
fileFolder <- dirname(rstudioapi::getSourceEditorContext()$path)
fileFolder <- gsub("/", "\\\\", fileFolder)
inputData <- paste(fileFolder, "RT-QuIC RLN Input Data.xlsx", sep="\\")
outputFile <- paste(fileFolder, "//RT-QuIC RLN Output Data.xlsx", sep="")

madeTabs <- excel_sheets(inputData)

#here, indicate how many RT-QuIC cycles you have used.
numberOfCycles <- 225

#here indicate your Ct cutoff, i.e. if a replicate doesn't cross the threshold, it has a Ct of 65 hours.
cutoffCycle <- 65


#if the ct of the raw fluoresence data is not already analysed, this section will analyse it and append to the input file
if(!('Analysed Data' %in% madeTabs)){
  myData <- read_excel(inputData,sheet = "Template")
  # colnames(myData) <- myData[1,]
  # myData <- myData[-c(1),]
  
  #remove clean ntc
  myData <- myData[!(myData$Tissue=="NA"),]
  
  #adjusted threshold to numeric
  myData$Threshold <- as.numeric(myData$Threshold)
  myData$Ct <- as.numeric(myData$Ct)
  myData$Concentration <- as.numeric(myData$Concentration)
  
  #determine the assay times
  assayTimes <- colnames(myData)
  assayTimes <- assayTimes[11:ncol(myData)]
  assayTimes <- as.numeric(assayTimes)
  
  #find Ct. If no Ct set to your cut-off
  for (currentRow in c(1:nrow(myData))){
    myData[currentRow, 10] <- cutoffCycle
    thresholdFound <- FALSE
    
    for (currentCycle in c(2:numberOfCycles)){
      if(!(is.na(myData[currentRow, currentCycle + 10]))){
        if(thresholdFound == FALSE){
          if(myData[currentRow, currentCycle + 10] > myData[currentRow, 9]){
            # if (y2-y1)/(x2-x1) = m, and yn = y1 + (xn-x1)/m
            # yn is threshold, then xn is ct
            slope <- (myData[currentRow, currentCycle + 10]-(myData[currentRow,  currentCycle + 9]))/(assayTimes[currentCycle]-assayTimes[currentCycle-1])
            b_value <- myData[currentRow, currentCycle + 10] - slope * assayTimes[currentCycle]
            myData[currentRow, 10] <- (myData[currentRow,9] - b_value)/slope
            thresholdFound <- TRUE
          } 
        }
      }
    }
  }
  
  #if ct > cutoff, change to cutoff
  myData$Ct[which(myData$Ct >cutoffCycle)] <- cutoffCycle
  
 # write.xlsx(as.data.frame(myData),file = myData, row.names=FALSE, sheetName = "Analysed Data", append=TRUE)
  write.xlsx(myData,inputData, sheetName = "Analysed Data", append = TRUE)
  
}else{
  myData <- read_excel(inputData,sheet = "Analysed Data")
  myData <- myData[,2:ncol(myData)]
  assayTimes <- colnames(myData[11:ncol(myData)])
}

################################################################################################################################################


#function that changes the format of scientific labels...can't rememember where I got this code from
scientific_10 <- function(x){
  x <-  gsub("e\\+*", "0^", scales::scientific_format()(x))
  parse(text=x)
}

#cleaning up the data labelling
myData[which(myData$IHC=="Positive"),8] <- "Pos"
myData[which(myData$IHC=="Negative"),8] <- "Neg"


##############################################################################################################
#this is a plot of RLN ct vs concentration. note that there will be a warning of missing points, this is because this figure only displays 10^-2 to 10^-8
#to include more dilutions, you can change the "scale_y_log10()" parameter

RLN <- ggplot(myData, aes(x = Ct, y = Concentration, colour = IHC)) +
  #geom_rect(aes(xmin=23.080398 , xmax=26.696104, ymin=0.002, ymax=0.000000005))+
  geom_point(size = 2)+
  scale_y_log10(breaks = c(0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001), expand = waiver(), label=scientific_10, limits = c(0.00000001, 0.01))+
  xlab("Cycle threshold (hours)")+
  ylab("Concentration (w/v)")+
  labs(colour="")+
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70), expand = waiver())+
  scale_colour_manual(values=c("brown2", "blue1"), breaks=c("Pos","Neg"))+
  theme_bw(base_size=12)+
  theme(axis.text.y=element_text(size = 12, hjust = 0.5), axis.text.x=element_text(size = 12), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("RLN Ct Distribution.png", dpi=1000, dev='png', height = 6, width = 6, units = "in")


###########################################################################################
#this section does stats to find out which dilutions are valid, i.e. positives have significantly lower ct than neg, p<0.05 by Mann-Whitney

dilutionPValues <- data.frame(matrix(nrow=10,ncol=4))
colnames(dilutionPValues) <- c("Concentration", "npos", "nneg", "p")

count <- 1
for(currentConcentration in c(10^-2, 10^-3,10^-4,10^-5,10^-6,10^-7,10^-8, 10^-9)){
  #below, this looks like a weird way of subsetting, it's just to avoid floating point errors
  currentSamples <- myData[(myData$Concentration>0.5*currentConcentration),]
  currentSamples <- currentSamples[(currentSamples$Concentration<1.5*currentConcentration),]
  
  dilutionPValues[count,1] <- currentConcentration
  dilutionPValues[count,2] <- length(which(currentSamples$IHC=="Pos"))

  dilutionPValues[count,3] <- length(which(currentSamples$IHC=="Neg"))

  dilutionPValues[count,4] <- wilcox.test(Ct ~ IHC, data=currentSamples)[3]
  count <- count + 1
}

write.xlsx(dilutionPValues, outputFile, sheetName = "Dilution p values", row.names = FALSE, append = TRUE)


#here, you can use the p values to guide your decision making about which dilutions to include in further analysis
#however, these p values are only a suggestion
#for example, we will do 10^-3 to 10^-8
  
#this will only get data from 10^-3 to 10^-8
includedSamples <- myData
includedSamples <- includedSamples[(includedSamples$Concentration<0.5*10^-2),]
includedSamples <- includedSamples[(includedSamples$Concentration>1.5*10^-9),]

#reworking data a little bit to have 0 and 1, for the sake of the logistic regression and roc curves
includedSamples[which(includedSamples$IHC=="Pos"),8] <- "1"
includedSamples[which(includedSamples$IHC=="Neg"),8] <- "0"
includedSamples$IHC <- as.numeric(includedSamples$IHC)

#this makes a logistic regression
logisticCurve <- glm(IHC ~ Ct, includedSamples, family=binomial)
includedSamples$Probability <- predict(logisticCurve, newdata = includedSamples, type = "response")


#this makes a logistic regression of ct vs positivity
classPredictions <- ggplot(includedSamples, aes(x = Ct, y = Probability)) +
  geom_point(size = 3)+
  geom_line()+
  xlab(expression("Cycle threshold (hours)"))+
  ylab("Probability of positivity")+
  xlim(0,65)+
  ylim(0,1)+
  labs(shape = "")+
  theme_bw(base_size=25)+
  theme(axis.text.y=element_text(size = 20, hjust = 0.5), axis.text.x=element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("RLN Prediction Graph.png", dpi=1000, dev='png', height = 8, width = 8, units = "in")

write.xlsx(as.data.frame(includedSamples), outputFile, sheetName = "Ct Regression Data", row.names = FALSE, append = TRUE)

################################################################################
#this next section does the roc analysis with ct
pred <- prediction(includedSamples$Ct, includedSamples$IHC, label.ordering = c(1,0))
perf <- performance(pred,"tpr","fpr")

aucROC <- performance(pred, measure = "auc")
aucCt <- aucROC@y.values[[1]]

#now let's find the threshold manually by calculating youden
## all the thresholds are in pred@cutoffs. assume a threshold of x, replicate is pos if ct < x

thresholdAnalysis <- data.frame(matrix(nrow=length(unlist(pred@cutoffs)), ncol=4))
colnames(thresholdAnalysis) <- c("Threshold", "Sens", "Spec", "Youden")

count <- 1
for(currentThreshold in unlist(pred@cutoffs)){
  thresholdAnalysis[count, 1] <- currentThreshold
  
  TP <- 0
  TN <- 0
  FP <- 0
  FN <- 0
  
  for(currentRow in c(1:nrow(thresholdAnalysis))){
    if(includedSamples[currentRow, 10] < currentThreshold){
      if(includedSamples[currentRow,8] ==1){                        ######## elisa pos
        TP <- TP + 1
      }else{                                                          ###eilsa neg
        FP <- FP + 1
      }
    }else{
      if(includedSamples[currentRow,8] ==1){
        FN <- FN + 1
      }else{
        TN <- TN + 1
      }
    }
  }
  
  thresholdAnalysis[count,1] <- currentThreshold
  thresholdAnalysis[count,2] <- TP/(TP+FN)
  thresholdAnalysis[count,3] <- TN/(TN+FP)
  thresholdAnalysis[count,4] <- thresholdAnalysis[count,2] + thresholdAnalysis[count,3] - 1
  
  count <- count + 1
}

#the optimal time here was defined as the ((longest time with highest youden) + (shortest time with highest youden))/2
optimizedCutoffs <- thresholdAnalysis[which(abs(thresholdAnalysis$Youden - maxCtYouden) < 0.005),]
  ctCutoff <- (max(optimizedCutoffs$Threshold)+min(optimizedCutoffs$Threshold))/2

write.xlsx(as.data.frame(thresholdAnalysis), outputFile, sheetName = "Ct Threshold Data", row.names = FALSE, append = TRUE)


################################################################################################################
#it's now time to do roc curves with MPR
#first, we will convert all the RFU data into MPR data

#now let's do mpr

MPRData <- data.frame(matrix(nrow=nrow(includedSamples), ncol=9+numberOfCycles-3))
MPRData[,1:9] <- includedSamples[,1:9]
colnames(MPRData) <- c(colnames(includedSamples[1:9]), paste(4:numberOfCycles))

#get initial value (4th reading)
MPRData[,10] <- 1

for(currentRow in c(1:nrow(MPRData))){
  
  baseline <- includedSamples[currentRow, 14]
  currentMPR <- 1
  
  for(currentCycle in c(5:numberOfCycles)){
    RFU <- includedSamples[currentRow,currentCycle + 10]
    
    if(RFU/baseline > currentMPR){
      currentMPR <- RFU/baseline
    }
    
    MPRData[currentRow, currentCycle + 6] <- currentMPR
  }
}


#now well do a bunch of roc curves over time and get tmpr and auc with youden index

analysedMPRData <- data.frame(matrix(nrow=numberOfCycles-4, ncol=5))
#analysedMPRData <- data.frame(matrix(nrow=10000, ncol=5))
colnames(analysedMPRData) <- c("Cycle", "Time", "Threshold", "AUC", "Youden")


for(currentROC in c(5:numberOfCycles)){
  
  newData <- cbind(MPRData[,1:8], MPRData[,currentROC+6])
  colnames(newData) <- c(colnames(MPRData)[1:8], "MPR")
  
  pred <- prediction(newData$MPR, newData$IHC)
  perf <- performance(pred,"tpr","fpr")
  
  count <- 1
  currentThresholdData <- data.frame(matrix(nrow=length(unlist(pred@cutoffs)),ncol=5))
  colnames(currentThresholdData) <- c("Cycle", "Threshold", "Sens", "Spec", "Youden")
  
  for(currentThreshold in unlist(pred@cutoffs)){
    
    
    TP <- 0
    TN <- 0
    FP <- 0
    FN <- 0
    
    for(currentRow in c(1:nrow(newData))){
      if(newData[currentRow, 9] > currentThreshold){
        if(newData[currentRow,8] ==1){                        ######## elisa pos
          TP <- TP + 1
        }else{                                                          ###eilsa neg
          FP <- FP + 1
        }
      }else{
        if(newData[currentRow,8] ==1){
          FN <- FN + 1
        }else{
          TN <- TN + 1
        }
      }
    }
    
    currentThresholdData[count,1] <- currentROC
    currentThresholdData[count,2] <- currentThreshold
    currentThresholdData[count,3] <- TP/(TP+FN)
    currentThresholdData[count,4] <- TN/(TN+FP)
    currentThresholdData[count,5] <- currentThresholdData[count,3] + currentThresholdData[count,4] - 1
    
    count <- count + 1
  }
  
  ##############################
  #now we have our all youdens, subset the highest youden values, then store the highest threshold (=highest spec)
  maxYouden <- max(currentThresholdData$Youden)
  idealThresholds <- currentThresholdData[which(abs(currentThresholdData$Youden - maxYouden) < 0.005),]     #did this cause of floating point precision errors
  mprcutoff <- ((max(idealThresholds$Threshold)+min(idealThresholds$Threshold))/2)
  
  if(currentROC==5){
    accumulatedTMPR <- 1
    mprcutoff <- 1
  }else{
    if(accumulatedTMPR > mprcutoff){
      mprcutoff <- accumulatedTMPR
    }else{
      accumulatedTMPR <- mprcutoff
    }
  }
  #now we need to get the median threshold and store values
  analysedMPRData[currentROC-4,1] <- currentROC
  analysedMPRData[currentROC-4,2] <- assayTimes[currentROC]
  analysedMPRData[currentROC-4,3] <- mprcutoff
  
  aucROCMPR <- performance(pred, measure = "auc")
  analysedMPRData[currentROC-4,4] <- aucROCMPR@y.values[[1]]
  analysedMPRData[currentROC-4,5] <- median(idealThresholds$Youden)
  
}  

MPRGraph <- ggplot(analysedMPRData, aes(x = Cycle, y = AUC)) +
  geom_point(size = 2)+
  geom_line()+
  ylab(expression("AUC"))+
  xlab("RT-QuIC cycle")+
  xlim(0,250)+
  ylim(0.4,1)+
  labs(shape = "")+
  theme_bw(base_size=25)+
  theme(axis.text.y=element_text(size = 20, hjust = 0.5), axis.text.x=element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank())


thresholdGraph <- ggplot(analysedMPRData, aes(x = Cycle, y = Threshold)) +
  geom_point(size = 2)+
  geom_line()+
  ylab(expression("T"[MPR]))+
  xlab("RT-QuIC cycle")+
  xlim(0,250)+
  ylim(1,4)+
  labs(shape = "")+
  theme_bw(base_size=25)+
  theme(axis.text.y=element_text(size = 20, hjust = 0.5), axis.text.x=element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank())


MPRGraph + thresholdGraph + plot_layout(ncol=1)
ggsave("RLN AUC over time, TMPR over time.png", dpi=1000, dev='png', height = 16, width = 8, units = "in")

###################################################################################
#now, we can do the dAUC/dT graph to help choose our assay duration with MPR

temp <- colnames(analysedMPRData)
temp <- c(temp, "Rate")

analysedMPRData$Time <- as.numeric(analysedMPRData$Time)

for (currentRow in c(1:nrow(analysedMPRData))){
  analysedMPRData[currentRow,6] <- (analysedMPRData[currentRow+1,4]-analysedMPRData[currentRow,4])/(analysedMPRData[currentRow+1,2]-analysedMPRData[currentRow,2])
}
colnames(analysedMPRData) <- temp

rateGraph <- ggplot(analysedMPRData, aes(x = Cycle, y = Rate)) +
  geom_point(size = 2)+
  geom_line()+
  ylab("dAUC/dt")+
  xlab("RT-QuIC cycle")+
  xlim(0,250)+
  ylim(-0.1, 0.1)+
  labs(shape = "")+
  theme_bw(base_size=25)+
  theme(axis.text.y=element_text(size = 20, hjust = 0.5), axis.text.x=element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("RLN dAUCdT.png", dpi=1000, dev='png', height = 8, width = 8, units = "in")

write.xlsx(as.data.frame(analysedMPRData), outputFile, sheetName = "MPR Threshold Data", row.names = FALSE, append = TRUE)
