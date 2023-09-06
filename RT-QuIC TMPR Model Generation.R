install.packages(setdiff(c("readxl", "xlsx", "ggplot2", "ROCR"), rownames(installed.packages())))  

library("readxl")
library("xlsx")
library("ggplot2")
library("ROCR")


fileFolder <- dirname(rstudioapi::getSourceEditorContext()$path)
fileFolder <- gsub("/", "\\\\", fileFolder)

inputData <- paste(fileFolder, "RT-QuIC Obex Input Data.xlsx", sep="\\")
outputFile <- paste(fileFolder, "//Obex TMPR Standards Output.xlsx", sep="")

madeTabs <- excel_sheets(inputData)

#here, indicate how many RT-QuIC cycles you have used.
numberOfCycles <- 300

#here is your assay duration
cutoffCycle <- 65

#if raw fluoresence data is not already analysed, this section will analyse it and append to the input file
if(!('Analysed Data' %in% madeTabs)){
  template <- read_excel(inputData,sheet = "Template")
  colnames(template) <- template[1,]
  template <- template[-c(1),]
  
  #remove clean ntc
  template <- template[!(template$Tissue=="NA"),]
  
  #adjusted threshold to numeric
  template$Threshold <- as.numeric(template$Threshold)
  template$Ct <- as.numeric(template$Ct)
  template$Concentration <- as.numeric(template$Concentration)
  
  #determine the assay times
  assayTimes <- colnames(template)
  assayTimes <- assayTimes[10:ncol(template)]
  assayTimes <- as.numeric(assayTimes)
  
  #find Ct. If no Ct set to your cut-off
  for (currentRow in c(1:nrow(template))){
    template[currentRow, 9] <- cutoffCycle
    thresholdFound <- FALSE
    
    for (currentCycle in c(2:numberOfCycles)){
      if(thresholdFound == FALSE){
        if(template[currentRow, currentCycle + 9] > template[currentRow, 8]){
          # if (y2-y1)/(x2-x1) = m, and yn = y1 + (xn-x1)/m
          # yn is threshold, then xn is ct
          slope <- (template[currentRow, currentCycle + 9]-(template[currentRow,  currentCycle + 8]))/(assayTimes[currentCycle-1]-assayTimes[currentCycle-2])
          b_value <- template[currentRow, currentCycle + 9] - slope * assayTimes[currentCycle-1]
          template[currentRow, 9] <- (template[currentRow,8] - b_value)/slope
          thresholdFound <- TRUE
        } 
      }
    }
  }
  
  template$Ct[which(template$Ct >cutoffCycle)] <- cutoffCycle
  
  write.xlsx(as.data.frame(template),file = inputData, row.names=FALSE, sheetName = "Analysed Data", append=TRUE)
}else{
  template <- read_excel(inputData,sheet = "Analysed Data")
  assayTimes <- colnames(template)
  assayTimes <- assayTimes[10:ncol(template)]
  assayTimes <- as.numeric(assayTimes)
}


##########################################################################################################
#let's try to determine which dilutions are valid
#here, we will pick an assay time, e.g. 28 hours and plot the distribution of MPRs

currentTime <- 28

MPRData <- data.frame(matrix(nrow=nrow(template), ncol=9))

#note background is defined here as the 4th reading
colnames(MPRData) <- c(colnames(template[,1:4]), "ELISA", "Hours", "Background", "Maximum", "MPR")

MPRData[,1:4] <- template[,1:4]
MPRData[,5] <- template[,7]
MPRData[,6] <- currentTime
MPRData[,7:8] <- template[,13]

#this code will determine all MPRs at the chosen assay time

for (currentRow in c(1:nrow(MPRData))){
  doneLoop <- FALSE
  for (currentCycle in c(4:length(assayTimes))){
    if(doneLoop == FALSE){
      if (assayTimes[currentCycle] > currentTime){
        #interpolate last RFU reading, this means that you find the mpr between this reading and the previous reading
        #y = mx + b, we have x, we need y which is the fluoresence reading
        slope <- (template[currentRow,currentCycle+9]-template[currentRow,currentCycle+8])/(assayTimes[currentCycle]-assayTimes[currentCycle-1])
        bValue <- template[currentRow,currentCycle+9]-slope*assayTimes[currentCycle]
        interpolatedRFU <- slope * currentTime + bValue
        
        if (interpolatedRFU > MPRData[currentRow, 7]){
          MPRData[currentRow, 8] <- interpolatedRFU
        }
        
        doneLoop <- TRUE
      }else if(template[currentRow, currentCycle + 9] > MPRData[currentRow, 7]){
        MPRData[currentRow, 8] <- template[currentRow, currentCycle + 9]
      }
    }
  }
}

MPRData[,9] <- MPRData[,8]/MPRData[,7]

scientific_10 <- function(x){
  x <-  gsub("e\\+*", "0^", scales::scientific_format()(x))
  parse(text=x)
}

ggplot(MPRData, aes(x = MPR, y = Concentration, colour = ELISA)) +
  geom_point(size = 3)+
  scale_y_log10(breaks=c(10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8, 10^-9, 10^-10, 10^-11))+
  labs(shape="")+
  scale_x_continuous(breaks = c(0,5,10,15,20,25), expand = waiver())+
  theme_bw(base_size=25)+
  theme(axis.text.y=element_text(size = 20, hjust = 0.5), axis.text.x=element_text(size = 20), panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#use whatever criteria you wish to remove dilutions. for our sake, we will use the 10^-3 to 10^-8 dilutions

#now we can drop dilutions that are not used
for (removeConcentration in c(10^-2, 10^-9, 10^-10, 10^-11)){
  template <- template[!(template$Concentration==removeConcentration),]
}
#################################################################################################################3
#now we can prepare a roc curve at each cycle

RFUDistributions <- data.frame(matrix(nrow=nrow(template), ncol=5+297))
colnames(RFUDistributions) <- c(colnames(MPRData)[1:5], assayTimes[4:300])

RFUDistributions[,1:4] <- template[,1:4]
RFUDistributions[,5] <- template[,7]
RFUDistributions[,6] <- template[,13]


for (currentRow in c(1:nrow(RFUDistributions))){
  for (currentCycle in 5:300){
    if (template[currentRow, currentCycle+9] > RFUDistributions[currentRow, currentCycle+1]){
      RFUDistributions[currentRow, currentCycle+2] <- template[currentRow, currentCycle+9]
    }else{
      RFUDistributions[currentRow, currentCycle+2] <- RFUDistributions[currentRow, currentCycle+1]
    }
  }
}

#now we will subset the mpr from each cycle then prepare a roc curve
durationAUC <- data.frame(matrix(nrow = 226, ncol = 5))
colnames(durationAUC) <- c("Time", "Cycle", "AUC", "Tissue", "TMPR")

currentData <- data.frame(matrix(nrow = 132, ncol = 7))
colnames(currentData) <- c(colnames(MPRData)[1:5], "Cycle", "MPR")

MPRDistributions <- data.frame(matrix(nrow=nrow(RFUDistributions), ncol=ncol(RFUDistributions)))
colnames(MPRDistributions) <- colnames(RFUDistributions)
MPRDistributions[1:5] <- RFUDistributions[1:5]

predictionData <- data.frame(matrix(nrow = nrow(template), ncol = 4))
colnames(predictionData) <- c("Concentration", "ELISA", "MPR", "Prediction")



for (rocCycle in c(4:226)){
  
  MPRDistributions[,rocCycle+2] <- RFUDistributions[,rocCycle+2]/RFUDistributions[,6]
  
  currentData[,1:5] <- MPRDistributions[,1:5]
  currentData[,6] <- rocCycle
  currentData[,7] <- MPRDistributions[,rocCycle+2]
  
  #class prediction model
  pred <- prediction(currentData$MPR, currentData$ELISA)
  
  #roc curve
  perf <- performance(pred,"tpr","fpr")
  
  integral <- performance(pred, measure = "auc")
  
  durationAUC[rocCycle-3, 1] <- assayTimes[rocCycle]
  durationAUC[rocCycle-3, 2] <- rocCycle
  durationAUC[rocCycle-3, 3] <- integral@y.values[[1]]
  durationAUC[rocCycle-3, 4] <- "Obex"
  
  predictionData$Concentration <- MPRDistributions$Concentration
  predictionData$ELISA <- MPRDistributions$ELISA
  predictionData$ELISA[which(predictionData$ELISA=="Negative")] <- 0
  predictionData$ELISA[which(predictionData$ELISA=="Positive")] <- 1
  predictionData$ELISA <- as.numeric(predictionData$ELISA)
  predictionData$MPR <- MPRDistributions[, rocCycle+2]
  
  logisticCurve <- glm(ELISA ~ MPR, predictionData, family=binomial)
  predictionData$Prediction <- predict(logisticCurve, newdata = predictionData, type = "response")
  
  ROCCurveData <- data.frame(matrix(nrow=length(unlist(perf@x.values)), ncol = 2))
  colnames(ROCCurveData) <- c("xValues", "yValues")
  ROCCurveData$xValues <- unlist(perf@x.values)
  ROCCurveData$yValues <- unlist(perf@y.values)
  
  
  sortedAFR <- unique(sort(predictionData$AFR, decreasing = TRUE))
  idealDistance <- 0
  idealDataPoint <- 0
  allDistances <- data.frame(matrix(nrow=length(ROCCurveData$xValues), ncol = 2))
  colnames(allDistances) <- c("Point", "Distance")
  
  #determine each points distance from the line y=x
  for (currentPoint in c(1:length(ROCCurveData$xValues))){
    currentDistance <- abs(ROCCurveData[currentPoint,1]-ROCCurveData[currentPoint,2])/sqrt(2)
    if (currentDistance > idealDistance){
      idealDistance <- currentDistance
      idealDataPoint <- currentPoint
    }
  }
  
  TMPR <- sort(predictionData$MPR)[idealDataPoint]

  if(idealDataPoint==0){
    TMPR <- NA
  }
  durationAUC[rocCycle-3, 5] <- TMPR

  
}

#an assay duration with the highest auc can be chosen as the ideal assay duration, and the calculated TMPR can be used.
write.xlsx(durationAUC,file = outputFile, row.names=TRUE, sheetName = "Results", append=TRUE)


