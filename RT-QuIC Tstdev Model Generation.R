install.packages(setdiff(c("readxl", "xlsx", "ggplot2", "ROCR"), rownames(installed.packages())))  

library("readxl")
library("xlsx")
library("ggplot2")
library("ROCR")


fileFolder <- dirname(rstudioapi::getSourceEditorContext()$path)
fileFolder <- gsub("/", "\\\\", fileFolder)

inputData <- paste(fileFolder, "RT-QuIC Obex Input Data.xlsx", sep="\\")
outputFile <- paste(fileFolder, "//Obex Tstdev Standards Output.xlsx", sep="")

madeTabs <- excel_sheets(inputData)

#here, indicate how many RT-QuIC cycles you have used.
numberOfCycles <- 300

#here indicate your Ct cutoff, i.e. if a replicate doesn't cross the threshold, it has a Ct of 65 hours.
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
  assayTimes <- assayTimes[11:ncol(template)]
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
}
  
#-------------------------------------------------------------------------------------------------------------------------------#

#. This next section will do compare the positive and negatives of each dilution, to see if the dilution is valid
results <- data.frame(nrow = 10, ncol = 2)
colnames(results) <- c("Concentration", "P-value")
count <- 1
#may encounter warning - cannot compute exact-pvalue with ties - this is fine.
for (concentration in unique(template$Concentration)){
  testedSamples <- template[template$Concentration==concentration, 1:9]
  #testedSamples[,11] <- 1/testedSamples[,9]

  wilcoxResult <- wilcox.test(unlist(testedSamples$Ct~testedSamples$ELISA, testedSamples))
  if(is.nan(unlist(wilcoxResult[3]))){
    wilcoxResult[3]<-1
  }
  results[count, 1] <- concentration
  results[count, 2] <- wilcoxResult[3]
  count <- count + 1
}
write.xlsx(as.data.frame(results),file = outputFile, row.names=FALSE, sheetName = "Dilution p-values", append=TRUE)




#now we can drop dilutions that are not used
for (removeConcentration in c(10^-2, 10^-9, 10^-10, 10^-11)){
  template <- template[!(template$Concentration==removeConcentration),]
}

#this data isn't strictly needed for analysis, but can be usual for visual
predictionData <- data.frame(matrix(nrow = nrow(template), ncol = 4))
colnames(predictionData) <- c("Concentration", "ELISA", "Ct", "PredictionCt")
predictionData$Concentration <- template$Concentration
predictionData$ELISA <- template$ELISA
predictionData$Ct <- template$Ct
predictionData$ELISA[which(predictionData$ELISA=="Negative")] <- 0
predictionData$ELISA[which(predictionData$ELISA=="Positive")] <- 1
predictionData$ELISA <- as.numeric(predictionData$ELISA)

logisticCurve <- glm(ELISA ~ Ct, predictionData, family=binomial)
predictionData$PredictionCt <- predict(logisticCurve, newdata = predictionData, type = "response")




#class prediction model
pred <- prediction(1/predictionData$Ct, predictionData$ELISA)

#roc curve
perf <- performance(pred,"tpr","fpr")


#roc curve, AUC
integral <- performance(pred, measure = "auc")
AUC <- integral@y.values[[1]]
ROCCurveData <- data.frame(matrix(nrow=length(unlist(perf@x.values)), ncol = 2))
colnames(ROCCurveData) <- c("xValues", "yValues")
ROCCurveData$xValues <- unlist(perf@x.values)
ROCCurveData$yValues <- unlist(perf@y.values)


write.xlsx(as.data.frame(predictionData),file = outputFile, row.names=FALSE, sheetName = "Regression Data", append=TRUE)
write.xlsx(ROCCurveData,file = outputFile, row.names=FALSE, sheetName = "ROC Points", append=TRUE)


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
  allDistances[currentPoint, 1] <- currentPoint
  allDistances[currentPoint, 2] <- currentDistance
}

#the ideal time is then the time range corresponding to the furthest distance from y =x to the next point
lowerTimeRange <- sort(predictionData$Ct)[idealDataPoint]
higherTimeRange <- sort(predictionData$Ct)[idealDataPoint+1]

timeRange <- data.frame(matrix(ncol=1, nrow=2))
timeRange[,1] <- c(lowerTimeRange, higherTimeRange)
rownames(timeRange) <- c("Lower estimate of ideal assay duration", "Higher estimate of ideal assay duration")
write.xlsx(timeRange,file = outputFile, row.names=TRUE, sheetName = "Results", append=TRUE)


