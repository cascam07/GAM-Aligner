library(magrittr)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(gam)
library(rvest)
library(xml2)

setwd("E://PVM_M//Alignments")

#Read list of anchor targets
anchors <- read.csv("Pos_anchors_GAM.csv", stringsAsFactors = F) #rows = files, cols = anchors
baseline <- "PVM_M_B_EZ_4_10Hr_L_Pos_90_30Jan17_Polaroid_HSST-A923" #file to align others to
mzXMLdocs <- list.files(path="E://PVM_M//MS_Files//mzXML//Pos",pattern="*.mzXML", full.names = T)
outputpath <- "E://PVM_M//MS_Files//mzXML/"
targetFiles <- list.files(path="E://PVM_M//Pos_Targets",pattern="*.csv", full.names = T)

rtResiduals <- GenResidualsFromAnchors(anchors, baseline)

#Visualize residuals vs basline retention time in lattice plots
rtResiduals %>% ggplot(aes(x = baseline_RT, y = diff_RT, color = Dataset)) + geom_point() + scale_color_discrete(guide = FALSE) + stat_smooth(colour = 'black') + facet_wrap(~ Dataset)

#split run into 3 time regions: 0-30 min, 30-50 min, 50-90 min
time_region <- function(x){ifelse(x <= 30, 1, ifelse(x>=50, 3, 2))}
#time_region <- function(x){ifelse(x <= 15, 1, ifelse(x>15 & x<=25, 2, ifelse(x>25 & x<=50, 3, ifelse(x>50 & x<=60, 4, 5))))}
#time_region <- function(x){ifelse(x >= 30, 1,2)}
rtResiduals %<>% mutate(timeRegion = time_region(baseline_RT)) %>% mutate(timeRegion = as.factor(timeRegion))

#model the retention time drifts to correct, df = 16
models <- rtResiduals %>% split(.$Dataset) %>% map(function(x) {
  gam(diff_RT ~ s(baseline_RT, 4) + timeRegion + timeRegion * s(baseline_RT, 4), data = x)
})

models <- rtResiduals %>% split(.$Dataset) %>% map(function(x) {
  gam(diff_RT ~ s(baseline_RT, 2), data = x)
})

AlignFiles(mzXMLdocs, models, outputpath)
CorrectTargetFiles(targetFiles, models, outputpath)




#################################
##### Large Functions Below #####
#################################

#Align the mzXML files using the provided models and saves aligned files to the output directory
AlignFiles <- function(mzXMLdocs, models, outputpath){

  for (d in mzXMLdocs) {
    file_root <- d %>% basename %>% sub(".mzXML", "",.)
    doc <- read_xml(d)
    print(paste(file_root,"    ", d))
    
    #Extract Scan Elements and retention times from XML
    scanNodes <- xml_children(doc)[[1]] %>% xml_children() %>% .[xml_name(.) == "scan"]
    
    #Cleaning and extraction 
    retentionTimes <- xml_attr(scanNodes, 'retentionTime') %>% gsub('^PT', '', .) %>% gsub('S$', '', .) %>% 
      as.numeric() %>% { . / 60 } #  minute/second conversion 
  
    #Creating a dataframe of values to predict on
    predDat <- data.frame(retentionTimes, time_region(retentionTimes)) %>% setNames(c('baseline_RT', 'timeRegion'))
    predDat$timeRegion %<>% factor(levels = levels(rtResiduals$timeRegion))
    
    model <- models[[file_root]]
    if(is.null(model)) { #Baseline file will have a null model. Save it an move on.
      write_xml(x=doc, file=paste(outputpath, file_root, "_Baseline.mzXML"))
      next
      }
    predDiffs <- predict(model, newdata = predDat) %>% as.vector()
    correctedRetentionTimes <- retentionTimes %>% `-`(predDiffs) %>% `*`(60)
    
    
    #Fix overcorrections
    overcorrected <- c(FALSE, diff(correctedRetentionTimes) <= 0)
    #print(overcorrected)
    #print(correctedRetentionTimes,max = 5000)
    while (any(overcorrected)) {
      #Index of first overcorrected time
      oc <- which.max(overcorrected)
      
      #Overcorrected time, and all times after it, are shifted such that the overcorrected time is delta greater than the time previous
      delta <- 0.0001
      
      correction <- (correctedRetentionTimes[(oc - 1)] + delta) - correctedRetentionTimes[(oc)]
      correctedRetentionTimes[oc:length(correctedRetentionTimes)] <- correctedRetentionTimes[oc:length(correctedRetentionTimes)] + correction
      
      overcorrected <- c(FALSE, diff(correctedRetentionTimes) <= 0)
    }
    
    negTime <- correctedRetentionTimes <= 0
    correctedRetentionTimes[negTime] <- 0
    
    write.csv(diff(correctedRetentionTimes),"diffs.csv")
    
    stopifnot(all(correctedRetentionTimes >= 0))
    stopifnot(all(diff(correctedRetentionTimes) >= 0))
    
    
    #Format corrected times and save them in new mzXML files
    correctedFormattedRetentionTimes <- correctedRetentionTimes %>% paste0("PT", ., "S")
    
    for (i in 1:length(correctedFormattedRetentionTimes)) {
      xml_attr(scanNodes[i], 'retentionTime') <- correctedFormattedRetentionTimes[i]
    }
    
    stopifnot(all(xml_attr(scanNodes, 'retentionTime') == correctedFormattedRetentionTimes))
    write_xml(x=doc, file=paste(outputpath, file_root, "_Aligned.mzXML",sep=""))
  }
}


CorrectTargetFiles <- function(targetFiles, models, outputpath, rtIndex = 2){
  
  i = rtIndex
  for(file in targetFiles){
    targets <- read.csv(file, stringsAsFactors = F)
    times <- as.data.frame(targets[,i])
    colnames(times) <- 'baseline_RT'
    
    file_root <- file %>% basename %>% sub(".csv", "",.)
    model <- models[[file_root]]
      
    predDiffs <- predict(model, newdata = times) %>% as.vector()
    correctedRetentionTimes <- times %>% `-`(predDiffs)
    targets[,i] <- correctedRetentionTimes
    write.csv(targets, paste(outputpath,file_root,"_Aligned.csv", sep=""), row.names = F)
  }
}


#Generate the retention time residuals for each file against the baseline file
GenResidualsFromIsos <- function(targets, baseline){

  #Read isos and scans files
  isos <- list.files(pattern = "*isos.csv")
  scans <- list.files(pattern = "*scans.csv")
  file_roots <- c(isos %>% sub("_isos.csv", "",.),scans %>% sub("_MS_scans.csv", "",.)) %>% unique
  
  #Define empty dataframe to populate with retention time information for each anchor
  rt_frame <- matrix(0, length(file_roots), dim(targets)[1]) %>% as.data.frame
  rownames(rt_frame) <- file_roots
  colnames(rt_frame) <- targets[,3]
  
  #For each file try to find each anchor target and add it to the dataframe
  for(c_file in file_roots){
    current_isos <- paste(c_file,"_isos.csv", sep="")
    current_scans <- paste(c_file,"_MS_scans.csv", sep="")
    
    current_isos <- read.csv(current_isos)
    current_scans <- read.csv(current_scans)
    
    for(i in 1:dim(targets)[1])
    {
      target <- targets[i,]
      mzTol <- 0.01 #0.01 m/z tolerance
      rtTol <- 0.3 #0.3 minutes retention time tolerance
      result <- FindTarget(target, current_isos, current_scans, mzTol, rtTol)
      if(!is.na(result)){
        rt_frame[c_file,target[[3]]] <- result$scan_time
      }
      else{
        rt_frame[c_file,target[[3]]] <- NA
      }
    }
  }
  
  #select the anchor points that did not have NA entries and computer the retention time residuals for all files from baseline
  complete_frame <- rt_frame[,complete.cases(t(rt_frame))]
  complete_frame <- t(complete_frame) %>% as.data.frame %>% select(starts_with(baseline), everything())
  rtResiduals <- t(complete_frame[,2:ncol(complete_frame)] - complete_frame[,1]) %>% as.data.frame
  
  x <- melt(as.matrix(rtResiduals))
  rts <- complete_frame[x$Var2,1]
  x$Var3 <- rts
  colnames(x) <- c("Dataset","Lipid","diff_RT","baseline_RT")
  
  return(x)
}

GenResidualsFromAnchors <- function(targets, baseline){

  mat <- targets[,-1]
  rownames(mat) <- targets[,1] %>% sub(".mzXML", "",.)
  #base <- matrix(mat[baseline,],nrow(mat), ncol(mat), byrow = T) %>% as.data.frame
  base <- mat[baseline,]
  rtResiduals <- sweep(mat, 2, as.numeric(base))
  
  x <- melt(as.matrix(rtResiduals))
  rts <- base[x$Var2] %>% as.numeric
  x$Var3 <- rts
  colnames(x) <- c("Dataset","Lipid","diff_RT","baseline_RT")
  return(x)
}


#Search for a target feature from the targets list in the current file
FindTarget <- function(target, current_isos, current_scans, mzTol, rtTol){
  current_file <- join(current_isos, current_scans, by="scan_num")
  feature_match <- filter(current_file, mz < target$row.m.z+mzTol & mz > target$row.m.z-mzTol & 
                scan_time < target$row.retention.time+rtTol & scan_time > target$row.retention.time-rtTol)
  feature_match <- feature_match[which.max(feature_match$abundance),]
  if(dim(feature_match)[1] == 0) {feature_match <- NA}
  return(feature_match)
}
