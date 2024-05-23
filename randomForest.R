# This script is for biomarker paper 2
library(AIMS)
library(org.Hs.eg.db)
library(randomForestSRC)
library(survival)
library(survminer)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggpubr)
library(GEOquery)
library(ComplexHeatmap)
library(dplyr)
library(circlize)

date <- Sys.Date()
project <- 'RandomForestSRC'

setwd("/mnt/Jason/HTG/data")

source('/mnt/Jason/HTG/scripts/mk.KM.plot.R')

dta <- read.csv("/mnt/Jason/HTG/data/all.norm.htg+clinical.csv", stringsAsFactors = TRUE)

# newSurvival <- read.csv('/mnt/Jason/HTG/data/Updated_PFS_04292024.csv')
# 
# # Update censored data
# newSurvival$pfs_time[newSurvival$true.progression == 0] <- newSurvival$pfs_time[newSurvival$true.progression == 0] + 7.594

rownames(dta) <- paste(dta$pcode,  dta$Timepoint, dta$Organ, dta$Block, sep = '.')

dta <- dta[dta$drug %in% c('AI', 'Fulvestrant'), ]

# P97.2.Soft.Tissue.A3 is missing in the raw data file. Must have a reason to be removed by Emily
#dta <- dta[!rownames(dta) %in% c('P97.2.Soft.Tissue.A3'), ]

# only use pre-treatment samples. Some pcode provided more than one biopsies for sequencing,
# we will keep them to increase our training dataset size
dtaPretx <- dta[grepl("2", dta$Timepoint), ]



# Set a color pallete for plotting use
colors <- c('#ffe119', '#3cb44b', '#e6194b', '#4363d8', '#f58231', 
            '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
            '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
            '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
            '#ffffff', '#000000')


# use published biomarkers for feature selection
biomarker <- list(
    PALOMA3 = './biomarkers/PALOMA3.txt',
    PALOMA23 = './biomarkers/PALOMA2-3_merged.txt',
    PEARL = './biomarkers/PEARL.txt',
    CDK2 = "./biomarkers/CDK2_signature.txt",
    OncotypeDX = "./biomarkers/OncotypeDX.txt",
    IHC4 = "./biomarkers/IHC4.txt",
    RBSig =  "./biomarkers/RBsig.txt",
    PAM50 = "./biomarkers/PAM50.txt",
    EndoPredict = "./biomarkers/EndoPredict.txt",
    Mammaprint = "./biomarkers/Mammaprint.txt",
    Mammostrat = "./biomarkers/Mammostrat.txt",
    Proliferation = "./biomarkers/ProliferationSignatureByCharlesPerou.txt",
    G2M_Checkpoint = "./biomarkers/G2M_Checkpoint.txt",
    ElasticNetSignature = "./biomarkers/ElasticNetSignature.txt",
    CDK46 = "./biomarkers/CDK4_6_GeneList.txt",
    BCI = "./biomarkers/BCI.txt"
)

biomarkerList <- list()

# read biomarkers
for (name in names(biomarker)){
    markerDf <- read.csv(biomarker[[name]], header = FALSE, stringsAsFactors = FALSE, sep = '\t')
    biomarkerList[[name]] <- markerDf$V1
}

setwd("/mnt/Jason/HTG/output/RandomForests")




vimp1 <- NULL

##=======================================================================================
## Round 1, initial biomarker evaluation
##=======================================================================================
# loop through each biomarker and perform feature selection


bestModelList.r1 <- list()


for (marker in names(biomarkerList)){
    genes.of.interest <- biomarkerList[[marker]]
    
    idx <- which(colnames(dtaPretx) %in% genes.of.interest)
    
    # scale the expression matrix
    tmp <- dtaPretx[, idx]
    tmp <- t(scale(t(tmp)))
    tmp <- as.data.frame(tmp)
    
    dtaSubset <- cbind(dtaPretx[, c(2,4)], tmp)
    
    set.seed(420)
    min.cindex <- 1
    bestModel.r1 <- NULL
    
    for (i in 1:10){
      # tune mtry and node size
      o1 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset)
      
      set.seed(420)
      obj.src <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset, ntree = 1000, importance = TRUE, mtry = o1[[2]][[2]], node.size = o1[[2]][[1]], na.action = "na.impute")
      
      sink(file = paste(date, project, marker, 'round1_run', i, 'traing_stats.txt', sep = '_'), append = FALSE, type = 'output')
      print(marker)
      print(obj.src)
      sink()
      
      # get OOB perforance error for best model selection
      cindex <- get.cindex(time = dtaSubset$pfs_time, censoring = dtaSubset$true.progression,
                           predicted = obj.src$predicted.oob)
      
      if (cindex < min.cindex){
        min.cindex <- cindex
        bestModel.r1 <- obj.src
      }
    }
    
    bestModelList.r1[[marker]] <- bestModel.r1
    
    # get variable importance
    print(predict(bestModel.r1, importance = TRUE)$importance)

    # plot OOB error rate and variable importance
    png(paste(date, project, marker, "VIM.png", sep = '_'),
        width = 8,
        height = 20,
        units = 'in',
        res = 200)
    plot(bestModel.r1)

    dev.off()


    png(paste(date, project, marker, "Survival.png", sep = '_'),
        width = 10,
        height = 10,
        units = 'in',
        res = 200)
    plot.survival(bestModel.r1, cens.model = "rfsrc", collapse = F)
    dev.off()
    
    
    # variable confidence interval
    set.seed(420)
    jk.obj <- subsample(bestModel.r1)
    
    # adjust plot height dynamicly
    height <- length(idx) * 0.22 + 1.3
    
    png(paste(date, project, marker, "Round1_VimpConfidenceInterval.png", sep = '_'), width = 5, height = height, units = 'in', res = 100)
    par(oma=c(1,1,1,1))
    par(mar=c(4,6,0,0))
    plot(jk.obj, cex = 0.6)
    dev.off()

    
    # VIMP
    vimpbytreatment <- as.data.frame(vimp(bestModel.r1, block.size = 10, joint = F)$importance)
    colnames(vimpbytreatment) <- "VIMP"
    vimpbytreatment <- vimpbytreatment[order(vimpbytreatment$VIMP, decreasing = TRUE), , drop = FALSE]
    
    # save to file
    write.table(vimpbytreatment,
                file = paste(date, project, marker, "Vimp.txt", sep = '_'),
                quote = FALSE,
                col.names = F,
                sep = "\t"
    )
    
    #
    # Test model performance with K-M plot
    #
    
    for (dg in c("Combined", "AI", "Fulvestrant")){
      survDf <- bestModel.r1$yvar
      survDf$OOB_Mortality <- bestModel.r1$predicted.oob
      survDf$Mortality <- bestModel.r1$predicted
      
      if (dg %in% c('AI', 'Fulvestrant')){
        survDf <- merge(survDf, dtaPretx[, 6, drop = F], by = 0)
        survDfSubset <- survDf %>% filter(drug == dg)
      } else {
        survDfSubset <- survDf
      }
      
      survDfSubset$Group <- ifelse(survDfSubset$OOB_Mortality > median(survDfSubset$OOB_Mortality), 'HighRisk', 'LowRisk')
      survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
      
      km <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                         event = 'true.progression', 
                         legendLabs = c('LowRisk', 'HighRisk'))
      
      png(paste(date, project, marker, dg, "R1_RSFModelPredictedOobSurvival_K-M_plot.png", sep = '_'),
          height = 6,
          width = 6,
          units = 'in',
          res = 100)
      print(km)
      dev.off()
      
      # IB data
      survDfSubset$Group <- ifelse(survDfSubset$Mortality > median(survDfSubset$Mortality), 'HighRisk', 'LowRisk')
      survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
      
      km2 <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                         event = 'true.progression', 
                         legendLabs = c('LowRisk', 'HighRisk'))
      
      png(paste(date, project, marker, dg, "R1_RSFModelPredictedIbSurvival_K-M_plot.png", sep = '_'),
          height = 6,
          width = 6,
          units = 'in',
          res = 100)
      print(km2)
      dev.off()
      
    }
    
    # save variables passing VIMP > 0.002 for next round
    tmp <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]
    vimp1 <- c(vimp1, rownames(tmp))
}

for (marker in names(biomarkerList)){
  bestModel.r1 <- bestModelList.r1[[marker]]
  
  survDf <- bestModel.r1$yvar
  survDf$Mortality <- bestModel.r1$predicted
  survDf$Group <- ifelse(survDf$Mortality > median(survDf$Mortality), 'HighRisk', 'LowRisk')
  survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))
  
  km2 <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time',
                      event = 'true.progression', 
                      legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, marker, "Combined_R1_RSFModelPredictedIbSurvival_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km2)
  dev.off()
  
  
}


##=======================================================================================
## Round 2, combine the important genes from round1 and repeat random survival modeling
## to further purify important genes
##=======================================================================================
important.genes <- unique(vimp1)

length(important.genes) 
# 140

# tune mtry and node size
idx2 <- which(colnames(dtaPretx) %in% important.genes)

tmp <- dtaPretx[, idx2]
tmp <- t(scale(t(tmp)))
tmp <- as.data.frame(tmp)
dtaSubset2 <- cbind(dtaPretx[, c(2,4)], tmp)

min.cindex <- 1
models.r2 <- NULL

# repeat 5 times
for (i in 1:10){

  o2 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset2)
  
  obj.src2 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset2, ntree = 1000, 
                    importance = TRUE, mtry = o2[[2]][[2]], node.size = o2[[2]][[1]], 
                    na.action = "na.impute")
  print(obj.src2)
  
  sink(file = paste(date, project, 'round2_run', i, 'traing_stats.txt', sep = '_'), append = F, type = 'output')
  print('Round2')
  print(obj.src2)
  sink()
  
  # get OOB perforance error for best model selection
  cindex <- get.cindex(time = dtaSubset2$pfs_time, censoring = dtaSubset2$true.progression,
                       predicted = obj.src2$predicted.oob)
  #rankings <- c(rankings, cindex)
  if (cindex < min.cindex){
    min.cindex <- cindex
    bestModel.r2 <- obj.src2
  }
}

png(paste(date, project, "Round2_Survival.png", sep = '_'),
    width = 10,
    height = 10,
    units = 'in',
    res = 200)
rv2 <- plot.survival(bestModel.r2, cens.model = "rfsrc", collapse = F)
dev.off()

# VIMP
vimpbytreatment <- as.data.frame(vimp(bestModel.r2, block.size = 10, joint = F)$importance)
colnames(vimpbytreatment) <- "VIMP"

# get variables passing VIMP > 0.002 for next round
vimp2 <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]

dim(vimp2) # 56 1

# variable confidence interval
set.seed(420)
jk.obj2 <- subsample(obj.src2)

height <- length(idx2) * 0.1 + 0.5
png(paste(date, project, "Round2_FilteredFeatureVimpConfidenceInterval.png", sep = '_'), width = 4.5, height = height, units = 'in', res = 100)
par(oma=c(1,1,1,1))
par(mar=c(4,6,0,0))
plot(jk.obj2, cex = 1.2)
dev.off()

important.genes2 <- rownames(vimp2)
length(important.genes2)

# # use the filtered genes to train a new model for later use
# dtaSubset2.1 <- dtaSubset2[, c('true.progression', 'pfs_time', important.genes2)]
# set.seed(420)
# o2.1 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset2.1)
# 
# set.seed(420)
# obj.src2.1 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset2.1, ntree = 1000, 
#                   importance = TRUE, mtry = o2.1[[2]][[2]], node.size = o2.1[[2]][[1]], 
#                   na.action = "na.impute")
# print(obj.src2.1)

##################################################################################################
# Make K-M plot to show the RSF performance by stratify patients with OOB predicted mortality
##################################################################################################

for (dg in c('Combined', 'AI', 'Fulvestrant')){
  survDf <- bestModel.r2$yvar
  survDf$OOB_Mortality <- bestModel.r2$predicted.oob
  survDf$Mortality <- bestModel.r2$predicted
  
  if (dg %in% c('AI', 'Fulvestrant')){
    survDf <- merge(survDf, dtaPretx[, 6, drop = F], by = 0)
    survDfSubset <- survDf %>% filter(drug == dg)
  } else {
    survDfSubset <- survDf
  }
  
  survDfSubset$Group <- ifelse(survDfSubset$OOB_Mortality > median(survDfSubset$OOB_Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "R2_RSFModelPredictedOobSurvival_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km)
  dev.off()
  
  #
  # All sample (IB + OOB)
  #
  survDfSubset$Group <- ifelse(survDfSubset$Mortality > median(survDfSubset$Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km2 <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "R2_RSFModelPredictedIbSurvival_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km2)
  dev.off()
}


# save the round2 model for testing in other dataset with no clinical varibles
#save(obj.src2.1, file = paste(date, project, 'Round2RSFModel.Rdata', sep = '_'))


##=======================================================================================
## Round 3, add clinical variables to the important genes selected from round2 and repeat
## random survival modeling to further purify important genes and clinical variables
##=======================================================================================
dtaSubset.cli <- cbind(dtaPretx[, c(2, 4, 6:12, 15)], tmp)

# Organ may not be selected in the important feature list, so handle it
if ('Organ' %in% colnames(dtaSubset.cli)){
    dtaSubset.cli$Organ <- as.character(dtaSubset.cli$Organ)
    dtaSubset.cli$Organ[!dtaSubset.cli$Organ %in% c('Bone', 'Lung', 'Breast', 'Liver', 'Skin')] <- 'Other'
    dtaSubset.cli$Organ <- as.factor(dtaSubset.cli$Organ)
}

min.cindex <- 1
models.r3 <- NULL

# repeat 5 times
for (i in 1:10){
  o3 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli)

  obj.src3 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli, ntree = 1000, 
                    importance = TRUE, mtry = o3[[2]][[2]], node.size = o3[[2]][[1]], 
                    na.action = "na.impute")
  print(obj.src3)
  sink(file = paste(date, project, 'round3_run', i, 'traing_stats.txt', sep = '_'), append = F, type = 'output')
  print('Round3')
  print(obj.src3)
  sink()
  
  # get OOB perforance error for best model selection
  cindex <- get.cindex(time = dtaSubset.cli$pfs_time, censoring = dtaSubset.cli$true.progression,
                       predicted = obj.src3$predicted.oob)
  #rankings <- c(rankings, cindex)
  if (cindex < min.cindex){
    min.cindex <- cindex
    bestModel.r3 <- obj.src3
  }
}


png(paste(date, project, "Round3_Survival.png", sep = '_'),
    width = 10,
    height = 10,
    units = 'in',
    res = 200)
rv3 <- plot.survival(bestModel.r3, cens.model = "rfsrc", collapse = F)
dev.off()

# VIMP
vimpbytreatment <- as.data.frame(vimp(bestModel.r3, block.size = 10, joint = F)$importance)
colnames(vimpbytreatment) <- "VIMP"

# get variables passing VIMP > 0.002 for next round use
vimp3 <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]
dim(vimp3)

important.features3 <- rownames(vimp3)

# 
# variable confidence interval
set.seed(420)
jk.obj3 <- subsample(bestModel.r3)
height <- length(idx2) * 0.07 + 0.2
png(paste(date, project, "R3_VimpConfidenceInterval.png", sep = '_'), width = 4.5, height = height, units = 'in', res = 100)
par(oma=c(1,1,1,1))
par(mar=c(4,6,0,0))
plot(jk.obj3, cex = 1.2)
dev.off()


##################################################################################################
# Make K-M plot to show the RSF performance by stratify patients with OOB predicted mortality
##################################################################################################

for (dg in c('Combined', 'AI', 'Fulvestrant')){
  survDf <- bestModel.r3$yvar
  survDf$OOB_Mortality <- bestModel.r3$predicted.oob
  survDf$Mortality <- bestModel.r3$predicted
  
  if (dg %in% c('AI', 'Fulvestrant')){
    survDf <- merge(survDf, bestModel.r3$xvar[, 1, drop = FALSE], by = 0)
    survDfSubset <- survDf %>% filter(drug == dg)
  } else {
    survDfSubset <- survDf
  }
  
  survDfSubset$Group <- ifelse(survDfSubset$OOB_Mortality > median(survDfSubset$OOB_Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "R3_RSFModelPredictedOobSurvival_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km)
  dev.off()
  
  # IB data 
  survDfSubset$Group <- ifelse(survDfSubset$Mortality > median(survDfSubset$Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "R3_RSFModelPredictedIbSurvival_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km)
  dev.off()
}


##=======================================================================================
## Round 4, repeat random survival modeling to verify features selected and further
## eliminate less important features
##=======================================================================================
dtaSubset.cli2 <- dtaSubset.cli[, c('true.progression', 'pfs_time', important.features3)]

# due to the algorithmic nature of random forest, each time a tree is built, a slightly
# different vimp list is returned. If multiple features have similar resolving power, they
# may have different vimp values. However, the strongest feature should always the same.
# We will run the below code black five times and pick the best model (with the lowest OOB
# error)

bestModel.r4 <- NULL
min.cindex <- 1

for (i in 1:10){ # repeat 5 times
   # set.seed(420)
    o4 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli2)
    
  #  set.seed(420)
    obj.src4 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli2, ntree = 1000, 
                      importance = TRUE, mtry = o4[[2]][[2]], node.size = o4[[2]][[1]], 
                      na.action = "na.impute")
    print(obj.src4)
    sink(file = paste(date, project, 'round4_run', i, 'traing_stats.txt', sep = '_'), append = F, type = 'output')
    print("Round4")
    print(obj.src4)
    sink()
    
    # get OOB perforance error for best model selection
    cindex <- get.cindex(time = dtaSubset.cli2$pfs_time, censoring = dtaSubset.cli2$true.progression,
               predicted = obj.src4$predicted.oob)
    #rankings <- c(rankings, cindex)
    if (cindex < min.cindex){
      min.cindex <- cindex
      bestModel.r4 <- obj.src4
    }
}

# 
# variable confidence interval
#set.seed(420)
jk.obj4 <- subsample(bestModel.r4)
height <- length(idx2) * 0.07 + 0.2
png(paste(date, project, "R4_VimpConfidenceInterval.png", sep = '_'), width = 4.5, height = height, units = 'in', res = 100)
par(oma=c(1,1,1,1))
par(mar=c(4,6,0,0))
plot(jk.obj3, cex = 1.2)
dev.off()

##############################################################################################
# Plot predicted survival function dichotomize wit OOB predicted mortality
##############################################################################################
pred.surv =  as.data.frame(t(bestModel.r4$survival.oob))
colnames(pred.surv) <- rownames(bestModel.r4$xvar)
pred.surv$Time <- bestModel.r4$time.interest

dd = melt(pred.surv, id=c("Time"))
colnames(dd) <- c('Time', 'PID', 'Survival')

mortality <- data.frame(Mortality = bestModel.r4$predicted.oob)
mortality$Risk <- ifelse(mortality$Mortality > median(mortality$Mortality), 'High', 'Low')
rownames(mortality) <- rownames(bestModel.r4$yvar)

dd$Risk <- 'Low'
dd$Risk[dd$PID %in% rownames(mortality[mortality$Risk == 'High', ])] <- 'High'

p <- ggplot(dd) + geom_line(aes(x=Time, y=Survival, color = Risk, group = PID), linetype = 2) + 
  labs(title="RSF Predicted Cohort Survival") + 
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12), 
        legend.position = c(0.8, 0.85)) +
  xlab('Time (Months)') + ylab('OOB Survival')

png(paste(date, project, "Round4PredictedSurvival.png", sep = '_'), 
    width = 5, height = 5, units = 'in', res = 200)
print(p)
dev.off() 

png(paste(date, project, "Round4_Survival.png", sep = '_'),
    width = 10,
    height = 10,
    units = 'in',
    res = 200)
rv4 <- plot.survival(bestModel.r4, cens.model = "rfsrc", collapse = F) #, subset = subset$Risk == 'High')
dev.off()



# VIMP
vimpbytreatment <- as.data.frame(vimp(bestModel.r4, block.size = 10, joint = F)$importance)
colnames(vimpbytreatment) <- "VIMP"

# save variables passing VIMP > 0.002. This should have no effect
vimp4 <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]
dim(vimp4)



# save to file
write.table(x = vimp4,
            file = paste(date, project, 'Round4_FilteredFeatures.csv', sep = '_'),
            quote = FALSE,
            sep = ",",
            row.names = TRUE)

##################################################################################################
# Make K-M plot to show the RSF performance by stratify patients with OOB predicted mortality
##################################################################################################

for (dg in c('Combined', 'AI', 'Fulvestrant')){
  survDf <- bestModel.r4$yvar
  survDf$OOB_Mortality <- bestModel.r4$predicted.oob
  survDf$Mortality <- bestModel.r4$predicted
  
  # join with OS data
  survDf <- merge(survDf, dtaPretx[, c(3,5)], by = 0)
  rownames(survDf) <- survDf$Row.names
  survDf$Row.names <- NULL
  
  if (dg %in% c('AI', 'Fulvestrant')){
    survDf <- merge(survDf, bestModel.r4$xvar[, 1, drop = FALSE], by = 0)
    survDfSubset <- survDf %>% filter(drug == dg)
  } else {
    survDfSubset <- survDf
  }
  
  survDfSubset$Group <- ifelse(survDfSubset$OOB_Mortality > median(survDfSubset$OOB_Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "Round4_ModelPredictedOobSurvival_PFS_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km)
  dev.off()
  
  # overall survival
  km2 <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'precdk_os_time',
                     event = 'os', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "Round4_ModelPredictedSurvival_OS_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km2)
  dev.off()
  
  # IB data
  survDfSubset$Group <- ifelse(survDfSubset$Mortality > median(survDfSubset$Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km3 <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "Round4_ModelPredictedIbSurvival_PFS_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km3)
  dev.off()
  
}


# Best model reduced, i.e. with only the VIMP slected features
#set.seed(420)
min.cindex <- 1
bestModel.reduced <- NULL
final.features <- rownames(vimp4)
dtaSubset.final <- dtaSubset.cli2[, c('true.progression', 'pfs_time', final.features)]

for (i in 1:10){
  ofn <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset.final)

  #  set.seed(420)
  obj.src <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset.final, ntree = 1000,
                   importance = TRUE, mtry = ofn[[2]][[2]], node.size = ofn[[2]][[1]],
                   na.action = "na.impute")
  print(obj.src)
  sink(file = paste(date, project, 'bestReduced_run', i, 'traing_stats.txt', sep = '_'), append = TRUE, type = 'output')
  print("BestReduced")
  print(obj.src)
  sink()

  #set.seed(420)
  jk.obj <- subsample(obj.src)

  height <- nrow(vimp4) * 0.24 + 1.3
  png(paste(date, project, "BestReducedModel_Run", i, "ConfidenceInterval.png", sep = '_'), width = 5, height = height, units = 'in', res = 100)
  par(oma=c(1,1,1,1))
  par(mar=c(4,6,2,2) + 0.1)
  plot(jk.obj, cex = 1.2)
  dev.off()

  # get OOB perforance error for best model selection
  cindex <- get.cindex(time = dtaSubset.final$pfs_time, censoring = dtaSubset.final$true.progression,
                       predicted = obj.src$predicted.oob)
  #rankings <- c(rankings, cindex)
  if (cindex < min.cindex){
    min.cindex <- cindex
    bestModel.reduced <- obj.src
  }
}

# plot OOB predicted survival function stratify patients by median dichotimze mortality score
pred.surv =  as.data.frame(t(bestModel.reduced$survival.oob))
colnames(pred.surv) <- rownames(bestModel.reduced$xvar)
pred.surv$Time <- bestModel.reduced$time.interest

dd = melt(pred.surv, id=c("Time"))
colnames(dd) <- c('Time', 'PID', 'Survival')

mortality <- data.frame(Mortality = bestModel.reduced$predicted.oob)
mortality$Risk <- ifelse(mortality$Mortality > median(mortality$Mortality), 'High', 'Low')
rownames(mortality) <- rownames(bestModel.reduced$yvar)

dd$Risk <- 'Low'
dd$Risk[dd$PID %in% rownames(mortality[mortality$Risk == 'High', ])] <- 'High'

p <- ggplot(dd) + geom_line(aes(x=Time, y=Survival, color = Risk, group = PID), linetype = 2) + 
  labs(title="RSF Predicted Cohort Survival") + 
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12), 
        legend.position = c(0.8, 0.85)) +
  xlab('Time (Months)') + ylab('OOB Survival')

png(paste(date, project, "FinalRoundPredictedSurvival.png", sep = '_'), 
    width = 5, height = 5, units = 'in', res = 200)
print(p)
dev.off() 

png(paste(date, project, "FinalRound_Survival.png", sep = '_'),
    width = 10,
    height = 10,
    units = 'in',
    res = 200)
rv5 <- plot.survival(bestModel.reduced, cens.model = "rfsrc", collapse = F) 
dev.off()


##################################################################################################
# Make K-M plot to show the RSF performance by stratify patients with OOB predicted mortality for
# bestModel reduced
##################################################################################################

for (dg in c('Combined', 'AI', 'Fulvestrant')){
  survDf <- bestModel.reduced$yvar
  survDf$OOB_Mortality <- bestModel.reduced$predicted.oob
  survDf$Mortality <- bestModel.reduced$predicted
  
  # join with OS data
  survDf <- merge(survDf, dtaPretx[, c(3,5)], by = 0)
  rownames(survDf) <- survDf$Row.names
  survDf$Row.names <- NULL
  
  if (dg %in% c('AI', 'Fulvestrant')){
    survDf <- merge(survDf, bestModel.reduced$xvar[, 1, drop = FALSE], by = 0)
    survDfSubset <- survDf %>% filter(drug == dg)
  } else {
    survDfSubset <- survDf
  }
  
  survDfSubset$Group <- ifelse(survDfSubset$OOB_Mortality > median(survDfSubset$OOB_Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "ReducedBestModelPredictedOobSurvival_PFS_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km)
  dev.off()
  
  # overall survival
  km2 <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'precdk_os_time',
                      event = 'os', 
                      legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "ReducedBestFModelPredictedSurvival_OS_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km2)
  dev.off()
  
  # IB data
  survDfSubset$Group <- ifelse(survDfSubset$Mortality > median(survDfSubset$Mortality), 'HighRisk', 'LowRisk')
  survDfSubset$Group <- factor(survDfSubset$Group, levels = c('LowRisk', 'HighRisk'))
  
  km3 <- make.KM.plot(survDf = survDfSubset, grouping.var = 'Group', surv.type = 'pfs_time',
                     event = 'true.progression', 
                     legendLabs = c('LowRisk', 'HighRisk'))
  
  png(paste(date, project, dg, "ReducedBestModelPredictedIbSurvival_PFS_K-M_plot.png", sep = '_'),
      height = 6,
      width = 6,
      units = 'in',
      res = 100)
  print(km3)
  dev.off()
  
}


##=======================================================================================
## Scatter plot to show correlation of mortality with actual PFS time
##=======================================================================================
data.for.scatterplot <- data.frame(Mortality = bestModel.reduced$predicted.oob, PFS_time = dtaSubset.final$pfs_time, 
                                   Status = ifelse(dtaSubset.final$true.progression == 0, 'Censored', 'Event'),
                                   Sample = rownames(bestModel.reduced$yvar))


# Some patients progressed since last follow up, so plot them on top of the above plot and see if
# the new data point are moving closer to the trend line
data.for.scatterplot$pcode <- sub('^(P\\d+).*', '\\1', data.for.scatterplot$Sample)


p <- ggscatter(data.for.scatterplot, x = "PFS_time", y = "Mortality",  shape = 'Status',
               color = "Status", bg = 'black', fill = 'Status', size = ifelse(data.for.scatterplot$Status == 'NewlyProgressed', 3, 1.5), # Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "pearson", label.x = 50, label.sep = "\n")
)
p <- p + scale_shape_manual(values = c(1, 2, 8)) # + geom_point(data = newlyProgressed, aes(x=PFS_time,y=Mortality),colour="purple", shape = 9, size = 3)
p <- p + ylab('OOB Mortality')
png(paste(date, project, "OOB_PredictedMortality_vs_PFS_on_RPCI_preTX.png", sep = '_'),
    width = 5,
    height = 5,
    units = 'in',
    res = 150
)
print(p)
dev.off()


# boxplot to show difference of predicted mortality for very short PFS (intrinsic resistance) and very long PFS patients
data.for.boxplot <- data.for.scatterplot %>% filter(PFS_time > 25 | Status == 'Event'& PFS_time < 6)
data.for.boxplot$Response <- ifelse(data.for.boxplot$PFS_time < 6, 'Resistant', 'Responsive')
data.for.boxplot$Response <- factor(data.for.boxplot$Response, levels = c('Resistant', 'Responsive'))

p <- ggboxplot(data.for.boxplot, x = "Response", y = "Mortality",
               color = "Response", palette = "jco",
               add = "jitter") + theme(legend.position="none", axis.title.x = element_blank()) +
  ylab('OOB Predicted Mortality') + stat_compare_means(method = "t.test")

png(paste(date, project, "OOB_PredictedMortalityResistantVsReponsivePfsRPCI_Boxplot.png", sep = '_'),
    width = 3,
    height = 4,
    units = 'in',
    res = 200
)
print(p)
dev.off()


# use (IB) predicted results on the whole ensemble

data.for.scatterplot2 <- data.frame(Mortality = bestModel.reduced$predicted, PFS_time = dtaSubset.final$pfs_time, 
                                   Status = ifelse(dtaSubset.final$true.progression == 0, 'Censored', 'Event'),
                                   #Status2 = bestModel.reduced$event.info$cens,
                                   Sample = rownames(bestModel.reduced$yvar))


p2 <- ggscatter(data.for.scatterplot2, x = "PFS_time", y = "Mortality",  shape = 'Status',
                color = "Status", bg = 'black', fill = 'Status', size =  1.5, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 50, label.sep = "\n")
)
p2 <- p2 + scale_shape_manual(values = c(1, 2, 8)) # + geom_point(data = newlyProgressed, aes(x=PFS_time,y=Mortality),colour="purple", shape = 9, size = 3)
p2 <- p2 + ylab('Mortality')
png(paste(date, project, "PredictedMortality_vs_PFS_on_RPCI_preTX.png", sep = '_'),
    width = 5,
    height = 5,
    units = 'in',
    res = 150
)
print(p2)
dev.off()

data.for.boxplot2 <- data.for.scatterplot2 %>% filter(PFS_time > 25 | Status2 == 'Event' & PFS_time < 6)
data.for.boxplot2$Response <- ifelse(data.for.boxplot2$PFS_time < 6, 'Resistant', 'Responsive')
data.for.boxplot2$Response <- factor(data.for.boxplot2$Response, levels = c('Resistant', 'Responsive'))

p2 <- ggboxplot(data.for.boxplot2, x = "Response", y = "Mortality",
               color = "Response", palette = "jco",
               add = "jitter") + theme(legend.position="none", axis.title.x = element_blank()) +
  ylab('Predicted Mortality') + stat_compare_means(method = "t.test", label.x = 1.5)

png(paste(date, project, "PredictedMortalityResistantVsResponsivePfsRPCI_Boxplot.png", sep = '_'),
    width = 3,
    height = 4,
    units = 'in',
    res = 200
)
print(p2)
dev.off()


##=======================================================================================
## Validation on NeoPalAna dataset which are treated with AI + Palbo. This dataset has no
## survival data available, however, it has 5 palbo resistant samples, so test if they 
## can all be predicted in the high risk patient group. One caveates is that this dataset
## is missing the clinical variables, so we will need to train a new module without them.
##=======================================================================================

#
# Get NeoPalAna dataset from GEO ussing GEOquery
#

# Avoid the following error by setting VROOM_CONNECTION_SIZE to a large enough value
# Error: The size of the connection buffer (262144) was not large enough                                   
# to fit a complete line:
#  * Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 6)
gse <- getGEO("GSE93204")
gse <- gse[[1]]
sampleInfo <- pData(gse)

# reduce to include only needed columns
sampleInfo <- sampleInfo[, c("characteristics_ch2.3", 'treatment:ch2')]
names(sampleInfo) <- c("PID", "Treatment")
sampleInfo$PID <- sub("patient id: ", '', sampleInfo$PID)

# replace probe id with gene names
features <- fData(gse)
features <- features[, "GENE_SYMBOL", drop = FALSE]
features$GENE_SYMBOL <- sub('CTSL2', 'CTSV', features$GENE_SYMBOL)

full_output <- cbind(features,exprs(gse))

samplesC0D1 <- sampleInfo[sampleInfo$Treatment == "Baseline", ]

c0d1 <- full_output[full_output$GENE_SYMBOL %in% colnames(dtaSubset.cli2)[-c(1:2)], c("GENE_SYMBOL", rownames(samplesC0D1))]

# remove rows with duplicate gene names
c0d1 <- c0d1[!duplicated(c0d1$GENE_SYMBOL), ]
rownames(c0d1) <- c0d1$GENE_SYMBOL
c0d1$GENE_SYMBOL <- NULL

# map patient IDs
tmp <- t(c0d1)
tmp <- merge(tmp, sampleInfo, by = 0)
rownames(tmp) <- tmp$PID
tmp <- subset(tmp, select=-c(PID, Treatment, Row.names))
c0d1 <- as.data.frame(t(scale(t(tmp))))

resistantPids <- c("102", "105", "106", "108", "113", "205")

# train a new model for this dataset by excluding the clinical variables (except for drug)
clinical.variables <- colnames(dtaPretx)[6:16]

dtaSubset.neopalana <- dtaSubset.cli2[, !colnames(dtaSubset.cli2) %in% clinical.variables]
dtaSubset.neopalana <- cbind(dtaSubset.neopalana[, 1:3], dtaSubset.neopalana[, colnames(dtaSubset.neopalana) %in% colnames(c0d1)])

o5 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset.neopalana)

obj.src5 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset.neopalana, ntree = 1000, 
                  importance = TRUE, mtry = o5[[2]][[2]], node.size = o5[[2]][[1]], 
                  na.action = "na.impute")

print(obj.src5)
# (OOB) Requested performance error: 0.38341323

#save(obj.src5, file = paste(date, project, 'RSFModelNoClinical.Rdata', sep = '_'))

# This dataset was treated with AI like drug, so add the drug column
c0d1$drug <- 'AI'
neopalana.predict <- predict(obj.src5, c0d1, na.action = 'na.impute')

pred.surv = t(neopalana.predict$survival)
colnames(pred.surv) <- rownames(neopalana.predict$xvar)
pred.surv <- as.data.frame(pred.surv)
pred.surv$Time <- neopalana.predict$time.interest

dd = melt(pred.surv, id=c("Time"))
colnames(dd) <- c('Time', 'PID', 'Survival')

dd$PID <- as.character(dd$PID)

dd$Response <- 'Responsive'
dd$Response[dd$PID %in% resistantPids] <- 'Resistant'

p <- ggplot(dd) + geom_line(aes(x=Time, y=Survival, color = Response, group = PID), linetype = 2) + 
    labs(title="RSF Predicted Survival on NeoPalAna Study") + 
    theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12), 
          legend.position = c(0.8, 0.8)) +
    xlab('Time (Months)') + ylab('PFS Probability')

png(paste(date, project, "NeoPalAnaPredictedSurvival.png", sep = '_'), 
    width = 5, height = 5, units = 'in', res = 200)
print(p)
dev.off() 


#
# Kolmogorov-Smirnov test to show CDK4/6i resistant tumors are skewed toward the short PFS end
#

# Use PFS probability values at 40 months (arbitrary, at this time point, we see good sevaration of the lines)
data.for.ks <- dd %>% filter(Time > 30, Time < 30.2)

# order by PFS time in increasing order
data.for.ks <- data.for.ks %>% arrange(Survival)
data.for.ks$rank <- 1:nrow(data.for.ks)

x = which(data.for.ks$Response == 'Resistant')
ks.test(x = x, y = data.for.ks$rank, alternative = 'greater')

# Exact two-sample Kolmogorov-Smirnov test
# 
# data:  x and data.for.ks$rank
# D^+ = 0.71875, p-value = 0.004857
# alternative hypothesis: the CDF of x lies above that of y


# Heatmap to show distribution
data.for.ks$Response2 <- ifelse(data.for.ks$Response == 'Resistant', 1, 0)
rownames(data.for.ks) <- data.for.ks$PID


response.color <- c('red', 'blue')
names(response.color) <- c('Resistant', 'Responsive')


dummy.data <- data.frame(dummy = rep(0, nrow(data.for.ks)))
rownames(dummy.data) <- data.for.ks$PID

ha <- HeatmapAnnotation(Response = data.for.ks$Response,
                        Survival = data.for.ks$Survival,
                        col = list(Response = response.color, Survival = colorRamp2(c(0, 0.25), c("white", "purple"))),
                        annotation_name_gp= gpar(fontsize = 12, fontface = 'bold'))
ht <- Heatmap(t(dummy.data),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              show_column_names = TRUE,
              show_row_names = FALSE,
              show_heatmap_legend = FALSE,
              bottom_annotation = ha,
              col = 'white',
              )

png(paste(date, project, "Neopalana_KS.png", sep = '_'),
    width = 7,
    height = 3,
    units = 'in',
    res = 100)
print(ht)
dev.off()


##===========================================================================================================
## Testing for predictive and prognostic power on Metabric dataset. If it has discriminating power, that 
## could indicate shared genes that contribute to both predictive and prognostic capabilities.
##===========================================================================================================

metabric.exp <- read.csv("/mnt/Jason/ExternalData/brca_metabric/data_mrna_illumina_microarray.txt", sep = '\t', row.names = 1)
metabric.exp$Entrez_Gene_Id <- NULL

# for AIMS subtype prediction
metabric.exp2 <- as.data.frame(t(metabric.exp))

# convert to z-scores
metabric.exp <- as.data.frame(t(scale(t(metabric.exp))))

# transpose to join with clinical data
metabric.exp <- as.data.frame(t(metabric.exp))

metabric.clinical.sample <- read.csv("/mnt/Jason/ExternalData/brca_metabric/data_clinical_sample.txt", sep = '\t', row.names = 1, skip = 4)
rownames(metabric.clinical.sample) <- sub('-', '.', rownames(metabric.clinical.sample))

hr.plus.her2.minus.sample <- metabric.clinical.sample %>% dplyr::filter(ER_STATUS == 'Positive', HER2_STATUS == 'Negative')

# survival data
metabric.survival <- read.csv('/mnt/Jason/ExternalData/brca_metabric/data_clinical_patient.txt', sep = '\t', row.names = 1, skip = 4)
rownames(metabric.survival) <- sub('-', '.', rownames(metabric.survival))

metabricSubset <- merge(metabric.survival, hr.plus.her2.minus.sample, by = 0)
metabricSubset <- merge(metabricSubset, metabric.exp, by.x = 'Row.names', by.y = 0)
rownames(metabricSubset) <- metabricSubset$Row.names

# for AIMS subtype prediction
metabric.for.aims <- metabric.exp2 %>% filter(rownames(metabric.exp2) %in% rownames(hr.plus.her2.minus.sample))
metabric.for.aims <- as.data.frame(t(metabric.for.aims))

metabricSubsetStandardized <- metabricSubset[, colnames(metabricSubset) %in% final.features[!final.features %in% colnames(dtaPretx)[5:16]]]

metabricOS <- metabricSubset[,  c("OS_STATUS", "OS_MONTHS")]
metabricRFS <- metabricSubset[,  c('RFS_STATUS', 'RFS_MONTHS')]

# fix first two column names. They need to match the variables used in the training ds
colnames(metabricOS) <- c('true.progression', 'pfs_time')
colnames(metabricRFS) <- c('true.progression', 'pfs_time')
metabricOS$true.progression <- sub('(\\d):.*', '\\1', metabricOS$true.progression)
metabricRFS$true.progression <- sub('(\\d):.*', '\\1', metabricRFS$true.progression)
metabricOS$true.progression <- as.numeric(metabricOS$true.progression)
metabricRFS$true.progression <- as.numeric(metabricRFS$true.progression)

metabricSubsetOS <- merge(metabricOS, metabricSubsetStandardized, by = 0)
metabricSubsetRFS <- merge(metabricRFS, metabricSubsetStandardized, by = 0)

rownames(metabricSubsetOS) <- metabricSubsetOS$Row.names
rownames(metabricSubsetRFS) <-metabricSubsetRFS$Row.names

metabricSubsetOS$Row.names <- NULL
metabricSubsetRFS$Row.names <- NULL

predList <- list()
predList[['RFS']] <- metabricSubsetRFS
predList[['OS']] <- metabricSubsetOS

# Since EBF4 and GTSE1 are not in metabric dataset, we will need to retrain a RFS
# model without these two genes

#
# Train a new model to exclude missing features in metabric ds
#

dtaSubset.for.metabric <- dtaSubset.cli2[, colnames(dtaSubset.cli2) %in% c('pfs_time',  'true.progression', colnames(metabricSubsetOS))]

o6 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset.for.metabric)

obj.src.metabric <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset.for.metabric, ntree = 3000, 
                          importance = TRUE, mtry = o6[[2]][[2]], node.size = o6[[2]][[1]],
                          na.action = "na.impute")
print(obj.src.metabric)

for (survType in names(predList)){
    survDf <- predList[[survType]]
 
    survDf <- survDf[, match(colnames(dtaSubset.for.metabric), colnames(survDf))]
    
    # predict mortality on new dataset
    obj.predict <- predict(obj.src.metabric, survDf)
    
    print(obj.predict)
    yhat <- obj.predict$predicted
    
    q.yhat <- quantile(yhat) #, probs = c(0, 1, 0.25))
    q.yhat[1] <- q.yhat[1] - .0001
    Group <- cut(yhat, q.yhat, labels = FALSE)
    
    survDf <- cbind(obj.predict$yvar, Group)
    survDf$Group[survDf$Group == 1] <- 'LowRisk'
    survDf$Group[survDf$Group == 2] <- 'mLowRisk'
    survDf$Group[survDf$Group == 3] <- 'mHighRisk'
    survDf$Group[survDf$Group == 4] <- 'HighRisk'
    survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'mLowRisk', 'mHighRisk', 'HighRisk'))
    
    colnames(survDf)[1] <- c(survType)
    p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = survType, event = 'true.progression')
    
    png(paste(date, project, survType, "Quartile_K-M_plot_metabric.png", sep = '_'),
        width = 8, 
        height = 8, 
        units = "in", 
        pointsize = 10, 
        res = 200
    )
    print(p)
    dev.off()
    
    # by median dichotomize
    survDf$Group2[survDf$Group == 'LowRisk'] <- 'LowRisk'
    survDf$Group2[survDf$Group == 'mLowRisk'] <- 'LowRisk'
    survDf$Group2[survDf$Group == 'mHighRisk'] <- 'HighRisk'
    survDf$Group2[survDf$Group == 'HighRisk'] <- 'HighRisk'
    survDf$Group2 <- factor(survDf$Group2, levels = c('LowRisk', 'HighRisk'))
    
    colnames(survDf)[1] <- c(survType)
    p <- make.KM.plot(survDf = survDf, grouping.var = 'Group2', surv.type = survType, event = 'true.progression')
    
    png(paste(date, project, survType, "Median_K-M_plot_metabric.png", sep = '_'),
        width = 6, 
        height = 6, 
        units = "in", 
        pointsize = 10, 
        res = 200
    )
    print(p)
    dev.off()
    
    # use subtype enrichment analysis guided cut point
    q.yhat <- quantile(yhat, probs = c(0:3)/3)
    q.yhat[1] <- q.yhat[1] - .0001
    Group <- cut(yhat, q.yhat, labels = FALSE)
    
    survDf <- cbind(obj.predict$yvar, Group)
    survDf$Group3[survDf$Group == 1] <- 'LowRisk'
    survDf$Group3[survDf$Group == 2] <- 'MediumRisk'
    survDf$Group3[survDf$Group == 3] <- 'HighRisk'
    survDf$Group3 <- factor(survDf$Group3, levels = c('LowRisk', 'MediumRisk', 'HighRisk'))
    
    colnames(survDf)[1] <- c(survType)
    p <- make.KM.plot(survDf = survDf, grouping.var = 'Group3', surv.type = survType, event = 'true.progression')
    
    png(paste(date, project, survType, "Tertile_K-M_plot_metabric.png", sep = '_'),
        width = 6, 
        height = 6, 
        units = "in", 
        pointsize = 10, 
        res = 200
    )
    print(p)
    dev.off()
}

##============================================================================================
## Test if this rsf model is able to separate LumA from other subtypes in POLOMA2 ds, etc.
##============================================================================================
# We need to train a new model with features available in Paloma dataset 
# read in the raw read counts from POLOMA2/3
expRawP2 <- read.csv("/mnt/Jason/HTG/externaldata/GSE133394_a5481008_raw_data.txt", sep = "\t", row.names = 1)
expRawP3 <- read.csv("/mnt/Jason/HTG/externaldata/GSE128500_a5481023_raw_data.csv", row.names = 1)
expRawPl <- read.csv("/mnt/Jason/HTG/externaldata/GSE223700_Raw_data.csv", row.names = 1)

na.sums <- apply(expRawPl, 2, function(x) sum(is.na(x)))
to.keep <- na.sums < 2500
expRawPl <- expRawPl[, to.keep]

expNormP2 <- read.csv("/mnt/Jason/HTG/externaldata/2023-07-31_PALOMA-2_edgeR.normalized.csv", row.names = 1)
expNormP3 <- read.csv("/mnt/Jason/HTG/externaldata/2023-07-31_PALOMA-3_edgeR.normalized.csv", row.names = 1)
expNormPl <- read.csv("/mnt/Jason/HTG/externaldata/2023-07-31_PEARL_edgeR.normalized.csv", row.names = 1)

expRawList <- list()
expRawList[['PALOMA-2']] <- expRawP2
expRawList[['PALOMA-3']] <- expRawP3
expRawList[['PEARL']] <- expRawPl 

# this is not raw counts, included here for comparition with subtype counts from the origianl data. 
# Will be replaced with the original data
expRawList[['Metabric']] <- metabric.for.aims


expNormList <- list()
expNormList[['PALOMA-2']] <- expNormP2
expNormList[['PALOMA-3']] <- expNormP3
expNormList[['PEARL']] <- expNormPl
expNormList[['Metabric']] <- as.data.frame(t(metabricSubset[, 39:ncol(metabricSubset)]))
#expNormList[['HTG']] <- as.data.frame(t(dtaPretx[, 18:ncol(dtaPretx)]))

#
# Map Entrez gene id to express matrix. This is needed by AIMS
#
id2symbol <- as.data.frame(org.Hs.egALIAS2EG)
id2symbol <- id2symbol[!duplicated(id2symbol$alias_symbol), ]
rownames(id2symbol) <- id2symbol$alias_symbol

# fix Entrez gene id problem
myDict <- list(CES2 = 8824, CCR9 = 10803, CCR10 = 2826, ADRA1A = 148, AHRR = 57491, RAD1 = 5810, MADD = 8567, AR = 367, SLK = 9748,
               MIF = 4282, NLK = 51701, MAGEE1 = 56792, RAD54L = 8438, TEP1 = 7011, SRA1 = 10011, TPO = 7173, TTF1 = 7270)

for (n in names(myDict)){
    id2symbol$gene_id <- ifelse(id2symbol$alias_symbol == n, myDict[[n]], id2symbol$gene_id)
}

##################################################################################################
# 1. Identify AIMS subtypes for all samples
##################################################################################################
subtypeList <- list()

for (ds in names(expRawList)){
    expRaw <- expRawList[[ds]]
    
    # add Entrez gene ID as required by AIMS
    expRaw <- merge(id2symbol, expRaw, by = 0)
    expRaw <- expRaw[!duplicated(expRaw$gene_id),]
    rownames(expRaw) <- expRaw$gene_id
    expRaw <- expRaw[, 4:ncol(expRaw)]
    
    dataForAIMS <- list()
    dataForAIMS[['D']] <- as.matrix(expRaw)
    dataForAIMS[['EntrezID']] <- rownames(expRaw)
    aimsRes <- applyAIMS(dataForAIMS$D, dataForAIMS$EntrezID)
    
    aimsSubtype <- as.data.frame(aimsRes$cl)
    colnames(aimsSubtype) <- "Subtype"
    aimsSubtype$Subtype <- factor(aimsSubtype$Subtype, levels = c("LumA", "LumB", "Her2", "Normal", "Basal"))
    
    subtypeList[[ds]] <- aimsSubtype
}


# Since clinical variables not present in testing ds (expcept for HTG), so train a new model without them
train.paloma <- dtaSubset.final[, !colnames(dtaSubset.final) %in% c('cdk_prior_et', 'Organ')]

# repeat 5 times and get the best model
bestModel.paloma <- NULL
min.cindex <- 1

for (i in 1:5){
  o7 <- tune(Surv(pfs_time, true.progression) ~ ., train.paloma)
  
  obj.src.paloma <- rfsrc(Surv(pfs_time, true.progression) ~ ., train.paloma, ntree = 3000, 
                          importance = TRUE, mtry = o7[[2]][[2]], node.size = o7[[2]][[1]],
                          na.action = "na.impute")
  print(obj.src.paloma)
  
  # get OOB perforance error for best model selection
  cindex <- get.cindex(time = train.paloma$pfs_time, censoring = train.paloma$true.progression,
                       predicted = obj.src.paloma$predicted.oob)
  #rankings <- c(rankings, cindex)
  if (cindex < min.cindex){
    min.cindex <- cindex
    bestModel.paloma <- obj.src.paloma
  }
}

for (ds in names(subtypeList)){
    subtype.df <- subtypeList[[ds]]
    
    test.df <- expNormList[[ds]]
    
    test.df <- test.df[rownames(test.df) %in% bestModel.paloma$xvar.names, ]
    test.df <- as.data.frame(t(scale(test.df)))
    
    if (ds == "PALOMA-2"){
        test.df$drug <- 'AI'
        ds.predict <- predict(bestModel.paloma, test.df, na.action = 'na.impute')
    } else if (ds %in% c("PALOMA-3", "PEARL")) {
        test.df$drug <- 'Fulvestrant'
        ds.predict <- predict(bestModel.paloma, test.df, na.action = 'na.impute')
    } else if (ds == "Metabric"){
        test.df$drug <- 'AI'
        ds.predict <- predict(obj.src.metabric, test.df, na.action = 'na.impute')
    } else if (ds == "HTG"){
        if (!all(rownames(dtaSubset.cli2) == rownames(test.df))){
          stop("dataset missmatch")
        }
        
        subtype.df <- subtype.df %>% filter(rownames(subtype.df) %in% rownames(test.df))
        
        if(!all(rownames(subtype.df) == rownames(test.df))){
          stop("dataset mismatch")
        }
        
        test.df$drug <- dtaSubset.cli2$drug
        test.df$Organ <- dtaSubset.cli2$Organ
        test.df$cdk_prior_et <- dtaSubset.cli2$cdk_prior_et
        ds.predict <- predict(obj.src4, test.df, na.action = 'na.impute')
    } else {
        stop('Unknown dataset')
    }
    
    Mortality <- ds.predict$predicted
    
    subtype.mortality <- cbind(subtype.df, Mortality)
    subtype.mortality <- subtype.mortality[order(subtype.mortality$Mortality), ]
    subtype.mortality$Order <- 1:nrow(subtype.mortality)
    
    # increment is 1/number of each subtype, and decreament is 1/number of other subtype(s)
    incrementList <- list()
    decreamentList <- list()
    
    for (subtype in levels(subtype.df$Subtype)){
        incrementList[[subtype]] <- 1/sum(subtype.mortality$Subtype == subtype)
        decreamentList[[subtype]] <- 1/sum(subtype.mortality$Subtype != subtype)
        
        subtype.mortality[[subtype]] <- 0
        subtype.mortality[[subtype]][1] <- ifelse(subtype.mortality$Subtype[1] == subtype, 
                                                  incrementList[[subtype]], -decreamentList[[subtype]])
    }
    
    for (subtype in levels(subtype.df$Subtype)){
        for (i in 2:nrow(subtype.mortality)){
            subtype.mortality[[subtype]][i] <- ifelse(subtype.mortality$Subtype[i] == subtype, 
                                                      subtype.mortality[[subtype]][i-1] + incrementList[[subtype]],
                                                      subtype.mortality[[subtype]][i-1] - decreamentList[[subtype]])
        }
    }
    
    # add a data column for drawing a base line at zero
    subtype.mortality$baseline <- 0
    
    enrichment.score <- list()
    
    num.simulation <- 1000
    
    for (subtype in levels(subtype.df$Subtype)){
        
        simu.res <- vector(length = num.simulation)
        
        for (j in 1:num.simulation){
            tmp <- vector(length = nrow(subtype.mortality))
            
            # randomly select the same number of subtypes
            subtype.random <- sample(1:nrow(subtype.mortality), sum(subtype.mortality$Subtype == subtype))
            tmp[1] <- ifelse(i %in% subtype.random, 
                             incrementList[[subtype]], -decreamentList[[subtype]])
            for (k in 2:nrow(subtype.mortality)){
                tmp[k] <- ifelse(k %in% subtype.random, tmp[k-1] + incrementList[[subtype]],
                                 tmp[k-1] - decreamentList[[subtype]])
            }
            
            simu.res[j] <- max(abs(tmp))
        }
        
        enrichment.score[[subtype]][['ES']] <- signif((max(abs(subtype.mortality[[subtype]]))), 2) 
        pvalue <- sum(simu.res >= (max(abs(subtype.mortality[[subtype]]))))/num.simulation
        padj <- p.adjust(pvalue, method = 'fdr', n = num.simulation)
        enrichment.score[[subtype]][['Padj']] <- signif(padj, 2)
        
        # Define normalized enrichment score as log2 ES to mean ES of 1000x simulations ratio
        enrichment.score[[subtype]][['NES']] <- signif(log2(max(abs(subtype.mortality[[subtype]]))/mean(simu.res)), 2)
        
    }
    
    
    subtype.col <- colors[1:length(levels(subtype.df$Subtype))]
    names(subtype.col) <- levels(subtype.df$Subtype)
    
    ha = HeatmapAnnotation(
        Scores = anno_lines(subtype.mortality[, 4:ncol(subtype.mortality)], height = unit(10, 'cm'),
                            gp = gpar(col = c(subtype.col, 'black'), 
                                      lty =  c(rep(1, length(levels(subtype.mortality$Subtype))), 2), lwd = 1.5), add_points = FALSE),
        Subtype = subtype.mortality$Subtype, col = list(Subtype = subtype.col),
        annotation_legend_param = list(),
        Mortality = anno_barplot(subtype.mortality$Mortality, height = unit(2, 'cm'), 
                                 border = FALSE, gp = gpar(fill = 'gray50', col = 'gray50')))
    
    dummy.data <- matrix(rep(0, nrow(subtype.mortality)), 1)
    
    ht <- Heatmap(dummy.data,
                  show_heatmap_legend = FALSE,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  col = 'white',
                  top_annotation = ha)
    
    # add enrichment scores and p values
    ens.label <- 'NES (Padj):'
    for (subtype in names(enrichment.score)){
        tmp <- paste0('\t', subtype, ': ', enrichment.score[[subtype]][['NES']], ' (', enrichment.score[[subtype]][['Padj']], ')')
        ens.label <- paste(ens.label, tmp, sep = '\n')
    }
    
    png(paste(date, project, ds, "SubtypeEnrichmentAnalysis.png", sep = '_'),
        width = 7,
        height = 6,
        units = 'in',
        res = 150)
    
    draw(ht)
    decorate_annotation("Scores", {
        grid.text(ens.label, x = unit(0.5, "cm"), y = unit(1.8, "cm"), 
                  just= "left")
    })
    dev.off()
}

# Plot OOR performance error for each marker
stat.files <- list.files('.', pattern = 'stats.txt')
stat.files <- stat.files[!grepl('ElasticNetSignature', stat.files)]
data.for.barplot <- NULL

for (file in stat.files){
  
  tmp <- read.csv(file, sep = ':')
  rownames(tmp) <- sub('^\\s+', '', rownames(tmp))
  colnames(tmp) <- sub('X\\.1\\.\\.', '', colnames(tmp))
  tmp[[1]] <- sub('\\s+', '', tmp[[1]])
  
  tmp2 <- data.frame(OOB_error = tmp['(OOB) Requested performance error', ], Marker = colnames(tmp)[1])
  if (is.null(data.for.barplot)){
    data.for.barplot <- tmp2
  } else {
  data.for.barplot <- rbind(data.for.barplot, tmp2)
  }
}

#============================================================================================
# Plot predicted mortality for PALOMA2/3 and PEARL and test for concordance
#============================================================================================
paloma2.exp <- expNormP2[rownames(expNormP2) %in% bestModel.paloma$xvar.names, ]
paloma2.exp <- as.data.frame(t(scale(paloma2.exp)))
paloma2.exp$drug <- 'AI'
rownames(paloma2.exp) <- paste('PALOMA2', rownames(paloma2.exp))


# PALOMA3
paloma3.exp <- expNormP3[rownames(expNormP3) %in% bestModel.paloma$xvar.names, ]
paloma3.exp <- as.data.frame(t(scale(paloma3.exp)))
paloma3.exp$drug <- 'Fulvestrant'
rownames(paloma3.exp) <- paste('PALOMA3', rownames(paloma3.exp))


# PEARL
pearl.exp <- expNormPl[rownames(expNormPl) %in% bestModel.paloma$xvar.names, ]
pearl.exp <- as.data.frame(t(scale(pearl.exp)))
pearl.exp$drug <- 'Fulvestrant'
rownames(pearl.exp) <- paste('PEARL', rownames(pearl.exp))

# boxplot to show difference predicted mortality
externalTesting <- rbind(paloma2.exp, paloma3.exp, pearl.exp)
external.predict <- predict(bestModel.paloma, externalTesting, na.action = 'na.impute')
palPearlMort <- external.predict$predicted

data.to.plot <- data.frame(Mortality =  palPearlMort,
                           Source = sub(' .*$', '', rownames(external.predict$xvar)))

data.to.plot$Source <- factor(data.to.plot$Source, levels = c('PALOMA2', 'PALOMA3', 'PEARL'))

my_comparisons <- list( c("PALOMA2", "PALOMA3"), c("PALOMA2", "PEARL"), c("PALOMA3", "PEARL") )
p <- ggboxplot(data.to.plot, x = "Source", y = "Mortality",
               color = "Source", 
               palette = "jco",
               
               add = "jitter") + stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(label.y = 50)  +
    theme(legend.position="none", axis.title.x = element_blank()) + scale_color_manual(values=c('PALOMA2' = 'red', 'PALOMA3' = 'green', PEARL = 'blue'))

png(paste(date, project, "PredictedMortalityPalomaPearlBoxplot.png", sep = '_'),
    width = 4,
    height = 5,
    units = 'in',
    res = 200
    )
print(p)
dev.off()

# Plot predicted survival functions for the three dataset together
externalTesting <- rbind(paloma2.exp, paloma3.exp, pearl.exp)
external.predict <- predict(bestModel.paloma, externalTesting, na.action = 'na.impute')

pred.surv = t(external.predict$survival)
colnames(pred.surv) <- rownames(external.predict$xvar)
pred.surv <- as.data.frame(pred.surv)
pred.surv$Time <- external.predict$time.interest

dd = melt(pred.surv, id=c("Time"))
colnames(dd) <- c('Time', 'PID', 'Survival')

dd$PID <- as.character(dd$PID)

dd$Study <- 'PEARL'
dd$Study[dd$PID %in% rownames(paloma2.exp)] <- 'PALOMA2'
dd$Study[dd$PID %in% rownames(paloma3.exp)] <- 'PALOMA3'


p <- ggplot() + geom_line(data = dd, aes(x=Time, y=Survival, color = Study, group = PID, alpha = Study), linetype = 2) +
  labs(title="RSF Predicted Survival on PALOMA2/3 and PEARL Study") +
  theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12),
        legend.position = c(0.8, 0.8)) + scale_alpha_manual(values=c(1, 0.2, 0.1)) + 
  xlab('Time (Months)') + ylab('PFS Probability')

png(paste(date, project, "PredictedSurvivalFuncitonForPaloma2.png", sep = '_'),
    width = 5.5, height = 5.5, units = 'in', res = 150)
print(p)
dev.off()




#============================================================================================
# Plot OOB Requested performance error for each iteration
#============================================================================================
data.for.barplot$OOB_error <- as.numeric(data.for.barplot$OOB_error)

# average values for each round
data.for.barplot <- data.for.barplot %>%  group_by(Marker) %>%  summarise(OOB_error = mean(OOB_error))

data.for.barplot$Marker[data.for.barplot$Marker == 'BestReduced'] <- 'Final'
data.for.barplot$Iteration <- "FirstRounds"
data.for.barplot$Iteration[grepl('Round|Final', data.for.barplot$Marker)] <- 'LaterRounds'
data.for.barplot$Marker <- factor(data.for.barplot$Marker, levels = c("BCI", "CDK2", "CDK46", "EndoPredict",
                                                                      "G2M_Checkpoint", "IHC4", "Mammaprint",  
                                                                      "Mammostrat", "OncotypeDX","PALOMA23",
                                                                      "PALOMA3","PAM50", "PEARL", "Proliferation", "RBSig",       
                                                                      "Round2", "Round3", "Round4", "Final"))

p <- ggplot(data.for.barplot, aes(x = Marker, y = OOB_error, fill = Iteration)) + geom_bar(width = 0.5, stat="identity")
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + ylab('(OOB) Requsted performance error')
png(paste(date, project, "OOB_error_over_iteration.png", sep = '_'),
    height = 4,
    width = 5,
    units = 'in',
    res = 200)
print(p)
dev.off()

# save workspace
save.image(paste(date, project, "Workspace.RData", sep = '_'))

# Need to fix the renaming probelm for rerun later on with the May 07 data
#bestModelList.r1 <- bestModels.r1

#### WRITE SESSION INFO TO FILE ####################################################################
writeLines(capture.output(sessionInfo()), paste(date, project, "SessionInfo.txt", sep = '_'))

