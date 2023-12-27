library(AIMS)
library(org.Hs.eg.db)
library(randomForestSRC)
library(survival)
library(survminer)
library("ggplot2")
library("tidyr")
library("reshape2")
library(ggpubr)
library(GEOquery)
library(ComplexHeatmap)
library(dplyr)

date <- Sys.Date()
project <- 'RandomForestSRC'

setwd("/mnt/Jason/HTG/data")

source('/mnt/Jason/HTG/scripts/mk.KM.plot.R')

dta <- read.csv("all.norm.htg+clinical.batch.1-5.edgeR.csv", stringsAsFactors = TRUE)

rownames(dta) <- paste(dta$pcode, dta$Block, dta$Timepoint, dta$Organ, sep = '|')
dta <- dta[dta$drug %in% c('AI', 'Fulvestrant'), ]

# only use pretreatment samples
dtaPretx <- dta[grepl("2", dta$Timepoint), ]

# remove rows with more than one biopsies from a single patient
dtaPretx <- dtaPretx[!duplicated(dtaPretx$pcode), ]



set.seed(420)

colors <- c('#ffe119', '#3cb44b', '#e6194b', '#4363d8', '#f58231', 
            '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
            '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
            '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
            '#ffffff', '#000000')

resList <- list()

# use published biomarkers for feature selection
biomarker <- list(
    CDK2 = "../data/biomarkers/CDK2_signature.txt",
    OncotypeDX = "../data/biomarkers/OncotypeDX.txt",
    IHC4 = "../data/biomarkers/IHC4.txt",
    RBSig =  "../data/biomarkers/RBsig.txt",
    PAM50 = "../data/biomarkers/PAM50.txt",
    EndoPredict = "../data/biomarkers/EndoPredict.txt",
    Mammaprint = "../data/biomarkers/Mammaprint.txt",
    Mammostrat = "../data/biomarkers/Mammostrat.txt",
    Proliferation = "../data/biomarkers/ProliferationSignatureByCharlesPerou.txt",
    G2M_Checkpoint = "../data/biomarkers/G2M_Checkpoint.txt",
    ElasticNetSignature = "../data/biomarkers/ElasticNetSignature.txt",
    CDK46 = "../data/biomarkers/CDK4_6_GeneList.txt",
    BCI = "../data/biomarkers/BCI.txt"
)

biomarkerList <- list()

# read biomarkers
for (name in names(biomarker)){
    markerDf <- read.csv(biomarker[[name]], header = FALSE, stringsAsFactors = FALSE, sep = '\t')
    biomarkerList[[name]] <- markerDf$V1
}

vimp <- NULL

setwd("../output/RandomForests")

##=======================================================================================
## Round 1, initial biomarker evaluation
##=======================================================================================
# loop through each biomarker and perform feature selection
for (marker in names(biomarkerList)){
    genes.of.interest <- biomarkerList[[marker]]
    
    idx <- which(colnames(dtaPretx) %in% genes.of.interest)
    
    # scale the expression matrix
    tmp <- dtaPretx[, idx]
    tmp <- t(scale(t(tmp)))
    tmp <- as.data.frame(tmp)
    
    dtaSubset <- cbind(dtaPretx[, c(2,4)], tmp)
    
    # tune mtry and node size
    o1 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset)
    
    obj.src <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset, ntree = 1000, importance = TRUE, mtry = o1[[2]][[2]], node.size = o1[[2]][[1]], na.action = "na.impute")
    sink(file = paste(date, project, 'myoutput1.txt', sep = '_'), append = TRUE, type = 'output')
    print(marker)
    print(obj.src)
    sink()
    
    # get variable importance
    print(predict(obj.src, importance = TRUE)$importance)
    
    # plot OOB error rate and variable importance
    png(paste(date, project, marker, "VIM.png", sep = '_'),
        width = 8,
        height = 20,
        units = 'in',
        res = 200)
    plot(obj.src)
    
    dev.off()
    
    
    png(paste(date, project, marker, "Survival.png", sep = '_'),
        width = 10,
        height = 10,
        units = 'in',
        res = 200)
    plot.survival(rfsrc(Surv(pfs_time, true.progression)~ ., dtaSubset), cens.model = "rfsrc", collapse = F)
    dev.off()
    
    
    # variable confidence interval
    jk.obj <- subsample(obj.src)
    
    # adjust plot height dynamicly
    height <- length(idx) * 0.22 + 1.3
    
    png(paste(date, project, marker, "Round1_VimpConfidenceInterval.png", sep = '_'), width = 5, height = height, units = 'in', res = 100)
    par(oma=c(1,1,1,1))
    par(mar=c(4,6,0,0))
    plot(jk.obj, cex = 0.6)
    dev.off()
    
    
    # VIMP
    vimpbytreatment <- as.data.frame(vimp(obj.src, block.size = 10, joint = F)$importance)
    colnames(vimpbytreatment) <- "VIMP"
    vimpbytreatment <- vimpbytreatment[order(vimpbytreatment$VIMP, decreasing = TRUE), , drop = FALSE]
    
    # save to file
    write.table(vimpbytreatment, 
                file = paste(date, marker, "Vimp.txt", sep = '_'),
                quote = FALSE,
                col.names = F,
                sep = "\t"
    )
    
    #
    # Test model performance with K-M plot
    #
    
    # predict mortality on training dataset
    obj.predict <- predict(obj.src, dtaSubset)
    
    sink(file = paste(date, project, 'myoutput2.txt', sep = '_'), append = TRUE, type = 'output')
    print(marker)
    print(obj.predict)
    sink()
    
    yhat <- obj.predict$predicted
    
    q.yhat <- quantile(yhat, probs = c(0, 0.5, 1))
    q.yhat[1] <- q.yhat[1] - .0001
    Group <- cut(yhat, q.yhat, labels = FALSE)
    
    survDf <- cbind(obj.predict$yvar, Group)
    survDf$Group[survDf$Group == 1] <- 'LowRisk'
    survDf$Group[survDf$Group == 2] <- 'HighRisk'
    survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))
    
    p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time', event = 'true.progression')
    
    png(paste(date, project, marker, "Round1_VerificationOfPublishedBioMarkerOnHtgPreTx_K-M_plot.png", sep = '_'),
        width = 6, 
        height = 6, 
        units = "in", 
        pointsize = 10, 
        res = 300
    )
    print(p)
    dev.off()
    
    # save variables passing VIMP > 0.002 for next round
    tmp <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]
    vimp <- c(vimp, rownames(tmp))
}

##=======================================================================================
## Round 2, combined the important genes from round1 and repeat random survival modeling
## to further purify important genes
##=======================================================================================
important.genes <- unique(vimp)

length(important.genes) 

# tune mtry and node size
idx2 <- which(colnames(dtaPretx) %in% important.genes)

tmp <- dtaPretx[, idx2]
tmp <- t(scale(t(tmp)))
tmp <- as.data.frame(tmp)
dtaSubset2 <- cbind(dtaPretx[, c(2,4)], tmp)

o2 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset2)

obj.src2 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset2, ntree = 1000, 
                  importance = TRUE, mtry = o2[[2]][[2]], node.size = o2[[2]][[1]], 
                  na.action = "na.impute")
print(obj.src2)

# VIMP
vimpbytreatment <- as.data.frame(vimp(obj.src2, block.size = 10, joint = F)$importance)
colnames(vimpbytreatment) <- "VIMP"

# get variables passing VIMP > 0.002 for next round
vimp <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]

dim(vimp) # 56 1

# variable confidence interval
jk.obj2 <- subsample(obj.src2)

height <- length(idx2) * 0.16 + 1.3
png(paste(date, project, "Round2_FilteredFeatureVimpConfidenceInterval.png", sep = '_'), width = 5, height = height, units = 'in', res = 100)
par(oma=c(1,1,1,1))
par(mar=c(4,6,0,0))
plot(jk.obj2, cex = 1.2)
dev.off()

important.genes2 <- rownames(vimp)
length(important.genes2)
#
# Test model performance with K-M plot
#

# predict mortality on training dataset
obj.predict2 <- predict(obj.src2, dtaSubset2)

print(obj.predict2)
yhat <- obj.predict2$predicted

q.yhat <- quantile(yhat, probs = c(0, 0.5, 1))
q.yhat[1] <- q.yhat[1] - .0001
Group <- cut(yhat, q.yhat, labels = FALSE)

survDf <- cbind(obj.predict2$yvar, Group)
survDf$Group[survDf$Group == 1] <- 'LowRisk'
survDf$Group[survDf$Group == 2] <- 'HighRisk'
survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))

p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time', event = 'true.progression')

png(paste(date, project, "Round2_VerificationOfFilteredPublishedBioMarkerOnHtgPreTx_K-M_plot.png", sep = '_'),
    width = 6, 
    height = 6, 
    units = "in", 
    pointsize = 10, 
    res = 300
)
print(p)
dev.off()

##=======================================================================================
## Round 3, add clinical variables to the important genes selected from round2 and repeat
## random survival modeling to further purify important genes and clinical variables
##=======================================================================================
dtaSubset.cli <- cbind(dtaPretx[, c(2, 4, 6:12, 15)], tmp)
dtaSubset.cli$Organ <- as.character(dtaSubset.cli$Organ)
dtaSubset.cli$Organ[!dtaSubset.cli$Organ %in% c('Bone', 'Lung', 'Breast', 'Liver', 'Skin')] <- 'Other'
dtaSubset.cli$Organ <- as.factor(dtaSubset.cli$Organ)

o3 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli)

obj.src3 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli, ntree = 1000, 
                  importance = TRUE, mtry = o3[[2]][[2]], node.size = o3[[2]][[1]], 
                  na.action = "na.impute")
print(obj.src3)

# VIMP
vimpbytreatment <- as.data.frame(vimp(obj.src3, block.size = 10, joint = F)$importance)
colnames(vimpbytreatment) <- "VIMP"

#
# Test model performance with K-M plot
#

# predict mortality on training dataset
obj.predict3 <- predict(obj.src3, dtaSubset.cli)

print(obj.predict)
yhat <- obj.predict$predicted

q.yhat <- quantile(yhat, probs = c(0, 0.5, 1))
q.yhat[1] <- q.yhat[1] - .0001
Group <- cut(yhat, q.yhat, labels = FALSE)

survDf <- cbind(obj.predict$yvar, Group)
survDf$Group[survDf$Group == 1] <- 'LowRisk'
survDf$Group[survDf$Group == 2] <- 'HighRisk'
survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))

p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time', event = 'true.progression')

png(paste(date, project, "Round3_VerificationOfFilteredPublishedBioMarkerPlusClinicalCovariatesOnHtgPreTx_K-M_plot.png", sep = '_'),
    width = 6, 
    height = 6, 
    units = "in", 
    pointsize = 10, 
    res = 300
)
print(p)
dev.off()

# get variables passing VIMP > 0.002 for next round use
vimp <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]
dim(vimp)

rownames(vimp)
# variable confidence interval
jk.obj3 <- subsample(obj.src3)
height <- length(idx2) * 0.16 + 1.3
png(paste(date, project, "Round3_FilteredFeatureVimpConfidenceInterval.png", sep = '_'), width = 5, height = height, units = 'in', res = 100)
par(oma=c(1,1,1,1))
par(mar=c(4,6,0,0))
plot(jk.obj3, cex = 1.2)
dev.off()

##=======================================================================================
## Round 4, repeat random survival modeling to verify features selected
##=======================================================================================
dtaSubset.cli2 <- cbind(dtaSubset.cli[, 1:2], dtaSubset.cli[, rownames(vimp)])

# due to the algorithmic nature of random forest, each time a tree is built, a slightly
# different vimp list is returned. If multiple features have similar resolving power, they
# may have different vimp values. However, the strongest feature should always the same.
# We repeat this process 100 times and then tabulate the frequencies of each feature being
# in the final vimp list
for (i in 1:2){
    o4 <- tune(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli2)
    
    obj.src4 <- rfsrc(Surv(pfs_time, true.progression) ~ ., dtaSubset.cli2, ntree = 1000, 
                      importance = TRUE, mtry = o4[[2]][[2]], node.size = o4[[2]][[1]], 
                      na.action = "na.impute")
    print(obj.src4)
    
    # VIMP
    vimpbytreatment <- as.data.frame(vimp(obj.src4, block.size = 10, joint = F)$importance)
    colnames(vimpbytreatment) <- "VIMP"
    
    # save variables passing VIMP > 0.002. This should have no effect
    vimp <- vimpbytreatment[vimpbytreatment$VIMP > 0.002, , drop = FALSE]
    dim(vimp)
    final.features <- rownames(vimp)
    
    # variable confidence interval
    jk.obj <- subsample(obj.src4)
    height <- nrow(vimp) * 0.24 + 1.3
    png(paste(date, project, i, "Round4_DoublelyFilteredFeatureVimpConfidenceInterval.png", sep = '_'), width = 5, height = height, units = 'in', res = 100)
    par(oma=c(1,1,1,1))
    par(mar=c(4,6,2,2) + 0.1)
    plot(jk.obj, cex = 1.2)
    dev.off()
    
    # save to file
    write.table(x = vimp,
                file = paste(date, project, i, 'FinalFeaturesFromRound4.txt', sep = '_'),
                quote = FALSE,
                sep = "\t",
                row.names = TRUE)
    
}

#
# Looks good and use obj.src4 as the final model for prediction
#

# save the final model
save(obj.src4, file = 'FinalRSFModel.Rdata')

##===========================================================================================
## verify on the same(training) dataset
##===========================================================================================

# predict mortality on training dataset
obj.predict4 <- predict(obj.src4, dtaSubset.cli2)

print(obj.predict4)
yhat <- obj.predict4$predicted

q.yhat <- quantile(yhat, probs = c(0, 0.5, 1))
q.yhat[1] <- q.yhat[1] - .0001
Group <- cut(yhat, q.yhat, labels = FALSE)

survDf <- cbind(obj.predict4$yvar, Group)
survDf$Group[survDf$Group == 1] <- 'LowRisk'
survDf$Group[survDf$Group == 2] <- 'HighRisk'
survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))

p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time', event = 'true.progression')

png(paste(date, project, "TestOfFinalModelOnHtgPreTx_PFS_K-M_plot.png", sep = '_'),
    width = 6, 
    height = 6, 
    units = "in", 
    pointsize = 10, 
    res = 300
)
print(p)
dev.off()

# OS
survDf2 <- merge(survDf, dta[, c('os', 'precdk_os_time')], by = 0)
p <- make.KM.plot(survDf = survDf2, grouping.var = 'Group', surv.type = 'precdk_os_time', event = 'os')

png(paste(date, project, "TestOfFinalModelOnHtgPreTx_OS_K-M_plot.png", sep = '_'),
    width = 6, 
    height = 6, 
    units = "in", 
    pointsize = 10, 
    res = 300
)
print(p)
dev.off()

png(paste(date, project, "FinalFeaturesSurvival.png", sep = '_'),
    width = 10,
    height = 10,
    units = 'in',
    res = 200)
plot.survival(rfsrc(Surv(pfs_time, true.progression)~ ., dtaSubset.cli2), cens.model = "rfsrc", collapse = F)
dev.off()



# test final model on primary/recurrence samples
feature.to.use <- colnames(dtaSubset.cli2)
dataPriRecur <- dta %>% filter(Timepoint %in% c('1', 'R')) %>% dplyr::select(all_of(feature.to.use)) 

obj.predict5 <- predict(obj.src4, dataPriRecur)

print(obj.predict5)
yhat <- obj.predict5$predicted

q.yhat <- quantile(yhat, probs = c(0, 0.5, 1))
q.yhat[1] <- q.yhat[1] - .0001
Group <- cut(yhat, q.yhat, labels = FALSE)

survDf <- cbind(obj.predict5$yvar, Group)
survDf$Group[survDf$Group == 1] <- 'LowRisk'
survDf$Group[survDf$Group == 2] <- 'HighRisk'
survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))

# test 
dta.pri.rec.AI <- dta %>% filter(drug == 'AI', Timepoint %in% c('1', 'R')) %>% dplyr::select('pcode')
survDf.Lowrisk <- survDf %>% filter(Group == 'LowRisk')

dta.pri.rec.Ful <- dta %>% filter(drug == 'Fulvestrant', Timepoint %in% c('1', 'R')) %>% dplyr::select('pcode')
survDf.Highrisk <- survDf %>% filter(Group == 'HighRisk')

setdiff(rownames(dta.pri.rec.Ful), rownames(survDf.Highrisk))
setdiff(rownames(dta.pri.rec.AI), rownames(survDf.Lowrisk))

p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time', event = 'true.progression')

png(paste(date, project, "TestOfFinalModelOnHtgPriRecurr_PFS_K-M_plot.png", sep = '_'),
    width = 6, 
    height = 6, 
    units = "in", 
    pointsize = 10, 
    res = 300
)
print(p)
dev.off()

# OS
survDf2 <- merge(survDf, dta[, c('os', 'precdk_os_time')], by = 0)
p <- make.KM.plot(survDf = survDf2, grouping.var = 'Group', surv.type = 'precdk_os_time', event = 'os')

png(paste(date, project, "TestOfFinalModelOnHtgPriRecurr_OS_K-M_plot.png", sep = '_'),
    width = 6, 
    height = 6, 
    units = "in", 
    pointsize = 10, 
    res = 300
)
print(p)
dev.off()


#
# Test for performance in AI and Fulvestrant treated cohorts separately
#
for (drug in c('AI', 'Fulvestrant')){
    test.df <- dtaSubset.cli2[dtaSubset.cli2$drug == drug, ]
    obj.predict.drug <- predict(obj.src4, test.df, na.action = "na.impute")
    
    print(obj.predict.drug)
    yhat <- obj.predict.drug$predicted
    
    q.yhat <- quantile(yhat, probs = c(0, 0.5, 1))
    q.yhat[1] <- q.yhat[1] - .0001
    Group <- cut(yhat, q.yhat, labels = FALSE)
    
    survDf <- cbind(obj.predict.drug$yvar, Group)
    survDf$Group[survDf$Group == 1] <- 'LowRisk'
    survDf$Group[survDf$Group == 2] <- 'HighRisk'
    survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))
    
    p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time', event = 'true.progression')
    
    png(paste(date, project, drug, "TestOfFinalModelOnHtgPreTx_PFS_K-M_plot.png", sep = '_'),
        width = 6, 
        height = 6, 
        units = "in", 
        pointsize = 10, 
        res = 300
    )
    print(p)
    dev.off()
    
    # also test on os
    survDf2 <- merge(survDf, dta[, c('os', 'precdk_os_time')], by = 0)
    
    if (drug == 'AI'){
        p <- make.KM.plot(survDf = survDf2, grouping.var = 'Group', surv.type = 'precdk_os_time', event = 'os', label.x = 0, label.y = 0.25)
    } else {
        p <- make.KM.plot(survDf = survDf2, grouping.var = 'Group', surv.type = 'precdk_os_time', event = 'os')
    }
    
    png(paste(date, project, drug, "TestOfFinalModelOnHtgPreTxPreCDK_OS_K-M_plot.png", sep = '_'),
        width = 6, 
        height = 6, 
        units = "in", 
        pointsize = 10, 
        res = 300
    )
    print(p)
    dev.off()
    
    png(paste(date, project, drug, "FinalFeaturesSurvival.png", sep = '_'),
        width = 10,
        height = 10,
        units = 'in',
        res = 200)
    plot.survival(rfsrc(Surv(pfs_time, true.progression)~ ., test.df), cens.model = "rfsrc", collapse = F)
    dev.off()
}


# OS using primary/recurrence samples
for (drug in c('AI', 'Fulvestrant')){
    test.df <- dataPriRecur[dataPriRecur$drug == drug, ]
    obj.predict.drug <- predict(obj.src4, test.df, na.action = "na.impute")
    
    print(obj.predict.drug)
    yhat <- obj.predict.drug$predicted
    
    q.yhat <- quantile(yhat, probs = c(0, 0.5, 1))
    q.yhat[1] <- q.yhat[1] - .0001
    Group <- cut(yhat, q.yhat, labels = FALSE)
    
    survDf <- cbind(obj.predict.drug$yvar, Group)
    survDf$Group[survDf$Group == 1] <- 'LowRisk'
    survDf$Group[survDf$Group == 2] <- 'HighRisk'
    survDf$Group <- factor(survDf$Group, levels = c('LowRisk', 'HighRisk'))
    
    p <- make.KM.plot(survDf = survDf, grouping.var = 'Group', surv.type = 'pfs_time', event = 'true.progression')
    
    png(paste(date, project, drug, "TestOfFinalModelOnHtgPriRecurr_PFS_K-M_plot.png", sep = '_'),
        width = 6, 
        height = 6, 
        units = "in", 
        pointsize = 10, 
        res = 300
    )
    print(p)
    dev.off()
    
    # also test on os
    survDf2 <- merge(survDf, dta[, c('os', 'precdk_os_time')], by = 0)
    
    if (drug == 'AI'){
        p <- make.KM.plot(survDf = survDf2, grouping.var = 'Group', surv.type = 'precdk_os_time', event = 'os', label.x = 0, label.y = 0.25)
    } else {
        p <- make.KM.plot(survDf = survDf2, grouping.var = 'Group', surv.type = 'precdk_os_time', event = 'os')
    }
    
    png(paste(date, project, drug, "TestOfFinalModelOnHtgPriRecurr_OS_K-M_plot.png", sep = '_'),
        width = 6, 
        height = 6, 
        units = "in", 
        pointsize = 10, 
        res = 300
    )
    print(p)
    dev.off()
}


##=======================================================================================
## Scatter plot to show correlation of mortality with actual PFS time
##=======================================================================================
data.for.scatterplot <- data.frame(Mortality = obj.predict4$predicted, PFS_time = dtaSubset.cli2$pfs_time, 
                                   Status = ifelse(dtaSubset.cli2$true.progression == 0, 'Censored', 'Event'))

p <- ggscatter(data.for.scatterplot, x = "PFS_time", y = "Mortality",  
               color = "Status", shape = 21, size = 2, # Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "pearson", label.x = 20, label.sep = "\n")
)

png(paste(date, project, "PredictedMortality_vs_PFS_on_HTG_preTX.png", sep = '_'),
    width = 5,
    height = 5,
    units = 'in',
    res = 200
)
print(p)
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

save(obj.src5, file = 'RSFModelNoClinical.Rdata')

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

dd$Response <- 'Sensitive'
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

# use raw count of HTG data for subtype prediction
htgRaw <- read.csv('/mnt/Jason/HTG/data/all.raw.htg+clinical.batch.1-5.csv')

# Remove P6, bone Run1 as it was sent twice. Keep the second one. There is no clear reason to keep the second one
htgRaw <- htgRaw %>% filter(!(pcode == 'P6' & Organ == 'Bone'& batch == 'Run1'))

rownames(htgRaw) <- paste(htgRaw$pcode, htgRaw$Block, htgRaw$Timepoint, htgRaw$Organ, sep = '|')

htgRaw <- htgRaw[htgRaw$drug %in% c('AI', 'Fulvestrant'), ]

# only use pretreatment samples
htgRawPreTx <- htgRaw %>% filter(Timepoint == '2')

expRawList[['HTG']] <- as.data.frame(t(htgRawPreTx[, 17:ncol(htgRawPreTx)]))

expNormList <- list()
expNormList[['PALOMA-2']] <- expNormP2
expNormList[['PALOMA-3']] <- expNormP3
expNormList[['PEARL']] <- expNormPl
expNormList[['Metabric']] <- as.data.frame(t(metabricSubset[, 39:ncol(metabricSubset)]))
expNormList[['HTG']] <- as.data.frame(t(dtaPretx[, 18:ncol(dtaPretx)]))

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
    
    # remov NA columns
    #na.col <- apply(expRaw, 2, sum)
    #expRaw <- expRaw[, !is.na(na.col)]
    
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
train.paloma <- dtaSubset.cli2[, !colnames(dtaSubset.cli2) %in% clinical.variables]

o7 <- tune(Surv(pfs_time, true.progression) ~ ., train.paloma)

obj.src.paloma <- rfsrc(Surv(pfs_time, true.progression) ~ ., train.paloma, ntree = 3000, 
                         importance = TRUE, mtry = o7[[2]][[2]], node.size = o7[[2]][[1]],
                         na.action = "na.impute")
print(obj.src.paloma)


for (ds in names(subtypeList)){
    subtype.df <- subtypeList[[ds]]
    
    test.df <- expNormList[[ds]]
    
    test.df <- test.df[rownames(test.df) %in% obj.src.paloma$xvar.names, ]
    test.df <- as.data.frame(t(scale(test.df)))
    
    if (ds == "PALOMA-2"){
        test.df$drug <- 'AI'
        ds.predict <- predict(obj.src.paloma, test.df, na.action = 'na.impute')
    } else if (ds %in% c("PALOMA-3", "PEARL")) {
        test.df$drug <- 'Fulvestrant'
        ds.predict <- predict(obj.src.paloma, test.df, na.action = 'na.impute')
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


#### WRITE SESSION INFO TO FILE ####################################################################
writeLines(capture.output(sessionInfo()), paste(date, project, "SessionInfo.txt", sep = '_'))
