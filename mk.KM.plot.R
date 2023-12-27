
# make a Kaplan-Meier plot
make.KM.plot <- function(survDf,  grouping.var = 'Group', surv.type = 'PFS', event = 'Status',  time.unit = 'Months', label.x = NULL, label.y = NULL, ...) {
  library("survival")
  library("survminer")
  
  fmla <- as.formula(paste('Surv(', surv.type, ', ', event, ') ~ ', grouping.var))
  survfit <- survfit(fmla, data = survDf)
  
  medianLabel <- paste0("\t", gsub('.*?=', '', rownames(median(survfit))), ': ', format(round(median(survfit), 2), nsmall = 2), collapse = "\n")
  medianLabel <- gsub('NA', 'Not reached', medianLabel)
  
  # get HR to put on K-M plot
  coxfit <- NULL
  
  tryCatch({
    coxfit <- coxph(fmla, data = survDf, ties = 'exact')
  },
  error = function(e) 
    stop(paste("Error: ", e))
  )
  
  x <- summary(coxfit)
  pvalue<-signif(x$sctest["pvalue"], digits=2)
  
  coefMatrix <- matrix(x$coefficients, ncol = 5)
  
  # Identify the most significant contrast
  ind <- which(coefMatrix[, 5] == min(coefMatrix[, 5]))
  
  hr <- signif(x$conf.int[ind, 'exp(coef)'], 2) 
  hr.low95 <- signif(x$conf.int[ind,"lower .95"], 2)
  hr.high95 <- signif(x$conf.int[ind,"upper .95"], 2)
  
  # Fit survival data using the Kaplan-Meier method
  survfit$call$formula <- fmla
  
  palette <- c("blue",  "orange", "red",   "purple", "coral1", "brown")
  
  if (length(unique(survDf[[grouping.var]])) == 2){
    palette <- c("blue", "red")
  } else if (length(unique(survDf[[grouping.var]])) == 3){
    palette <- c("blue",  "orange", "red")
  } else if (length(unique(survDf[[grouping.var]])) == 4){
    palette <- c("blue", "mediumturquoise", "orange", "red")
  }
  
  p <- ggsurvplot(survfit, data = survDf, 
                  pval = FALSE,
                  legend.title = '',
                  risk.table=T, 
                  fontsize=4,
                  font.tickslab = 14,
                  font.legend=13, 
                  palette=palette[1:length(unique(survDf[[grouping.var]]))],
                  xlab = time.unit, 
                  ylab = paste(surv.type, "Probability")
  )
  
  if (is.null(label.x)){
    label.x <- max(survDf[[surv.type]], na.rm = TRUE) * 0.5
  }
  if (is.null(label.y)){
    label.y <- 0.96
  }
  p$plot <- p$plot + ggplot2::annotate(
    "text",
    x = label.x, y = label.y,
    vjust  = 0.8, hjust = 0,
    
    label = paste0('HR = ', hr, ' (', hr.low95, '-', hr.high95, ')\nP = ', pvalue, ' (Logrank test)\nMedian ', surv.type, ':\n', medianLabel),
    size = 4
  ) + theme(legend.position = "none")     
  p
}

