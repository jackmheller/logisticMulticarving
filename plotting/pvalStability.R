

library(BBmisc)
library(readxl)
library(stringr)

setwd('~/KU/Thesis/Thesis')

# adjust depending on folder structure
source("inference/hdi_adjustments.R")
source("inference/carving.R")
source("inference/sample_from_truncated.R")
source("inference/tryCatch-W-E.R")

setwd('~/KU/Thesis/Thesis/simulation_setups/multi_carve')

simDf <- read_excel('simulations.xlsx')[7,]
simDf$Time <-str_pad(simDf$Time, width=5, side="left", pad="0")

# Calculate data across B
summStatsB <- function(data) {
  fullVar <- as.data.frame(matrix(NA, nrow = 0, ncol = dim(data[[1]])[2]))
  fullStd <- as.data.frame(matrix(NA, nrow = 0, ncol = dim(data[[1]])[2]))
  fullMean <- as.data.frame(matrix(NA, nrow = 0, ncol = dim(data[[1]])[2]))
  for (sim in data) {
    variance <- sapply(as.data.frame(sim), var)
    standardDev <- sapply(as.data.frame(sim), sd)
    meanVal <- sapply(as.data.frame(sim), mean)
    
    fullVar <- rbind(fullVar, variance)
    fullStd <- rbind(fullStd, standardDev)
    fullMean <- rbind(fullMean, meanVal)
    hist(as.data.frame(sim)[,1])
  }
  return(list(fullVar, fullStd, fullMean))
}

# Calculate data across Sims
summStatsSim <- function(data, B = 50) {
  aggregatedList <- as.data.frame(matrix(NA, nrow = 0, ncol = dim(data[[1]])[2]))
  for (sim in data) {
    pvals.aggregated <- pval.aggregator(list(sim),
                                        round(seq(ceiling(0.05 * B)/B, 1 - 1/B, by = 1/B), 2), cutoff = TRUE)
    aggregatedList <- rbind(aggregatedList, pvals.aggregated[[1]])
  }
  variance <- sapply(aggregatedList, var)
  standardDev <- sapply(aggregatedList, sd)
  meanVal <- sapply(aggregatedList, mean)
}

mainFunc <- function() {
  B.vec <- c(1, 5, 10, 20, 50) # number of splits
  frac.vec <- c(0.5, 0.75, 0.9, 0.95, 0.99) # selection fraction
  sparsity <- 5
  intCols <- c('carve5', 'split5', 'carve30', 'split30')
  outCols <- c('Variance', 'Mean')
  finMat <- matrix(NA, nrow = 0, ncol = length(outCols))
  colnames(finMat) <- outCols
  boxData <- as.data.frame(matrix(NA, nrow = 100, ncol = 0))
  PvalTable <- as.data.frame(matrix(NA, nrow = length(frac.vec), ncol = length(B.vec)))
  rownames(PvalTable) <- frac.vec
  colnames(PvalTable) <- B.vec
  for (frac in 0.75) {
    for (B in B.vec) {
      for (row in 1:nrow(simDf)) {
        # curDat <- tryCatch(expr = {
        #   date <- format(as.Date(toString(simDf[row, 1]), format = "%d.%m.%y"), "%d-%b-%Y")
        #   time <- toString(simDf[row, 2])
        #   folder <- paste0('./Binomial_', date, ' ', time, '/')
        #   # file <- list.files(folder, pattern=paste0('results.* split=0.5 B=',B,' seed=.*'))
        #   file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
        #   curDat <- list(load2(paste0(folder, file))$pvals.aggregated$`1_carveFWER`, load2(paste0(folder, file))$sel.index)
        #   return(curDat)
        # },
        # error = function(e) {
        #   date <- format(as.Date(toString(simDf[row, 1]), format = "%d.%m.%y"), "%d-%b-%Y")
        #   date <- gsub("Mar", "Mrz", date)
        #   time <- toString(simDf[row, 2])
        #   folder <- paste0('./Binomial_', date, ' ', time, '/')
        #   file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
        #   # file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
        #   curDat <- list(load2(paste0(folder, file))$pvals.aggregated$`1_carveFWER`, load2(paste0(folder, file))$sel.index)
        #   return(curDat)
        # })
        # 
        # selIndex <- curDat[[2]]
        # curDat <- curDat[[1]]
        # statsB <- summStatsB(curDat)
        # boxVar <- statsB[[1]][, selIndex]
        # boxStd <- statsB[[2]][, selIndex]
        # boxMean <- statsB[[3]][, selIndex]
        # 
        # boxplot(boxVar)
        # boxplot(boxStd)
        # boxplot(boxMean)
        # hist(boxVar[,1])
        
        curDat <- tryCatch(expr = {
          date <- format(as.Date(toString(simDf[row, 1]), format = "%d.%m.%y"), "%d-%b-%Y")
          time <- toString(simDf[row, 2])
          folder <- paste0('./Binomial_', date, ' ', time, '/')
          # file <- list.files(folder, pattern=paste0('results.* split=0.5 B=',B,' seed=.*'))
          file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
          curDat <- list(load2(paste0(folder, file))$results)
          return(curDat)
        },
        error = function(e) {
          date <- format(as.Date(toString(simDf[row, 1]), format = "%d.%m.%y"), "%d-%b-%Y")
          date <- gsub("Mar", "Mrz", date)
          time <- toString(simDf[row, 2])
          folder <- paste0('./Binomial_', date, ' ', time, '/')
          file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
          # file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
          curDat <- list(load2(paste0(folder, file))$results)
          return(curDat)
        })
        
        print(B)
        print(frac)
        print(colnames(curDat[[1]])[7])
        PvalTable[which(frac == frac.vec), which(B == B.vec)] <- sum(curDat[[1]][,8] < 0.05)
        # hist(curDat[[1]][,7], breaks = 20)
        browser()
        hist(curDat[[1]][,7], breaks = 20, main = paste('P values across Sims with B = ', B, ' Selec Frac = ', frac, ' ', expression(gamma), ' = 0.05'),
             cex = 2, xlab = "p-value (adjusted for multiplicity)", ylab = "Number of Simulations (out of 100)", cex.lab = 2, cex.main = 2, cex.axis = 2)
        boxplot(curDat[[1]][,7])
        if (frac == 0.9) {
          boxData[,as.character(B)] <- curDat[[1]][,9]
          # colnames(boxData)[dim(boxData)[2]] <- B
        }
        
        print(B)
        print(frac)
        print(colnames(curDat[[1]])[19])
        hist(curDat[[1]][,19], breaks = 20)
        
        # statsSim <- summStatsSim(curDat)
      }
    }
  }
  browser()
  
  boxplot(boxData)
}

finalMat <- mainFunc()


