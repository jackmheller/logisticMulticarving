

library(BBmisc)
library(readxl)
library(stringr)

setwd('/Users/jackheller/KU/Thesis/Thesis/simulation_setups/multi_carve')

calculatePower <- function(subres, sparsity) {
  names<-c("carve5", "carvefw5", "split5", "splitfw5", "carve30",
           "carvefw30", "split30", "splitfw30", "carvefix5",
           "carvefwfix5", "splitfix5", "splitfwfix5", "carvefix30",
           "carvefwfix30", "splitfix30", "splitfwfix30")
  allrej <- matrix(NA, nrow = 1,ncol = length(names))
  fwer <- numeric()
  colnames(allrej) <- names
  for (name in names) {
    nameind <- which(colnames(subres) == name)
    mat <- subres[, nameind]
    fwer[as.character(name)] <- mean(mat[,sparsity + 1] < 0.05, na.rm = TRUE)
    rejmat <- mat[, 1:sparsity] < 0.05
    allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
  }
  return(list(allrej, fwer))
}

calculatePowerB1 <- function(subres, sparsity) {
  names<-c("carve", "carvefw", "split", "splitfw")
  allrej <- matrix(NA, nrow = 1,ncol = length(names))
  fwer <- numeric()
  colnames(allrej) <- names
  for (name in names) {
    nameind <- which(colnames(subres) == name)
    mat <- subres[, nameind]
    fwer[as.character(name)] <- mean(mat[,sparsity + 1] < 0.05, na.rm = TRUE)
    rejmat <- mat[, 1:sparsity] < 0.05
    allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
  }
  return(list(allrej, fwer))
}

calculatePScreen <- function(subres, sparsity) {
  names<-c("carve5", "carvefw5", "split5", "splitfw5", "carve30",
           "carvefw30", "split30", "splitfw30", "carvefix5",
           "carvefwfix5", "splitfix5", "splitfwfix5", "carvefix30",
           "carvefwfix30", "splitfix30", "splitfwfix30")
  allrej <- matrix(NA, nrow = 1,ncol = length(names))
  colnames(allrej) <- names
  for (name in names) {
    nameind <- which(colnames(subres) == name)
    mat <- subres[, nameind]
    
    binMat <- mat[, 1:sparsity] < 0.05
    pscreen <- mean(rowMeans(binMat) == 1)
    allrej[, as.character(name)] <- pscreen
  }
  return(allrej)
}

calculateExpectations <- function(subres, sparsity) {
  names<-c("carve5", "carvefw5", "split5", "splitfw5", "carve30",
           "carvefw30", "split30", "splitfw30", "carvefix5",
           "carvefwfix5", "splitfix5", "splitfwfix5", "carvefix30",
           "carvefwfix30", "splitfix30", "splitfwfix30")
  EV <- matrix(NA, nrow = 1, ncol = length(names))
  colnames(EV) <- names
  ERV <- matrix(NA, nrow = 1, ncol = length(names))
  colnames(ERV) <- names
  FDR <- matrix(NA, nrow = 1, ncol = length(names))
  colnames(FDR) <- names
  for (name in names) {
    nameind <- which(colnames(subres) == name)
    mat <- subres[, nameind]
    
    EV[, as.character(name)] <- mean(subres$`R-V`)
    ERV[, as.character(name)] <- mean(subres$V)
    FDR[, as.character(name)] <- mean(subres$`R-V` / pmax(subres$R, 1))
  }
  return(list(EV, ERV, FDR))
}

calculateExpectationsB1 <- function(subres, sparsity) {
  names<-c("carve", "carvefw", "split", "splitfw")
  EV <- matrix(NA, nrow = 1, ncol = length(names))
  colnames(EV) <- names
  ERV <- matrix(NA, nrow = 1, ncol = length(names))
  colnames(ERV) <- names
  FDR <- matrix(NA, nrow = 1, ncol = length(names))
  colnames(FDR) <- names
  for (name in names) {
    nameind <- which(colnames(subres) == name)
    mat <- subres[, nameind]
    
    EV[, as.character(name)] <- mean(subres$`R-V`)
    ERV[, as.character(name)] <- mean(subres$V)
    FDR[, as.character(name)] <- mean(subres$`R-V` / pmax(subres$R, 1))
  }
  return(list(EV, ERV, FDR))
}

simDf <- read_excel('simulations.xlsx')[c(1, 4, 5, 11:13),]
simDf$Time <-str_pad(simDf$Time, width=5, side="left", pad="0")

mainFunc <- function() {
  B.vec <- c(1, 5, 10, 20, 50) # number of splits
  frac.vec <- c(0.5, 0.75, 0.9, 0.95, 0.99) # selection fraction
  sparsity <- 5
  intCols <- c('carve5', 'split5', 'carve30', 'split30')
  outCols <- c('EV', 'ERV', 'FWER', 'FDR', 'Power')
  finMat <- matrix(NA, nrow = 0, ncol = length(outCols))
  colnames(finMat) <- outCols
  for (B in B.vec) {
    for (frac in frac.vec) {
      for (row in 1:nrow(simDf)) {
        curDat <- tryCatch(expr = {
          date <- format(as.Date(toString(simDf[row, 1]), format = "%d.%m.%y"), "%d-%b-%Y")
          time <- toString(simDf[row, 2])
          folder <- paste0('./Binomial_', date, ' ', time, '/')
          # file <- list.files(folder, pattern=paste0('results.* split=0.5 B=',B,' seed=.*'))
          file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
          curDat <- load2(paste0(folder, file))$results
          #return(curDat)
        },
        error = function(e) {
          date <- format(as.Date(toString(simDf[row, 1]), format = "%d.%m.%y"), "%d-%b-%Y")
          date <- gsub("Mar", "Mrz", date)
          time <- toString(simDf[row, 2])
          folder <- paste0('./Binomial_', date, ' ', time, '/')
          file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
          # file <- list.files(folder, pattern=paste0('results.* split=', frac, ' B=',B,' seed=.*'))
          curDat <- load2(paste0(folder, file))$results
          return(curDat)
        })
        
        if (B == 1) {
          retList <- calculatePowerB1(curDat, sparsity)
        }
        else {
          retList <- calculatePower(curDat, sparsity)
        }

        power <- retList[[1]]
        fwer <- retList[[2]]
        
        if (B == 1) {
          # create dummy matrix to load
          curMat <- matrix(NA, nrow = 2, ncol = ncol(finMat))
          rownames(curMat) <- c(paste0('carve_B=', B, 'select=', frac, 'n=', simDf[row,]$n, 'p=', simDf[row,]$p), 
                                paste0('split_B=', B, 'select=', frac, 'n=', simDf[row,]$n, 'p=', simDf[row,]$p))
          
          curMat[1, 5] <- power[2]
          curMat[2, 5] <- power[4]
          
          curMat[1, 3] <- fwer[2]
          curMat[2, 3] <- fwer[4]
        }
        else {
          # create dummy matrix to load
          curMat <- matrix(NA, nrow = 4, ncol = ncol(finMat))
          rownames(curMat) <- c(paste0('carve_gam5_B=', B, 'select=', frac, 'n=', simDf[row,]$n, 'p=', simDf[row,]$p), 
                                paste0('split_gam5_B=', B, 'select=', frac, 'n=', simDf[row,]$n, 'p=', simDf[row,]$p),
                                paste0('carve_gam30_B=', B, 'select=', frac, 'n=', simDf[row,]$n, 'p=', simDf[row,]$p), 
                                paste0('split_gam30_B=', B, 'select=', frac, 'n=', simDf[row,]$n, 'p=', simDf[row,]$p))
          # load power
          curMat[1, 5] <- power[2]
          curMat[2, 5] <- power[4]
          curMat[3, 5] <- power[6]
          curMat[4, 5] <- power[8]
          
          # load fwer
          curMat[1, 3] <- fwer[2]
          curMat[2, 3] <- fwer[4]
          curMat[3, 3] <- fwer[6]
          curMat[4, 3] <- fwer[8]
        }
        
        if (B == 1) {
          expectations <- calculateExpectationsB1(curDat, sparsity)
        }
        else {
          # capture EV and ERV
          expectations <- calculateExpectations(curDat, sparsity)
        }
        EV <- expectations[[1]]
        ERV <- expectations[[2]]
        FDR <- expectations[[3]]
        
        if (B == 1) {
          curMat[1, 1] <- EV[2]
          curMat[2, 1] <- EV[4]
          
          curMat[1, 2] <- ERV[2]
          curMat[2, 2] <- ERV[4]
          
          curMat[1, 4] <- FDR[2]
          curMat[2, 4] <- FDR[4]
        }
        else {
          curMat[1, 1] <- EV[2]
          curMat[2, 1] <- EV[4]
          curMat[3, 1] <- EV[6]
          curMat[4, 1] <- EV[8]
          
          curMat[1, 2] <- ERV[2]
          curMat[2, 2] <- ERV[4]
          curMat[3, 2] <- ERV[6]
          curMat[4, 2] <- ERV[8]
          
          curMat[1, 4] <- FDR[2]
          curMat[2, 4] <- FDR[4]
          curMat[3, 4] <- FDR[6]
          curMat[4, 4] <- FDR[8]
          
        }
        finMat <- rbind(finMat, curMat)
      }
    }
  }
  return(finMat)
}

finalMat <- mainFunc()

