

library(BBmisc)
library(readxl)
library(stringr)

setwd('/Users/jackheller/KU/Thesis/Thesis/simulation_setups/multi_carve')

calculatePowerAdj <- function(subres, sparsity) {
  names<-c("carve5", "carvefw5", "split5", "splitfw5", "carve30",
           "carvefw30", "split30", "splitfw30", "carvefix5",
           "carvefwfix5", "splitfix5", "splitfwfix5", "carvefix30",
           "carvefwfix30", "splitfix30", "splitfwfix30")
  allrej <- matrix(NA, nrow = 1,ncol = length(names))
  colnames(allrej) <- names
  for (name in names) {
    nameind <- which(colnames(subres) == name)
    mat <- subres[, nameind]
    rej <- quantile(mat[, sparsity + 1], 0.05, na.rm = TRUE)
    rejmat <- mat[, 1:sparsity] < rej
    allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
  }
  return(allrej)
}

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

plotVals <- function(df, B, frac, fwer = NA) {
  browser()
  legendNames <- c("Beta = 1", "Beta = 5", "Beta = 2", "Beta = 7.5", "Beta = 0.5")
  mypal <- colorRampPalette(c("black", "red", "green", "blue", "pink") )(5)
  # jpeg(file=paste0("Power with B = ", B, " Selection Frac = ", frac, ".jpeg"), height = 10, width = 10, units = 'in', res = 300)
  plot(jitter(as.numeric(rownames(df)), amount = 0.01), df[,1], type = "o", col = mypal[1], ylim = c(0, 1), xlim = c(0, 1),
       xlab = "Selection Fraction", ylab = "Power (-) and FWER (o)", main = paste0("Power of Carving (gamma = 0.05) with B = ", B), 
       lty = 1, pch = (1 + 14), lwd = 3, cex.lab = 1.5, cex.axis = 1.5)
  for (col in 2:ncol(df)) {
    lines(jitter(as.numeric(rownames(df)), amount = 0.01), df[,col], type = "b", col = mypal[col], lty = col, pch = (col + 14),
          lwd = 3)
  }
  if (!is.na(fwer)) {
    for (col in 1:ncol(fwer)) {
      # check for fill vs unfilled
      if (col != 3) {
        points(jitter(as.numeric(rownames(fwer)), amount = 0.01), fwer[,col], col = mypal[col], pch = (col - 1),
               lwd = 2)
      }
      else {
        points(jitter(as.numeric(rownames(fwer)), amount = 0.1), fwer[,col], col = mypal[col], pch = (5),
               lwd = 2)
      }
    }
  }
  abline(h = 0.05, col = "grey")
  # legend('topleft',c('','name'),lty=c(1,NA),pch=c(NA,'X'),bg='white',ncol=2)
  legend(0, 1, legend=legendNames, col = mypal, lty = 1:4, pch = 1:4 + 14, ncol = 2)
  # dev.off()
}

simDf <- read_excel('simulations.xlsx')[c(1, 6:9),]
simDf$Time <-str_pad(simDf$Time, width=5, side="left", pad="0")

mainFunc <- function() {
  intCols <- c(1, 5, 2, 7.5, 0.5)
  B.vec <- c(5, 10, 20, 50) # number of splits
  frac.vec <- c(0.5, 0.75, 0.9, 0.95, 0.99) # selection fraction
  sparsity <- 5
  powerMat <- matrix(NA, nrow = length(frac.vec), ncol = length(intCols))
  fwerMat <- matrix(NA, nrow = length(frac.vec), ncol = length(intCols))
  rownames(powerMat) <- frac.vec
  rownames(fwerMat) <- frac.vec
  colnames(powerMat) <- intCols
  colnames(fwerMat) <- intCols
  for (B in 50) {
    for (frac in frac.vec) {
      fracInd <- match(frac, frac.vec)
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
        retList <- calculatePower(curDat, sparsity)
        # retList <- calculatePowerAdj(curDat, sparsity)
        # power <- retList
        power <- retList[[1]]
        fwer <- retList[[2]]
        subPower <- power[1, 1]
        powerMat[fracInd, row] <- subPower
        fwerMat[fracInd, row] <- fwer[1]
      }
    }
  }
  # powerMat <- powerMat[order(row.names(powerMat)), ]
  # fwerMat <- fwerMat[order(row.names(fwerMat)), ]
  plotVals(powerMat, B, frac, fwerMat)
}

mainFunc()

# plotVals(powerMat)


