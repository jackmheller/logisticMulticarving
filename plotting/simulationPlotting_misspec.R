

library(BBmisc)
library(readxl)
library(stringr)
library(ggplot2)

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
  # legendNames <- c(expression(paste("Carve (", gamma, " = 0.05)")), expression(paste("Split (", gamma, " = 0.05)")), expression(paste("Carve (", gamma, " = 0.30")), 
  #                  expression(paste("Split (", gamma, " = 0.30)")))
  mypal <- colorRampPalette(c("black", "red", "green", "blue") )(4)
  # jpeg(file=paste0("Power with B = ", B, " Selection Frac = ", frac, ".jpeg"), height = 10, width = 10, units = 'in', res = 300)
  # Grouped
  barplot(df, beside = TRUE)
  library(reshape)
  library(plyr)
  df2 <- melt(df)
  fwer2 <- melt(fwer)
  # colnames(df2) <- c('Link', 'Selection Fraction', 'Power')
  df2[,1] <- mapvalues(df2[,1], from=c("cloglog", "Logit", "bbinom03"), to=c("Cloglog", "Logit", "Beta-Binomial"))
  fwer2[,1] <- mapvalues(fwer2[,1], from=c("cloglog", "Logit", "bbinom03"), to=c("Cloglog", "Logit", "Beta-Binomial"))
  # ggplot(df2, aes(x = as.factor(X2), y = value, fill = as.factor(X1))) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 0.05) + geom_point(data = fwer2, aes(x = as.factor(X2), y = value, fill = as.factor(X1)), position = position_dodge(width = 1))
  ggplot(df2, aes(x = as.factor(X2), y = value, fill = as.factor(X1))) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 0.05) + xlab('Selection Fraction') + ylab('Power') + ggtitle('Power of Different Link Functions with B = 50') + 
    theme(plot.title = element_text(hjust = 0.5, size = 22), axis.text=element_text(size=18), axis.title=element_text(size=18), legend.title=element_text(size=18), legend.text=element_text(size=16)) + guides(fill=guide_legend(title="Link")) + 
    geom_point(data = fwer2, aes(x = as.factor(X2), y = value, fill = as.factor(X1)), position = position_dodge(width = 1))
  print('here')
}

simDf <- read_excel('simulations.xlsx')[c(1, 3),]
simDf$Time <-str_pad(simDf$Time, width=5, side="left", pad="0")

mainFunc <- function() {
  intCols <- c('carve5', 'split5', 'carve30', 'split30')
  B.vec <- c(5, 10, 20, 50) # number of splits
  frac.vec <- c(0.5, 0.75, 0.9, 0.95, 0.99) # selection fraction
  powerMat <- matrix(NA, nrow = nrow(simDf), ncol = length(frac.vec))
  fwerMat <- matrix(NA, nrow = nrow(simDf), ncol = length(frac.vec))
  rownames(powerMat) <- simDf$Link
  rownames(fwerMat) <- simDf$Link
  colnames(powerMat) <- frac.vec
  colnames(fwerMat) <- frac.vec
  sparsity <- 5
  for (B in 50) {
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
        retList <- calculatePower(curDat, sparsity)
        # retList <- calculatePowerAdj(curDat, sparsity)
        # power <- retList
        power <- retList[[1]]
        fwer <- retList[[2]]
        subPower <- power[1, 2]
        powerMat[row, match(frac, frac.vec)] <- subPower
        fwerMat[row, match(frac, frac.vec)] <- fwer[2]
      }
    }
    # powerMat <- powerMat[order(row.names(powerMat)), ]
    # fwerMat <- fwerMat[order(row.names(fwerMat)), ]
    plotVals(powerMat, B, frac, fwerMat)
  }
}

mainFunc()

# plotVals(powerMat)


