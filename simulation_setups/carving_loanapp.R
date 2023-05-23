rm(list = ls(all = TRUE))
save <- TRUE

if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("simulation_setups/multi_carve/Binomial_", newdir, sep="")) 
}

require(MASS)
require(glmnet)
require(Matrix)
require(tictoc)
require(hdi)
require(selectiveInference)
require(doSNOW)
require(parallel)
require(doRNG)
require(tmg)
require(git2r)

commit <- revparse_single(revision = "HEAD")
print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))


# adjust depending on folder structure
source("inference/hdi_adjustments.R")
source("inference/carving.R")
source("inference/sample_from_truncated.R")
source("inference/tryCatch-W-E.R")

# toeplitz

# VK 08.04.23 Including code for loanapp (until B.vec)

load(file = "loanappInter.Rdata")

x <- as.matrix(loanapp.inter[, -1])
x <- x[, colSums(is.na(x)) == 0]
lincomb <- caret::findLinearCombos(x)
x <- x[, -lincomb$remove]

# Find the Duplicated Columns
duplicated_columns <- duplicated(t(x))

# for (i in colnames(x)) {
#   ident <- identical(x[,'bd1:mortlat11'],x[,i])
#   if (ident == TRUE) {
#     print(i)
#   }
# }

# Remove the Duplicated Columns
x <- x[,!duplicated_columns]

n <- dim(x)[1]
p <- dim(x)[2]
report.sigma <- FALSE
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
sparsity <- length(sel.index)
level <- 0.05
# table than the output 


B.vec <- c(1) # c(1, (1:5) * 10) # number of splits
frac.vec <- c(0.5) # selection fraction

nsim <- 2
ntasks <- nsim
progress <- function(n, tag) {
  mod <- 16
  if (n %% mod == 0 ) {
    cat(sprintf('tasks completed: %d; tag: %d\n', n, tag))
  }
  if (n %% mod == 0 ) {
    toc()
    tic()
  }
}
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed.v <- sample(1:10000, length(frac.vec))
print(seed.v) # 3588 3052 2252 5257 8307
seed.n <- 0
B <- max(B.vec)

for (frac in frac.vec) {
  seed.n <- seed.n + 1
  set.seed(seed.v[seed.n])
  print(frac)
  print(B)
  opts <- list(progress = progress)
  
  # parallelization
  # choose different number of cores if wished
  cl<- makeSOCKcluster(16)
  
  rseed <- seed.v[seed.n]
  clusterSetRNGStream(cl, iseed = rseed) #make things reproducible
  registerDoSNOW(cl)
  tic()
  # res<-foreach(gu = 1:nsim, .combine = rbind,
  #              .packages = c("MASS", "selectiveInference", "glmnet", "Matrix",
  #                            "hdi", "tmg", "truncnorm", "tictoc") ,.options.snow=opts) %dorng%{
                               # alternative if sequential computation is preferred
                               res<-foreach(gu = 1:nsim, .combine = rbind) %do%{

                               browser()
                               y <- as.matrix(loanapp.inter[,1])
                               c100try <- tryCatch_W_E(carve100(x, y, model.selector = lasso.cvcoef,
                                                                args.model.selector = list(standardize = FALSE, intercept = TRUE, tol.beta = 1e-5, use.lambda.min = TRUE),
                                                                verbose = FALSE, FWER = FALSE, return.selmodels = TRUE,
                                                                family = "binomial"), 0)

                               mcrtry <- tryCatch_W_E(multi.carve(x, y, B = B, fraction = frac, model.selector = lasso.cvcoef, classical.fit = glm.pval.pseudo,
                                                                  parallel = FALSE, ncores = getOption("mc.cores", 2L), gamma = 1, skip.variables = FALSE,
                                                                  args.model.selector = list(standardize = FALSE, intercept = TRUE, tol.beta = 1e-15, use.lambda.min = TRUE),
                                                                  args.classical.fit =list(), verbose = FALSE, FWER = FALSE, split.pval = TRUE, family = "binomial",
                                                                  return.selmodels = TRUE, return.nonaggr = TRUE), 0)
                               
                               # mcrtry <- c100try
                               # mcrtry$value[[1]] <- c100try$value
                               # mcrtry$value[[2]] <- c100try$value
                               # mcrtry$error <- NULL
                               # mcrtry$warning <- NULL

                               out.list <- list()
                               out.list$y <- y
                               if (!is.null(mcrtry$error) || !is.null(c100try$error)) {
                                 # error handling
                                 err <- paste("mcr:", mcrtry$error, "carve100:", c100try$error)
                                 war <- if (is.null(mcrtry$warning) || is.null(c100try$warning)) NA
                                 else c(mcrtry$warning, c100try$warning)
                                 
                                 for (b in B.vec) {
                                   if (b > 1) {
                                     out.list[[as.character(b)]] <- rep(NA, (sparsity + 1) * 16)
                                   } else {
                                     out.list[[as.character(b)]] <- rep(NA, (sparsity + 1) * 6 + 9)
                                   }
                                 }
                                 out.list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
                                 out.list
                               } else {
                                 mcr <- mcrtry$value
                                 pcarve.nofwer <- mcr[[1]]$pvals.nonaggr
                                 psplit.nofwer <- mcr[[2]]$pvals.nonaggr
                                 model.size <- apply(mcr[[1]]$sel.models, 1, sum)
                                 model.size[model.size == 0]  <- 1 # if no variable is selected, p-values are 1
                                 # 4/27/23 JMH add betas in
                                 betas.carve100 <- c100try$value$coefs
                                 betas.carve100.vars <- c100try$value$carCoefs
                                 betasList <- list()
                                 betasList[[paste0('carve_coefs')]] <- betas.carve100
                                 betasList[[paste0('carve_vars')]] <- betas.carve100.vars
                                 
                                 out.list$betas <- betasList
                                 
                                 # ommit the clipping to calculate adjusted power
                                 # 20/3/23 JMH/VK change from below to cap as in Meinshausen 2.1
                                 if (B > 1) {
                                   pcarve.fwer <- pmin(pcarve.nofwer * model.size, 1)
                                   psplit.fwer <- pmin(psplit.nofwer * model.size, 1)
                                 }
                                 else {
                                   # 20/3/23 JMH/VK applying the single-split method from Meinshausen (2.1)
                                   pcarve.fwer <- pcarve.nofwer * model.size
                                   psplit.fwer <- psplit.nofwer * model.size
                                 }
                                 # pcarve.fwer <- pcarve.nofwer * model.size
                                 # psplit.fwer <- psplit.nofwer * model.size
                                 c100 <- c100try$value
                                 pc100.nofwer <- c100$pval.corr
                                 model.size100 <- sum(c100$sel.models)
                                 model.size100[model.size100 == 0]  <- 1
                                 # ommit the clipping to calculate adjusted power
                                 # 15/2/23 JMH/VK change from below to cap as in Meinshausen 2.1
                                 if (B > 1) {
                                   pc100.fwer <- pmin(pc100.nofwer * model.size100, 1)
                                 }
                                 else {
                                   pc100.fwer <- pc100.nofwer * model.size100
                                 }
                                 # pc100.fwer <- pc100.nofwer * model.size100
                                 for (B in B.vec) {
                                   if (B > 1) {
                                     use <- 1:B #JMH/ VK 20/3/23 cutoff from FALSE to TRUE
                                     pvals.aggregated <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                                                         round(seq(ceiling(0.05 * B)/B, 1 - 1/B, by = 1/B), 2), cutoff = TRUE)
                                     pvals.aggregated2 <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                                                          round(seq(ceiling(0.3 * B)/B, 1 - 1/B, by = 1/B), 2), cutoff = TRUE)
                                     pvals.aggregated3 <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                                                          round(ceiling(0.05 * B)/B, 2), cutoff = TRUE)
                                     pvals.aggregated4 <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                                                          round(ceiling(0.3 * B)/B, 2), cutoff = TRUE)
                                   } else {
                                     pvals.aggregated <- list(pcarve.nofwer[1,], pcarve.fwer[1,], psplit.nofwer[1,], psplit.fwer[1,])
                                   }
                                   
                                   run.res <- vector(length = 0) # store important quantities
                                   np <- length(pvals.aggregated)
                                   for (i in 1:np) {
                                     true.pv <- pvals.aggregated[[i]][ind] # p-values of active variables
                                     bad.pv <- min(pvals.aggregated[[i]][-ind]) # lowest p-value of inactive variables to check FWER
                                     run.res <- c(run.res, true.pv, bad.pv)
                                   }
                                   if (B > 1) {
                                     for (i in 1:np) {
                                       true.pv <- pvals.aggregated2[[i]][ind]
                                       bad.pv <- min(pvals.aggregated2[[i]][-ind])
                                       run.res <- c(run.res, true.pv, bad.pv)
                                     }
                                     for (i in 1:np) {
                                       true.pv <- pvals.aggregated3[[i]][ind]
                                       bad.pv <- min(pvals.aggregated3[[i]][-ind])
                                       run.res <- c(run.res, true.pv, bad.pv)
                                     }
                                     for (i in 1:np) {
                                       true.pv <- pvals.aggregated4[[i]][ind]
                                       bad.pv <- min(pvals.aggregated4[[i]][-ind])
                                       run.res <- c(run.res, true.pv, bad.pv)
                                     }
                                   }
                                   # 20/3/23 JMH/VK adapt and add R, TS, V for all values
                                   # 23/3/23 JMH/VK add row subset to B rows to make calculations correct
                                   R <- length(which(mcr[[1]]$sel.models[1:B, ])) / B
                                   if (B == 1) {
                                     TS <- sum(which(mcr[[1]]$sel.models[1:B, ], arr.ind = TRUE) %in% ind) / B
                                   }
                                   else {
                                     TS <- sum(which(mcr[[1]]$sel.models[1:B, ], arr.ind = TRUE)[,2] %in% ind) / B
                                   }
                                   V <- R - TS
                                   run.res <- c(run.res, R, TS, V)
                                   if (B == 1) {
                                     # analyse first split specially for B = 1 and analyse carve100
                                     R <- length(which(mcr[[1]]$sel.models[1, ])) # number of variables selected in first split
                                     # Does TS need to be a float and not binary?
                                     TS <- sum(ind %in% which(mcr[[1]]$sel.models[1, ])) # number of active variables selected
                                     V <- R - TS # number of inactive variables selected
                                     carve.err <- sum(pvals.aggregated[[1]][-ind] < level) # number of false rejection from single-carving #17/02/23 VK, setting significance level only once
                                     split.err <- sum(pvals.aggregated[[3]][-ind] < level) # number of false rejection from single-splitting #17/02/23 VK, setting significance level only once
                                     carve100.err <- sum(pc100.nofwer[-ind] < 0.05)
                                     # run.res <- c(run.res, R, V, TS) #20/3/23 JMH, VK commented out 
                                     true.pv <- pc100.nofwer[ind]
                                     bad.pv <- min(pc100.nofwer[-ind])
                                     R100 <- length(which(c100$sel.models))
                                     TS100 <- sum(ind %in% which(c100$sel.models))
                                     V100 <- R100 - TS100
                                     run.res <- c(run.res, true.pv, bad.pv)
                                     true.pv <- pc100.fwer[ind]
                                     bad.pv <- min(pc100.fwer[-ind])
                                     run.res <- c(run.res, true.pv, bad.pv, R100, V100, TS100,
                                                  carve.err, split.err, carve100.err)
                                   }
                                   out.list[[as.character(B)]] <- run.res
                                   # 22/3/23 JMH add pvals aggregated
                                   out.list[[paste0(as.character(B), '_carveNoFWER')]] <- pcarve.nofwer
                                   out.list[[paste0(as.character(B), '_carveFWER')]] <- pcarve.fwer
                                   out.list[[paste0(as.character(B), '_splitNoFWER')]] <- psplit.nofwer
                                   out.list[[paste0(as.character(B), '_splitFWER')]] <- psplit.fwer
                                 }
                                 err <- if (is.null(mcrtry$error) && is.null(c100try$error)) NA
                                 else c(mcrtry$error, c100try$error) # should not happen due to earlier check
                                 war <- if (is.null(mcrtry$warning) && is.null(c100try$warning)) NA
                                 else c(mcrtry$warning, c100try$warning)
                                 out.list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
                                 out.list
                                 # end of simulation run
                               }
                               # end of simulation for given fraction
                             }
  toc()
  stopCluster(cl)
  # analyse results for given fraction
  expmatr <- matrix(unlist(res[, "exception"]), nrow = dim(res)[1], ncol = 2, byrow = TRUE)
  print(sum(is.na(expmatr[, 1])))
  succ = which(is.na(expmatr[, 1]))
  print("succesful runs")
  listPvals <- list()

  # 4/27/23 JMH add betas to return data
  betasList <- list()
  counter <- TRUE
  for (i in res[, "betas"]) {
    if (counter) {
    t <- matrix(unlist(i), nrow = 2, byrow = TRUE)
    colnames(t) <- names(i$carve_vars)
    
    betasList[[counter]] <- t
    counter <- FALSE
    }
  }
  
  all.y <- matrix(unlist(res[,"y"]), nrow = dim(res), byrow = TRUE)
  sd <- attr(res, "rng")
  for (B in B.vec) {
    if (B == 1) {
      names1 <- c("carve","carvefw", "split", "splitfw")
      names2 <- c("carve100","carve100fw")
      subres <- matrix(unlist(res[,as.character(B)]), nrow = dim(res)[1],
                       ncol = 6 * (sparsity + 1) + 9, byrow = TRUE)
      if (any(!is.na(subres[-succ, ]))) print("not as it should be") # sanity check
      subres <- subres[succ,]
      colnames(subres) <- c(rep(names1, each = (sparsity + 1)), "R", "V", "R-V",
                            rep(names2, each = (sparsity + 1)), "R100", "V100",
                            "R-V100", "carve.err", "split.err", "carve100.err")
      names <- c(names1, names2)
      # 22/3/23 JMH add pvals aggregated
      listPvals[[paste0(as.character(B), '_carveNoFWER')]] <- res[, paste0(as.character(B), '_carveNoFWER')]
      listPvals[[paste0(as.character(B), '_carveFWER')]] <- res[, paste0(as.character(B), '_carveFWER')]
      listPvals[[paste0(as.character(B), '_splitNoFWER')]] <- res[, paste0(as.character(B), '_splitNoFWER')]
      listPvals[[paste0(as.character(B), '_splitFWER')]] <- res[, paste0(as.character(B), '_splitFWER')]
    } else {
      names<-c("carve5", "carvefw5", "split5", "splitfw5", "carve30",
               "carvefw30", "split30", "splitfw30", "carvefix5",
               "carvefwfix5", "splitfix5", "splitfwfix5", "carvefix30",
               "carvefwfix30", "splitfix30", "splitfwfix30")
      subres <- matrix(unlist(res[,as.character(B)]), nrow = dim(res)[1],
                       ncol = length(names) * (sparsity + 1) + 3, byrow = TRUE) #17/02/23 VK, replacing 16, +3 for R, V, R-V
      if (any(!is.na(subres[-succ, ]))) print("not as it should be")
      subres <- subres[succ,]
      # 20/3/23 JMH/VK add R, V, R-V cols
      colnames(subres) <- c(rep(names, each = (sparsity + 1)), "R", "V", "R-V")
      # 22/3/23 JMH add pvals aggregated
      listPvals[[paste0(as.character(B), '_carveNoFWER')]] <- res[, paste0(as.character(B), '_carveNoFWER')]
      listPvals[[paste0(as.character(B), '_carveFWER')]] <- res[, paste0(as.character(B), '_carveFWER')]
      listPvals[[paste0(as.character(B), '_splitNoFWER')]] <- res[, paste0(as.character(B), '_splitNoFWER')]
      listPvals[[paste0(as.character(B), '_splitFWER')]] <- res[, paste0(as.character(B), '_splitFWER')]
    }
    subres <- as.data.frame(subres)
    # 20/3/23 VK add selection index, 22/3/23 JMH/VK add pvals.aggregated
    simulation <- list("results" = subres, "exceptions" = expmatr, "y" = all.y, "B" = B, "split" = frac,
                       "nsim" = nsim, "seed" = rseed, "All used B" = B.vec, "sd" = sd, "commit" = commit, "sparsity"=sparsity, "sel.index"=sel.index,
                       "pvals.aggregated" = listPvals, "betas" = betasList)
    print(paste("results using fraction ", frac, " and B=", B, sep = ""))
    # if (B == 1) {
    #   print(mean(subres$`R-V` == sparsity)) # probability of screening
    #   good <- which(subres$`R-V` == sparsity)
    #   print(apply(subres[, c("R", "V", "R-V")], 2, mean)) # totally, active, inactive selected
    #   # probability of screening using all data for selection
    #   print(mean(subres$`R-V100` == sparsity)) 
    #   good100 <- which(subres$`R-V100` == sparsity)
    #   # totally, active and inacitve selected using all data for selection
    #   print(apply(subres[, c("R100", "V100", "R-V100")], 2, mean))
    #   print(c(c(sum(subres$carve.err), sum(subres$split.err)) / sum(subres$V),
    #           sum(subres$carve100.err) / sum(subres$V100))) # Rejection amongst falsely selected
    #   # Rejection amongst falsely selected conditioned on screening, should be below 0.05
    #   print(c(c(sum(subres$carve.err[good]), sum(subres$split.err[good])) / sum(subres$V[good]),
    #           sum(subres$carve100.err[good100]) / sum(subres$V100[good100]))) 
    # } 
    # allrej <- matrix(NA, nrow = 1,ncol = length(names))
    # colnames(allrej) <- names
    # for (name in names) {
    #   nameind <- which(colnames(subres) == name)
    #   mat <- subres[, nameind]
    #   rej <- quantile(mat[, sparsity + 1], level, na.rm = TRUE) #17/02/23 VK, setting significance level only once
    #   rejmat <- mat[, 1:sparsity] < rej
    #   allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
    #   
    # }
    # print("Adjusted")
    # print(allrej) # adjusted power
    # print("Unadjusted")
    # fwer <- numeric(length(names))
    # names(fwer) <- names
    # for (name in names) {
    #   nameind <- which(colnames(subres) == name)
    #   mat <- subres[, nameind]
    #   fwer[as.character(name)] <- mean(mat[,sparsity + 1] < level, na.rm = TRUE)#20/3/23 0.05 to level
    #   rejmat <- mat[, 1:sparsity] < level #VK 20/3/23 level 
    #   allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
    # }
    # print(fwer) # FWER
    # print(allrej) # power
    resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"),
                      " split=", frac, " B=", B, " seed=", rseed)
    # adjust depending on folder structure
    if (save) save(simulation,file = paste("simulation_setups/multi_carve/Binomial_", newdir, "/", resname, ".RData", sep = ""))
    
  }
}

print("Finale")
