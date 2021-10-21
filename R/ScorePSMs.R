library(statmod)
library(mzR)

###Input:
#Spectra and Headers: the output from my de-isotope function.
#data.F: the database. This is either the forward or reverse database
#data.R: due to poor database management on my part, this should be left as null
#n_cores: should be left as 1
#match.charge: should be left as TRUE
#append.file: a character string that specifies the parameters to use. This should be kept as "_Train.Q9.lowpH" when searching the high pH data
#scale: specifies how to normalize intensities. Defaults to Quant90 (normalization by the 90th percentile). This should be Quant90 for now.
#path.spline,...,Noise.MA.path: the paths to the hyperparameters directories
#ppm.precursor: Delta mass (in ppm) of the precursor
#min.peaks: The minimum number of peaks a spectrum must contain for it to be scored
#mu.ppm: a hyperparameter that should be left as 0
#delta.ppm: Delta mass (in ppm) of matched b/y MS2 fragments
#delta.ppm.prec: Delta mass (in ppm) of the precursor in the MS2 spectrum. Should be set to 20

###Output: A list of length length #MS2_spectra
#$MS2.indices: The scan number
#$Results: The results of the search, which are contained in $Results$FM. If NA, it means that no peptides were scored.
#          Otherwise, it a #PeptidesSearched x 10 matrix. The score for each peptide is the sum of the first four columns

#' Search MS/MS spectra with MSeQUiP
#' 
#' Computes MSeQUiP Bayes factors for each candidate peptide
#' 
#' @param Spectra A list of MS and MS/MS spectra. Can be obtained using mzR, or \code{Deisotope}.
#' @param Headers A data frame of headers. Can be obtained using mzR, or \code{Deisotope}.
#' @param data.F Spectral library. This can be obtained from Prosit.
#' @param clean.database Should the spectral library be cleaned. Defaults to \code{F}.
#' @param match.charge Should inferred fragment charges obtained from \code{Deisotope} be used. Defaults to \code{T}.
#' @param append.file Name of training parameters to use. Detaults to \code{"_Train.Q9.lowpH"}.
#' @param append.noise Name of training parameters to use. Detaults to \code{append.file}.
#' @param scale Type of intensity normalization to use. Can be one of \code{"Quant90"} (normalize by the 0.9 quantile), \code{"base-peak"} (normalize by the base peak), or \code{"none"} (no normalization). Defaults to \code{"Quant90"}.
#' @param path.spline Path to parameters to define lambda(m), the probability a noise peak was generated around a predicted m/z m
#' @param rt.path Path to retention time parameters.
#' @param byInt.path Path to signal b/y intensity parameters.
#' @param byMA.path Path to signal mass accuracy parameters.
#' @param byObs.path Path to parameters determine the probability a b/y ion is observed.
#' @param NoiseInt.path Path to parameters describing the noise intensity.
#' @param Noise.MA.path Path to parameters describing the noise mass accuracy.
#' @param ppm.precursor Precursor mass accuracy, on the ppm scale. Defaults to 10.
#' @param min.peaks Minimum number of peaks in the MS/MS spectrum. Defaults to 20.
#' @param delta.ppm MS/MS fragment mass accuracy, on the ppm scale. Defaults to 20.
#' @param delta.ppm.prec Maximum difference between the precursors observed in MS and MS/MS spectra.
#' @param scale.obs The intensity scale, if available. The default of \code{NA} will almost always be appropriate.
#' @param include.scale If \code{scale.obs} is available, should it be used. Defaults to \code{T}.
#' @param include.pep.IDs Should peptide IDs be included in the results. Defaults to \code{T}.
#' 
#' @return A list. \item{MS2.indices}{The indices of \code{Spectra} that are MS/MS spectra.} \item{Results}{A list. List element \code{FM} contains the results, and is a list of matrices giving the MSeQUiP scores and other statistics for each scored peptide. The sum of the first 4 columns gives the log Bayes factor for the peptide. The element \code{RM} is NA.}
#' @export
PerformSearch <- function(Spectra, Headers=NULL, data.F, clean.database=F, match.charge=T, append.file="_Train.Q9.lowpH",
                          scale=c("Quant90", "base-peak", "none"), path.spline="<path to>/NoiseObserved",
                          rt.path="<path to>/RetentionTime", byInt.path="<path to>/BYInt",
                          byMA.path="<path to>/BYMassAccuracy", byObs.path="<path to>/BYObserved",
                          NoiseInt.path="<path to>/NoiseInt", Noise.MA.path="<path to>/NoiseMassAccuracy", Int0.path=NULL,
                          ppm.precursor=10, min.peaks=20, delta.ppm = 20, delta.ppm.prec = 20, scale.obs=NA, include.scale=T, append.noise=NULL, include.pep.IDs=T) {
  
  n_cores <- 1
  data.R <- NULL
  mu.ppm <- 0
  if (is.null(append.noise)) {append.noise <- append.file}
  if (is.character(Spectra)) {
    tmp <- Deisotope(mzml.file.path = Spectra, n_cores = 3, ppm = 20, ppm.filter = 30, write.mzml = F, max.addition = 2, remove.precursor = F)
    Spectra <- tmp$Spectra
    Headers <- tmp$Headers
    rm(tmp)
  }
  if (!is.na(scale.obs[1])) {scale <- "none"}
  if (clean.database) {
    cat("Cleaning database...")
    data.F <- lapply(data.F,function(x){x$mat<-Clean.prosit(Prosit = x$mat, ppm.filter = 40, add = T); return(x)})
    if (!is.null(data.R)) {data.R <- lapply(data.R,function(x){x$mat<-Clean.prosit(Prosit = x$mat, ppm.filter = 40, add = T); return(x)})}
    cat("done!\n")
  }
  charge.bin <- list(2,3,4:10)
  scale <- match.arg(scale, c("Quant90", "base-peak", "none"))
  z.F <- sapply(data.F, function(x){x$z})
  mz.F <- sapply(data.F, function(x){x$mz})
  tmp <- order(mz.F); z.F <- z.F[tmp]; mz.F <- mz.F[tmp]; data.F <- data.F[tmp]; rm(tmp)
  Ind.F.Theory <- sapply(1:10,function(zz){which(z.F==zz)})
  if (include.pep.IDs) {Pep.IDs.F <- paste(sapply(data.F, function(x){x$Pep}),z.F,sep=".")}
  if (!is.null(data.R)) {
    z.R <- sapply(data.R, function(x){x$z})
    mz.R <- sapply(data.R, function(x){x$mz})
    tmp <- order(mz.R); z.R <- z.R[tmp]; mz.R <- mz.R[tmp]; data.R <- data.R[tmp]; rm(tmp)
    Ind.R.Theory <- sapply(1:10,function(zz){which(z.R==zz)})
    if (include.pep.IDs) {Pep.IDs.R <- paste(sapply(data.R, function(x){x$Pep}),z.R,sep=".")}
  }
  
  #Get spline model#
  tmp.spline <- Read.NoiseObs(path.dir = path.spline, append.file = ifelse(length(path.spline)>1,NA,append.noise))
  Spline.Fit <- tmp.spline$Spline.Fit; Spline.cutoff <- tmp.spline$Spline.cutoff
  
  #Int 0 model#
  if (!is.null(Int0.path)) {
    Int0.model <- Read.AbsentIons(dir.path = Int0.path, append.file=ifelse(length(Int0.path)>1,NULL,append.noise), parametrization="binom")
  } else {
    Int0.model <- NULL
  }
  
  #RT#
  param.RT <- Read.RT(path.dir = rt.path, ifelse(length(rt.path)>1,NA,append.file))
  
  #BY Int#
  param.BYInt <- Read.BYint(path.dir = byInt.path, append.file = ifelse(length(byInt.path)>1,NA,append.file))
  
  #BY Observed#
  param.BYObs <- Read.BYobs(path.dir = byObs.path, append.file = ifelse(length(byObs.path)>1,NA,append.file))
  
  #MA Observed#
  param.MA <- Read.BYMA(path.dir = byMA.path, append.file = ifelse(length(byMA.path)>1,NA,append.file))
  
  #MA Noise#
  param.MA.noise <- Read.NoiseMA(path.dir = Noise.MA.path, append.file = ifelse(length(Noise.MA.path)>1,NA,append.noise))
  
  #Noise Int#
  param.Noise.Int <- Read.NoiseInt(path.dir = NoiseInt.path, append.file = ifelse(length(NoiseInt.path)>1,NA,append.noise), include.k = T)
  
  #Score PSMs#
  MS2.indices <- which(Headers$msLevel==2)
  n.peptides <- length(MS2.indices)
  cat(paste0("Scoring ", n.peptides, " MS2 spectra..."))
  if (n_cores <= 1) {
    Index.mat <- matrix(NA, nrow=length(MS2.indices), ncol=2)
    for (zz in 1:length(Ind.F.Theory)) {
      Ind.obs.zz <- which(Headers$precursorCharge[MS2.indices] == zz)
      if (length(Ind.F.Theory[[zz]]) > 0 && length(Ind.obs.zz) > 0) {Index.mat[Ind.obs.zz,] <- Match.Prec(mz.prec = Headers$precursorMZ[MS2.indices][Ind.obs.zz], mz.theory = mz.F[Ind.F.Theory[[zz]]], delta.ppm = ppm.precursor, mu.ppm = mu.ppm)}
    }

    out <- lapply(1:length(MS2.indices), function(i){
      if (round(i/1e4)==i/1e4) {cat(paste0(i, "/", n.peptides, "..."))}
      scan <- MS2.indices[i]
      if (length(Spectra[[scan]])<min.peaks*3 || is.na(Index.mat[i,1])) {return(list(FM=NA,RM=NA))}
      ind.forward <- Ind.F.Theory[[Headers$precursorCharge[scan]]][Index.mat[i,1]:Index.mat[i,2]]
      #ind.forward <- which(abs(Headers$precursorMZ[scan]/mz.F-1) <= ppm.precursor/10^6 & z.F==Headers$precursorCharge[scan])
      Forward.match <- FindMatch(S.obs = Spectra[[scan]], rt.obs = Headers$retentionTime[scan]/60, S.theo.list = lapply(data.F[ind.forward],function(x){x$mat}), rt.theo.vec = sapply(data.F[ind.forward],function(x){x$irt}), prec.mz = Headers$precursorMZ[scan], z = Headers$precursorCharge[scan], min.match = 0, delta.ppm = delta.ppm, mu.ppm = mu.ppm, delta.ppm.prec = delta.ppm.prec, return.all=T, match.charge = match.charge, Spline.Fit = Spline.Fit, Spline.cutoff = Spline.cutoff, Int0.model=Int0.model, scale = scale, param.RT = param.RT, param.BYInt = param.BYInt, param.BYObs = param.BYObs, param.MA = param.MA, param.Noise.Int = param.Noise.Int, param.MA.noise = param.MA.noise, scale.obs = ifelse(is.na(scale.obs[1]),NA,scale.obs[scan]), include.scale = include.scale)
      if (include.pep.IDs) {rownames(Forward.match) <- Pep.IDs.F[ind.forward]}
      
      if (!is.null(data.R)) {
        ind.decoy <- which(abs(Headers$precursorMZ[scan]/mz.R-1) <= ppm.precursor/10^6 & z.R==Headers$precursorCharge[scan])
        Reverse.match <- FindMatch(S.obs = Spectra[[scan]], rt.obs = Headers$retentionTime[scan]/60, S.theo.list = lapply(data.R[ind.decoy],function(x){x$mat}), rt.theo.vec = sapply(data.R[ind.decoy],function(x){x$irt}), prec.mz = Headers$precursorMZ[scan], z = Headers$precursorCharge[scan], min.match = 0, delta.ppm = delta.ppm, mu.ppm = mu.ppm, delta.ppm.prec = delta.ppm.prec, return.all=T, match.charge = T, return.BF = F, return.BF.full = T, null.noise.mean = T, Spline.Fit = Spline.Fit, Spline.cutoff = Spline.cutoff)
      } else {
        Reverse.match <- NA
      }
      
      return(list(FM=Forward.match,RM=Reverse.match))
    })
    cat("done!\n")
    return(list(Results=out,MS2.indices=MS2.indices))
  }
  
  ##Create cluster##
  cl <- makeCluster(ifelse(is.null(n_cores),max(detectCores()-1,1),max(detectCores()-1,n_cores)))
  clusterEvalQ(cl = cl, expr = {
    source("../R/Scoring/ScorePSMs.R")
    library(splines)
  })
  clusterExport(cl = cl, varlist = c("data.F", "data.R", "z.F", "z.R", "mz.F", "mz.R", "Spline.Fit", "Spline.cutoff",
                                     "ppm.precursor", "min.peaks", "mu.ppm", "delta.ppm", "delta.ppm.prec"), envir = environment())
  out <- parLapply(cl = cl, X = lapply(MS2.indices,function(x){list(s=Spectra[[x]],h=Headers[x,])}), function(x){
    if (length(x$s) < 3*min.peaks) {return(list(FM=NA,RM=NA))}
    ind.forward <- which(abs(x$h$precursorMZ/mz.F-1) <= ppm.precursor/10^6 & z.F==x$h$precursorCharge)
    if (length(ind.forward) == 0) {return(list(FM=NA,RM=NA))}
    ind.decoy <- which(abs(x$h$precursorMZ/mz.R-1) <= ppm.precursor/10^6 & z.R==x$h$precursorCharge)
    Forward.match <- FindMatch(S.obs = x$s, rt.obs = x$h$retentionTime/60, S.theo.list = lapply(data.F[ind.forward],function(m){m$mat}), rt.theo.vec = sapply(data.F[ind.forward],function(m){m$irt}), prec.mz = x$h$precursorMZ, z = x$h$precursorCharge, Spline.Fit = Spline.Fit, Spline.cutoff = Spline.cutoff, min.match = 0, delta.ppm = delta.ppm, mu.ppm = mu.ppm, delta.ppm.prec = delta.ppm.prec, return.all=T, match.charge = T, return.BF = F, return.BF.full = T)
    Reverse.match <- FindMatch(S.obs = x$s, rt.obs = x$h$retentionTime/60, S.theo.list = lapply(data.R[ind.decoy],function(m){m$mat}), rt.theo.vec = sapply(data.R[ind.decoy],function(m){m$irt}), prec.mz = x$h$precursorMZ, z = x$h$precursorCharge, Spline.Fit = Spline.Fit, Spline.cutoff = Spline.cutoff, min.match = 0, delta.ppm = delta.ppm, mu.ppm = mu.ppm, delta.ppm.prec = delta.ppm.prec, return.all=T, match.charge = T, return.BF = F, return.BF.full = T)
    return(list(FM=Forward.match,RM=Reverse.match))
  })
  stopCluster(cl = cl)
  return(list(Results=out,MS2.indices=MS2.indices))
}

FindMatch <- function(S.obs, rt.obs, S.theo.list, rt.theo.vec, prec.mz, z, scale=c("Quant90", "base-peak", "none"), param.RT, param.BYInt, param.BYObs, param.MA, param.MA.noise, param.Noise.Int, Spline.Fit, Spline.cutoff, Int0.model=NULL, min.match=3, delta.ppm=20, mu.ppm=0, delta.ppm.prec=20, return.all=F, match.charge=T, scale.obs=NA, include.scale=F, charge.bins=list(2,3,4:10)) {
  scale <- match.arg(scale, c("Quant90", "base-peak", "none"))
  
  #Remove precursor#
  ind.prec <- which(abs(S.obs[,1]/prec.mz-1) <= delta.ppm.prec/10^6)
  if (length(ind.prec) > 0) {S.obs <- S.obs[-ind.prec,]}
  
  #Re-scale intensity#
  tmp <- Find.Range(mz = prec.mz, z = z)
  S.obs <- S.obs[S.obs[,1]>=tmp[1] & S.obs[,1]<=tmp[2],]
  if (scale == "base-peak") {scale.obs <- ifelse(is.na(scale.obs),max(S.obs[,2]),scale.obs); S.obs[,2] <- S.obs[,2]/scale.obs}
  if (scale == "Quant90") {scale.obs <- ifelse(is.na(scale.obs),quantile(S.obs[,2],0.9),scale.obs); S.obs[,2] <- S.obs[,2]/scale.obs}
  if (is.null(scale.obs) || is.na(scale.obs)) {scale.obs <- 1}
  S.obs <- S.obs[S.obs[,2]>=1e-8,]
  
  #Charge bin and precursor mass#
  c <- which(sapply(charge.bins,function(x){z%in%x}))
  prec.mass <- prec.mz*z - 1.007825032*z
  spline.obs <- Spline.Fit[[c]]
  if (c > 1) {spline.obs <- spline.obs[[ifelse(prec.mass<Spline.cutoff[c],1,2)]]}
  if (!is.null(Int0.model)) {Int0.model <- Int0.model[[c]]}
  
  #Score PSMs#
  out <- sapply(1:length(S.theo.list), function(i){
    Score.PSM(S.obs = S.obs, rt.obs = rt.obs, z = z, S.theo = S.theo.list[[i]], rt.theo = rt.theo.vec[i],
              min.match = min.match, delta.ppm = delta.ppm, mu.ppm = mu.ppm, match.charge = match.charge,
              return.all = return.all, spline.fit = spline.obs, Int0.model=Int0.model, Omega = param.BYObs[[c]]$Omega,
              scale.obs = scale.obs, include.scale=include.scale, mu.AD = param.BYObs[[c]]$ad,
              mu.obs = as.numeric(c(param.BYInt[[c]]$mu[1],param.BYInt[[c]]$beta[1])),
              se.obs = as.numeric(c(param.BYInt[[c]]$mu[2],param.BYInt[[c]]$beta[2])),
              mu.noise = param.BYInt[[c]]$mu.noise[1], se.noise = param.BYInt[[c]]$mu.noise[2],
              d.var = param.BYInt[[c]]$var[1], phi.var = param.BYInt[[c]]$var[2], poly.noise = param.Noise.Int[[c]]$beta,
              coef.noise = list(norm2=param.Noise.Int[[c]]$norm2,alpha=param.Noise.Int[[c]]$alpha),
              min.noise = param.Noise.Int[[c]]$bound[1], max.noise = param.Noise.Int[[c]]$bound[2],
              rho = param.BYInt[[c]]$mu[3], k.noise = param.Noise.Int[[c]]$k, beta.rt = param.RT$beta,
              v.rt = param.RT$sigma^2, lambda.rt = param.RT$lambda, sigma.MA = param.MA$sigma, beta.MA = param.MA$beta,
              sigma.MA.noise = param.MA.noise$sigma, beta.MA.noise = param.MA.noise$beta, bound.mz.MA.noise = sort(param.MA.noise$bound), bound.int.MA = sort(param.MA$bound))
  })
  if (return.all) {return(t(out))}
  return(colSums(out))
}

#####Score PSM#####
#S.obs: #peaks x 2 matrix
  #1st column: m/z
  #2nd column: base-peak normalized intensity
Score.PSM <- function(S.obs, rt.obs, z, S.theo, rt.theo, min.match = 3, delta.ppm = 20, mu.ppm = 0, match.charge = T,
                      return.all = F, spline.fit, Int0.model=NULL, Omega, mu.AD, mu.obs, se.obs, mu.noise, se.noise, d.var, phi.var,
                      poly.noise, coef.noise, min.noise, max.noise, rho, scale.obs, include.scale=F, k.noise, beta.rt, v.rt, lambda.rt, sigma.MA,
                      beta.MA, sigma.MA.noise, beta.MA.noise, bound.mz.MA.noise, bound.int.MA) {
  
  #Find matches#
  if (match.charge) {
    D.mass <- Delta.Mass(mass.find = as.vector(S.theo[,1]), masses = S.obs[,1], sorted = T, charge.find = as.vector(S.theo[,3]), charge.masses = S.obs[,3])
  } else {
    D.mass <- Delta.Mass(mass.find = as.vector(S.theo[,1]), masses = S.obs[,1], sorted = T)
  }
  S.match <- cbind(D.mass$delta, S.obs[D.mass$ind,2], D.mass$ind)
  y.prosit <- S.theo[,2]
  x.prosit <- as.numeric(S.match[,2])
  tmp.ind <- abs(S.match[,1])>delta.ppm; if (sum(tmp.ind)>0){x.prosit[tmp.ind] <- rep(0,sum(tmp.ind))}
  ind.theo.tmp <- S.theo[,2]>1e-8
  score.Int0 <- !is.null(Int0.model) & sum(!ind.theo.tmp) > 0
  if (score.Int0) {
    Int0.data <- GetInt0Info(S.match = S.match[!ind.theo.tmp,], S.theo = S.theo[!ind.theo.tmp,], ppm.ms2 = delta.ppm)
    n.Int0 <- length(Int0.data$int); n.Int0.obs <- sum(Int0.data$int >= 1e-8)
  } else {
    n.Int0 <- 0; n.Int0.obs <- 0
  }
  S.match <- S.match[ind.theo.tmp,]
  S.theo <- S.theo[ind.theo.tmp,]
  ind.obs <- abs(S.match[,1]) <= delta.ppm
  n.obs <- sum(ind.obs); n.total <- length(ind.obs)
  if (n.obs < min.match && !return.all) {return(-Inf)}
  if (n.obs < length(ind.obs)) {S.match[!ind.obs,2] <- 0}
  
  #Get noise#
  if (n.obs > 0) {
    Noise <- S.obs[-S.match[ind.obs,3],]
  } else {
    Noise <- S.obs
  }
  tmp.v <- var(log(Noise[,2]))/length(Noise[,2])
  mean.noise <- (mean(log(Noise[,2]))/tmp.v + mu.noise/se.noise^2)/(1/tmp.v + 1/se.noise^2)
  tmp.v <- var(log(S.obs[,2]))/length(S.obs[,2])
  mean.noise.2 <- (mean(log(S.obs[,2]))/tmp.v + mu.noise/se.noise^2)/(1/tmp.v + 1/se.noise^2)
  
  #RT score#
  like.rt <- RT.score(rt = rt.obs, prosit.rt = rt.theo, M = 90, lambda = lambda.rt, v = v.rt, beta = beta.rt)
  
  #Mass accuracy#
  like.ma <- 0
  #if (n.obs > 0 && !is.na(sigma.MA.noise[1])) {like.ma <- Mass.Accuracy.match(delta = S.match[ind.obs,1], mu = 0, mz.theo = NULL, log.int = log(S.match[ind.obs,2])+ifelse(include.scale&!is.na(scale.obs),log(scale.obs),0), sigma = sigma.MA, beta.int = beta.MA, bound.int.MA = bound.int.MA) - Mass.Accuracy.noise(delta = S.match[ind.obs,1], mz.theo = S.theo[ind.obs,1], ppm.max = delta.ppm, theta = beta.MA.noise, sigma = sigma.MA.noise, bound = bound.mz.MA.noise)}
  if (n.obs > 0 && is.na(sigma.MA.noise[1])) {like.ma <- Mass.Accuracy.score(delta = S.match[ind.obs,1], log.int = log(S.match[ind.obs,2])+ifelse(include.scale&!is.na(scale.obs),log(scale.obs),0), mz.theo = S.theo[ind.obs,1], beta.obs = beta.MA, beta.noise = beta.MA.noise, sigma.obs = sigma.MA, bound.log.int = bound.int.MA, bound.log.mz = bound.mz.MA.noise, ppm.max = delta.ppm)}
  
  #pr(b/y observed | T)#
  like.obs <- by.score.observed(I.obs = S.match[,2], I.pred = S.theo[,2], mu = mu.AD, Omega = Omega, min.match = 0)
  if (!is.list(like.obs)) {like.obs <- list(log.laplace=like.obs)}
  if (n.obs > 0) {like.obs$log.laplace <- like.obs$log.laplace - Prob.random.match(mz = S.theo[ind.obs,1], n.total = nrow(S.obs), spline.fit = spline.fit, min.mz = min(S.obs[,1]), max.mz = max(S.obs[,1]), ppm.ms2 = delta.ppm, dens = spline.fit$density, breaks = spline.fit$Breaks)}
  
  #pr(b/y Int0 observed | T)#
  if (score.Int0) {
    like.Obs.Int0 <- Int0.Obs.score(mz = Int0.data$mz, int = Int0.data$int, param.beta = Int0.model, n.total = nrow(S.obs)-n.obs, spline.fit = spline.fit, min.mz = min(S.obs[,1]), max.mz = max(S.obs[,1]), ppm.ms2 = delta.ppm, dens = spline.fit$density, breaks = spline.fit$Breaks)
  } else {
    like.Obs.Int0 <- 0
  }
  
  #pr(b/y abundance | b/y observed, T)#
  if (n.obs==0) {
    like.abund <- list(log.like=0)
  } else {
    mean.by <- c( mu.obs[1] + rho*se.obs[1]/se.noise*(mean.noise-mu.noise), mu.obs[2] )
    se.by <- c(sqrt(1-rho^2)*se.obs[1], se.obs[2])
    like.abund <- by.score.abundance2(I.obs = S.match[,2], I.pred = S.theo[,2], mu.beta = mean.by, se.beta = se.by, d = d.var, phi = phi.var, min.match = 0) 
    like.abund$log.like <- as.numeric(like.abund$log.like) + GLM.noise.map(y = log(Noise[,2])-mean.noise, mu.0 = mu.noise, sigma.0 = se.noise, beta = poly.noise, coefs = coef.noise, min.noise = min.noise, max.noise = max.noise, center.y = F) - GLM.noise.map(y = log(S.obs[,2])-mean.noise.2, mu.0 = mu.noise, sigma.0 = se.noise, beta = poly.noise, coefs = coef.noise, min.noise = min.noise, max.noise = max.noise, k = k.noise, center.y = F)
  }
  
  #return results#
  if (!return.all) {return( like.rt + like.ma + like.obs$log.laplace + like.abund$log.like )}
  if (sum(ind.obs)>=3) {
    cor.naive1 <- cor(x.prosit,y.prosit)
    cor.naive2 <- cor(x.prosit[y.prosit>1e-8],y.prosit[y.prosit>1e-8])
    cor.sqrt <- cor(sqrt(x.prosit),sqrt(y.prosit))
    Sim.Score <- sum(sqrt(x.prosit)*sqrt(y.prosit))/sqrt(sum(x.prosit)*sum(y.prosit))
    #Rank.x <- rep(0,length(x.prosit)); Rank.x[x.prosit>1e-8] <- rank(x.prosit[x.prosit>1e-8])
  } else {
    cor.naive1 <- NA; cor.naive2 <- NA; cor.sqrt <- NA; Sim.Score <- NA
  }
  out <- c(like.rt, like.ma, like.obs$log.laplace, like.abund$log.like, like.Obs.Int0, length(ind.obs), sum(ind.obs), n.Int0, n.Int0.obs, cor.naive1, cor.naive2, cor.sqrt, Sim.Score)
  #names(out) <- c("RT", "MA", "Obs", "Int", "Obs.Int0", "N", "N.match", "N.Int0", "N.Int0.obs", "Corr.Prosit", "Corr.Obs", "Corr.Sqrt", "SimScore")
  return(out)
}

#####Intensity 0 peak functions#####
Int0.Obs.score <- function(mz, int, param.beta, n.total=NULL, spline.fit=NULL, min.mz=NULL, max.mz=NULL, ppm.ms2=20, dens=NULL, breaks=NULL) {
  n <- length(mz)
  ind.use <- int >= 1e-8
  k <- sum(ind.use)
  out <- lbeta(param.beta[1]+k, n-k+param.beta[2]) - lbeta(param.beta[1], param.beta[2])
  if (k == 0) {return(out)}
  return( out - Prob.random.match(mz = mz[ind.use], n.total = n.total, spline.fit = spline.fit, min.mz = min.mz, max.mz = max.mz, ppm.ms2 = ppm.ms2, dens = dens, breaks = breaks) )
}

GetInt0Info <- function(S.match, S.theo, ppm.ms2=20) {
  out <- list()
  out$mz <- S.theo[,1]
  out$int <- S.match[,2]; out$int[abs(S.match[,1])>ppm.ms2] <- 0
  return(out)
}

#####pr(b/y ion abundance | b/y ion is observed, T)#####
by.score.abundance <- function(I.obs, I.pred, mu.beta, se.beta, alpha.sigma, gamma.sigma, min.match=3, n.quad=50, out.gauss=NULL, return.log=T, full.integral=F) {
  ind.use <- I.obs > 1e-8
  if (sum(ind.use) < min.match) {return(ifelse(return.log,-Inf,0))}
  I.obs <- I.obs[ind.use]; I.pred <- I.pred[ind.use]
  y <- log(I.obs) - mean(log(I.obs))
  x <- log(I.pred) - mean(log(I.pred))
  SXX <- sum(x^2)
  n <- length(y) - 1
  beta.hat <- sum(x*y)/SXX
  if (beta.hat < 0) {return(ifelse(return.log,-Inf,0))}
  n.sigma2.hat <- sum( (y-beta.hat*x)^2 )
  out <- -n/2*log(2*pi) + lgamma(n/2+alpha.sigma)
  if (!full.integral) {
    if (is.null(out.gauss)) {out.gauss <- statmod::gauss.quad(n = n.quad, kind = "hermite")}
    denom <- gamma.sigma + SXX/2*(mu.beta + se.beta*sqrt(2)*out.gauss$nodes - beta.hat)^2 + n.sigma2.hat/2
    min.denom <- min(denom)
    out <- out - (n/2+alpha.sigma)*log(min.denom) + log( sum(out.gauss$weights/(denom/min.denom)^(n/2+alpha.sigma)) )
    return( ifelse(return.log,out,exp(out)) )
  }
  x.seq <- seq(-5,5,length=1e6)
  y.int <- 1/(gamma.sigma + SXX/2*(mu.beta + se.beta*sqrt(2)*x.seq - beta.hat)^2 + n.sigma2.hat/2)^(n/2+alpha.sigma)
  return( out + log(sum(y.int*exp(-x.seq^2))*(x.seq[2]-x.seq[1])) )
}

#mu.beta is (beta0, beta1)
by.score.abundance2 <- function(I.obs, I.pred, mu.beta, se.beta, d, phi, min.match=3, full.integral=F) {  #I.obs has been normalized by the base peak
  ind.use <- I.obs > 1e-8
  if (sum(ind.use) < min.match) {return(-Inf)}
  out.abund <- list()
  I.obs <- I.obs[ind.use]; I.pred <- I.pred[ind.use]
  y <- log(I.obs)
  alpha.sigma <- d/2; gamma.sigma <- d/phi/2
  if (sum(ind.use) > 1) {
    x <- log(I.pred) - mean(log(I.pred))
    sxx <- sum(x^2)
    n <- length(y)
    beta.hat <- c(mean(y),sum(x*y)/sxx)
    syy <- sum(y^2)
    sigma2.hat <- (syy - sum( beta.hat^2*c(n,sxx) ))/(n-2)
    w1 <- 1/se.beta[1]^2; w2 <- 1/se.beta[2]^2
    a <- mu.beta*c(w1,w2)
    out.abund$post.mean <- ifelse(sum(ind.use)>2,(beta.hat[2]/(sigma2.hat/sxx)+mu.beta[2]*w2)/(1/(sigma2.hat/sxx)+w2),NA)
    if (!full.integral) {
      #Posterior for 1/sigma2 is (s^2 n + d/phi)^{-1}chi^2_{n + d}#
      out.optim <- optim(par = ifelse(n>2, log((n-2+d)/(sigma2.hat*(n-2)+d/phi)), log(phi)), fn = Log.Func, gr = Grad.Log.Func, method = "BFGS", control = list(fnscale=-1), hessian = F, beta=beta.hat, n=n, w=c(w1,w2), a=a, alpha=alpha.sigma, gamma=gamma.sigma, syy=syy, sxx=sxx)
      phi.map <- exp(out.optim$par)
      log.like.map <- n*out.optim$value
      hess <- Hess.Log.Func(phi = phi.map, beta=beta.hat, n=n, w=c(w1,w2), a=a, alpha=alpha.sigma, gamma=gamma.sigma, syy=syy, sxx=sxx)
      s <- ifelse(hess<0, sqrt(-hess), 1)
      int.f <- integrate(f = Like.int, lower = -phi.map*s+1e-8, upper = Inf, C=log.like.map, s=s, phi.0=phi.map, beta=beta.hat, n=n, w=c(w1,w2), a=a, alpha=alpha.sigma, gamma=gamma.sigma, syy=syy, sxx=sxx)
      out.abund$log.like <- -n/2*log(2*pi) + 1/2*sum(log(c(w1,w2))) - 1/2*sum(mu.beta*a) + log.like.map - log(s) + log(int.f$value)
      return(out.abund)
    }
    out <- -n/2*log(2*pi) + lgamma(n/2+alpha.sigma) - (n/2+alpha.sigma)*log(gamma.sigma+syy/2)
    #x.min <- qgamma(1e-8, shape = n/2+alpha.sigma, rate = syy/2+gamma.sigma)
    #x.max <- qgamma(1-1e-8, shape = n/2+alpha.sigma, rate = syy/2+gamma.sigma)
    x.min <- 0; x.max <- 5
    x.seq <- seq(x.min,x.max,length=1e6)
    tmp1 <- x.seq*n + w1; tmp2 <- x.seq*sxx + w2
    log.func <- 1/2*x.seq^2*(n^2/tmp1*beta.hat[1]^2 + sxx^2/tmp2*beta.hat[2]^2) + x.seq*(n/tmp1*beta.hat[1]*a[1] + sxx/tmp2*beta.hat[2]*a[2]) + 1/2*(a[1]^2/tmp1 + a[2]^2/tmp2)
    log.func <- log.func - 1/2*( log(tmp1) + log(tmp2) ) + dgamma(x.seq,n/2+alpha.sigma,syy/2+gamma.sigma,log=T)
    max.func <- max(log.func)
    return( out + max.func + log(sum( exp(log.func-max.func)*(x.seq[2]-x.seq[1]) )) )
  }
  out.optim <- optim(par = log(phi), fn = Log.Func.y1, gr = Grad.Func.y1, method = "BFGS", control = list(fnscale=-1), hessian = F, y=y, tau=1/se.beta[1]^2, mu.0=mu.beta[1], alpha=alpha.sigma, gamma=gamma.sigma)
  phi.map <- exp(out.optim$par)
  log.like.map <- out.optim$value
  hess <- Hess.Func.y1(phi = phi.map, y = y, tau = 1/se.beta[1]^2, mu.0 = mu.beta[1], alpha = alpha.sigma, gamma = gamma.sigma)
  s <- ifelse(hess<0, sqrt(-hess), 1)
  int.f <- integrate(f = Like.int.y1, lower = -phi.map*s+1e-8, upper = Inf, C=log.like.map, s=s, phi.0=phi.map, y = y, tau = 1/se.beta[1]^2, mu.0 = mu.beta[1], alpha = alpha.sigma, gamma = gamma.sigma)
  out.abund$log.like <- -1/2*log(2*pi) - 1/2*log(se.beta[1]) - 1/2*mu.beta[1]^2/se.beta[1]^2 + log.like.map - log(s) + log(int.f$value)
  return(out.abund)
}

Log.Func <- function(theta, beta, n, w, a, alpha, gamma, syy, sxx) {
  phi <- exp(theta)
  out <- (n/2 + alpha - 1)*log(phi) - phi*(gamma + syy/2)
  tmp.n <- c(n, sxx)
  tmp <- phi*tmp.n + w
  g <- phi^2*sum( beta^2*tmp.n^2/tmp ) + 2*phi*sum( beta*a*tmp.n/tmp ) + sum( a^2/tmp )
  h <- sum(log(tmp))
  return( 1/n*(out + 1/2*g - 1/2*h) )
}

Grad.Log.Func <- function(theta, beta, n, w, a, alpha, gamma, syy, sxx) {
  phi <- exp(theta)
  tmp.n <- c(n, sxx)
  tmp <- phi*tmp.n + w
  out <- (n/2 + alpha - 1)/phi - (gamma + syy/2) - 1/2*sum(tmp.n/tmp)
  g1 <- phi*sum( beta^2*tmp.n^2/tmp ) - phi^2/2*sum( beta^2*tmp.n^3/tmp^2 )
  g2 <- sum( beta*a*tmp.n/tmp ) - phi*sum( beta*a*tmp.n^2/tmp^2 )
  g3 <- -1/2*sum( a^2*tmp.n/tmp^2 )
  return( phi/n*(out + g1 + g2 + g3) )
}

Hess.Log.Func <- function(phi, beta, n, w, a, alpha, gamma, syy, sxx) {
  tmp.n <- c(n, sxx)
  tmp <- phi*tmp.n + w
  out <- -(n/2 + alpha - 1)/phi^2 + 1/2*sum(tmp.n^2/tmp^2)
  g1 <- sum( beta^2*tmp.n^2/tmp ) - 2*phi*sum( beta^2*tmp.n^3/tmp^2 ) + phi^2*sum( beta^2*tmp.n^4/tmp^3 )
  g2 <- -2*sum( beta*a*tmp.n^2/tmp^2 ) + 2*phi*sum( beta*a*tmp.n^3/tmp^3 )
  g3 <- sum( a^2*tmp.n^2/tmp^3 )
  return( out+g1+g2+g3 )
}

Log.Func.y1 <- function(theta, y, tau, mu.0, alpha, gamma) {
  phi <- exp(theta)
  return( -1/2*log(phi+tau) + (alpha-1/2)*theta + 1/2*(phi*y+tau*mu.0)^2/(phi+tau) - gamma*phi - phi/2*y^2 )
}

Grad.Func.y1 <- function(theta, y, tau, mu.0, alpha, gamma) {
  phi <- exp(theta)
  out <- -1/2/(phi+tau) + (alpha-1/2)/phi + y*(phi*y + tau*mu.0)/(phi+tau) - 1/2*(phi*y + tau*mu.0)^2/(phi+tau)^2 - gamma - y^2/2
  return( phi*out )
}

Hess.Func.y1 <- function(phi, y, tau, mu.0, alpha, gamma) {
  denom <- phi+tau; post <- phi*y+tau*mu.0
  g1 <- y^2/denom - (y*post - 1/2)/denom^2
  g2 <- -(alpha-1/2)/phi^2
  g3 <- post^2/denom^3 - y*post/denom^2
  return(g1+g2+g3)
}

Like.int <- function(x, C, s, phi.0, beta, n, w, a, alpha, gamma, syy, sxx) {
  phi <- x/s+phi.0
  out <- (n/2 + alpha - 1)*log(phi) - phi*(gamma + syy/2)
  g1 <- 1/2*phi^2*( beta[1]^2*n^2/(phi*n+w[1]) + beta[2]^2*sxx^2/(phi*sxx+w[2]) )
  g2 <- phi*( beta[1]*a[1]*n/(phi*n+w[1]) + beta[2]*a[2]*sxx/(phi*sxx+w[2]) )
  g3 <- 1/2*( a[1]^2/(phi*n+w[1]) + a[2]^2/(phi*sxx+w[2]) )
  h <- -1/2*( log(phi*n+w[1]) + log(phi*sxx+w[2]) )
  return(exp(out + g1 + g2 + g3 + h - C))
}

Like.int.y1 <- function(x, C, s, phi.0, y, mu.0, tau, alpha, gamma) {
  phi <- x/s+phi.0
  out <- -1/2*log(phi+tau) + (alpha-1/2)*log(phi) + 1/2*(phi*y+tau*mu.0)^2/(phi+tau) - gamma*phi - phi/2*y^2
  return(exp(out-C))
}

#####pr(b/y ion is observed | T)#####
by.score.observed <- function(I.obs, I.pred, mu, Omega, min.match=3, positive=F) { #mu is (alpha, delta)
  ind.obs <- I.obs > 1e-8
  if (sum(ind.obs) < min.match) {return(-Inf)}
  out.obs <- list()
  x <- log(I.pred)
  out.obs$flag.opp <- ifelse(sum(ind.obs)==0,0,ifelse(max(x[ind.obs]) < min(x[!ind.obs]),1,0))
  out.obs$flag.neg <- ifelse(sum((x-mean(x))*as.numeric(ind.obs)) < 0,1,0)
  
  y <- as.numeric(ind.obs)
  out <- optim(par = mu, fn = Penalized.log.like, gr = Penalized.grad, method = "BFGS", control = list(fnscale=-1), x=x, y=y, Omega=Omega, mu=mu, positive=positive)
  theta <- out$par
  X <- cbind(x,rep(1,length(x)))
  if (positive) {
    p <- Expit(exp(theta[1])*x + theta[2])
    grad.h <- diag(c(exp(theta[1]),1))
    A <- grad.h%*%(t(X*(p*(1-p)))%*%X)%*%grad.h
    A <- A + Omega; A[1,1] <- A[1,1] - exp(theta[1])*sum(x*(y-p))
  } else {
    p <- Expit(theta[1]*x + theta[2])
    #A <- t(X*(p*(1-p)))%*%X + diag(2)
    A <- t(X*(p*(1-p)))%*%X + Omega
  }
  out.obs$theta <- theta
  out.obs$log.laplace <- out$value - 1/2*log(A[1,1]*A[2,2]-A[1,2]^2) + 1/2*log(Omega[1,1]*Omega[2,2]-Omega[1,2]^2)
  return(out.obs)
}

Expit <- function(x) {
  return( 1/(1+exp(-x)) )
}

Penalized.log.like <- function(theta, Omega, mu, x, y, positive) {
  if (positive) {
    p <- Expit(exp(theta[1])*x + theta[2])
  } else {
    p <- Expit(theta[1]*x + theta[2])
  }
  theta <- theta - mu
  ifelse(sum(y==1)==0,0,sum(log(p[y==1]))) + ifelse(sum(y==0)==0,0,sum(log(1-p[y==0]))) - 1/2*sum( theta*as.vector(Omega%*%theta) )
}

Penalized.grad <- function(theta, Omega, mu, x, y, positive) {
  if (positive) {
    p <- Expit(exp(theta[1])*x + theta[2])
    g <- c( sum(x*(y-p))*exp(theta[1]), sum(y-p) )
  } else {
    p <- Expit(theta[1]*x + theta[2])
    g <- c( sum(x*(y-p)), sum(y-p) )
  }
  return( g - as.vector(Omega%*%(theta-mu)) )
}

#######Noise model#######
GLM.noise.map <- function(y, mu.0, sigma.0, beta, coefs, min.noise, max.noise, k=1, center.y=T) {
  if (center.y) {
    mu.hat <- mean(y); v.hat <- var(y)/length(y)
    mu.map <- (mu.hat/v.hat + mu.0/sigma.0^2)/(1/v.hat + 1/sigma.0^2)
    y <- y - mu.map
  }
  y[y < min.noise] <- min.noise; y[y > max.noise] <- max.noise
  sum(as.vector(poly(x = y, degree = length(beta)-1, coefs = coefs)%*%beta[2:length(beta)]) + beta[1] - log(k))
}

########Mass accuracy########
#Matched ions
Mass.Accuracy.match <- function(delta=NULL, mz.obs=NULL, mz.theo=NULL, log.int=NULL, mu=0, sigma=c(7.889799, 1.126632, 3.035488), lambda=c(0.3766329, 0.2205460, 0.4028211), beta.int=NULL, bound.int.MA=NULL) {
  if (is.null(delta)) {delta <- (mz.obs/mz.theo - 1)*10^6 - mu}
  if (length(delta) > 1) {
    if (!is.null(log.int) && !is.null(beta.int) && !is.null(bound.int.MA)) {
      delta[delta<bound.int.MA[1]] <- bound.int.MA[1]; delta[delta>bound.int.MA[2]] <- bound.int.MA[2]
      sigma <- -sort(-sigma)
      pi.large <- 1/(1 + exp(-( poly(x = log.int, degree = length(beta.int)-1, raw = T)%*%beta.int[2:length(beta.int)] + beta.int[1] )))
      k <- as.vector(cbind(pi.large,1-pi.large)%*%(1-2*pnorm(q = -20, mean = 0, sd = sigma)))
      return(sum(log( pi.large*dnorm(delta,0,sigma[1])+(1-pi.large)*dnorm(delta,0,sigma[2]) )-log(k)))
    }
    k <- sum(lambda*(1-2*pnorm(q=-20,mean=0,sd=sigma)))
    tmp <- sum(log(rowSums( sapply(1:length(sigma),function(i){lambda[i]*dnorm(delta,0,sigma[i])}) ))-log(k))
    if (is.null(mz.theo)) {return(tmp)}
    return( tmp - sum(log(mz.theo/10^6)) )
  }
  if (!is.null(log.int) && !is.null(beta.int) && !is.null(bound.int.MA)) {
    delta <- min(max(bound.int.MA[1],delta),bound.int.MA[2])
    sigma <- -sort(-sigma)
    pi.large <- 1/(1 + exp(-( sum(poly(x = log.int, degree = length(beta.int)-1, raw = T)*beta.int[2:length(beta.int)]) + beta.int[1] )))
    k <- sum(c(pi.large,1-pi.large) * (1-2*pnorm(q = -20, mean = 0, sd = sigma)))
    return(log(sum(c(pi.large,1-pi.large)*dnorm(delta,0,sigma))) - log(k))
  }
  k <- sum(lambda*(1-2*pnorm(q=-20,mean=0,sd=sigma)))
  tmp <- log(sum( lambda*dnorm(delta,0,sigma) ))-log(k)
  if (is.null(mz.theo)) {return(tmp)}
  return(tmp - log(mz.theo/10^6))
}

Mass.Accuracy.noise <- function(delta, mz.theo=NULL, ppm.max=20, theta=c(-29.808050, 5.242391), sigma=3.097397, bound=NULL) {
  if (is.null(mz.theo)) {return(-length(delta)*log(2*ppm.max))}
  x <- log(mz.theo); x[x<bound[1]] <- bound[1]; x[x>bound[2]] <- bound[2]
  if (length(x)>1) {
    pi.prob <- 1/(1+exp(-as.vector(poly(x = x, degree = length(theta)-1, raw = T)%*%theta[2:length(theta)] + theta[1])))
  } else {
    pi.prob <- 1/(1+exp(-(sum(poly(x = x, degree = length(theta)-1, raw = T)*theta[2:length(theta)]) + theta[1])))
  }
  return( sum(log( pi.prob/2/ppm.max + (1-pi.prob)*dnorm(delta,0,sigma) )) )
}

Mass.Accuracy.score <- function(delta, log.int, mz.theo, beta.obs, beta.noise, sigma.obs, bound.log.int, bound.log.mz, ppm.max=20) {
  ind.small <- log.int < bound.log.int[1]; ind.large <- log.int > bound.log.int[2]
  if (sum(ind.small)>0) {log.int[ind.small] <- bound.log.int[1]}; if (sum(ind.large)>0) {log.int[ind.large] <- bound.log.int[2]}
  log.mz.theo <- log(mz.theo); rm(mz.theo)
  ind.small <- log.mz.theo < bound.log.mz[1]; ind.large <- log.mz.theo > bound.log.mz[2]
  if (sum(ind.small)>0) {log.mz.theo[ind.small] <- bound.log.mz[1]}; if (sum(ind.large)>0) {log.mz.theo[ind.large] <- bound.log.mz[2]}
  
  sigma.obs <- -sort(-sigma.obs)
  if (length(delta)==1) {
    pi.large <- 1/(1+exp(-( sum(poly(x = log.int, degree = length(beta.obs)-1, raw = T, simple = T)*beta.obs[2:length(beta.obs)]) + beta.obs[1] )))
    pi.unif <- 1/(1+exp(-( sum(poly(x = log.mz.theo, degree = length(beta.noise)-1, raw = T, simple = T)*beta.noise[2:length(beta.noise)]) + beta.noise[1] )))
  } else {
    pi.large <- 1/(1 + exp(-as.vector( poly(x = log.int, degree = length(beta.obs)-1, raw = T, simple = T)%*%beta.obs[2:length(beta.obs)] + beta.obs[1] )))
    pi.unif <- 1/(1 + exp(-as.vector( poly(x = log.mz.theo, degree = length(beta.noise)-1, raw = T, simple = T)%*%beta.noise[2:length(beta.noise)] + beta.noise[1] )))
  }
  dens.obs <- (pi.large*dnorm(delta,0,sigma.obs[1]) + (1-pi.large)*dnorm(delta,0,sigma.obs[2])) / (pi.large*(1-2*pnorm(-ppm.max,0,sigma.obs[1])) + (1-pi.large)*(1-2*pnorm(-ppm.max,0,sigma.obs[2])))
  return(sum(log(dens.obs) - log(pi.unif/2/ppm.max + (1-pi.unif)*dens.obs)))
}

########Random match with noise########

Prob.random.match <- function(mz, n.total, spline.fit, min.mz, max.mz, ppm.ms2=20, dens=NULL, breaks=NULL) {  #mz is theoretical m/z
  n.mz <- length(mz)
  out <- log(2*ppm.ms2/10^6*mz/(max.mz-min.mz))
  out.n <- sum(log(n.total:(n.total-n.mz+1)))
  if (is.null(dens) || is.null(breaks)) {
    ind.spline <- which(mz>=spline.fit$Boundary.knots[1] & mz<=spline.fit$Boundary.knots[2])
    if (length(ind.spline)==0) {return(out.n+sum(out))}
    X <- splines::bs(mz[ind.spline], degree = spline.fit$degree, knots = spline.fit$knots, Boundary.knots = spline.fit$Boundary.knots)
    out[ind.spline] <- pmax(as.vector(spline.fit$beta[1] + X%*%spline.fit$beta[2:length(spline.fit$beta)]),out[ind.spline])
    return(out.n+sum(out))
  }
  ind.small <- which(mz <= max(breaks))
  if (length(ind.small) > 0) {out[ind.small] <- pmax(out[ind.small], dens[findInterval(mz[ind.small],breaks)])}
  ind.spline <- which(mz>=spline.fit$Boundary.knots[1] & mz<=spline.fit$Boundary.knots[2] & mz > max(breaks))
  if (length(ind.spline) > 0) {
    X <- splines::bs(mz[ind.spline], degree = spline.fit$degree, knots = spline.fit$knots, Boundary.knots = spline.fit$Boundary.knots)
    out[ind.spline] <- pmax(as.vector(spline.fit$beta[1] + X%*%spline.fit$beta[2:length(spline.fit$beta)]),out[ind.spline])
  }
  return(out.n+sum(out))
}


########RT########
RT.score <- function(rt, prosit.rt, lambda=c(0.82309646, 0.17690354), v=c(0.02856569, 0.71997677), M=90, beta=c(-2.458869436302792088611113285879,0.023796245986352743129188525018,-0.000021477680631739091239237882)) {
  delta <- log(rt/M) - log(1-rt/M) - sum(beta*c(1,prosit.rt,prosit.rt^2))
  lambda <- lambda/sum(lambda)
  log.dens <- dnorm(x = delta, mean = 0, sd = sqrt(v), log = T)
  C <- max(log.dens)
  return( C + log( sum(lambda*exp( log.dens-C )) ) )
}

####Precursor matching####
Match.Prec <- function(mz.prec, mz.theory, delta.ppm=10, mu.ppm=0) {
  delta <- delta.ppm/10^6
  vec <- c(-Inf,mz.theory,Inf)
  n <- length(vec) - 2
  Int <- findInterval(x = mz.prec, vec = vec)
  return(t(sapply(1:length(Int),function(i){
    ind <- Int[i]
    m.hat <- mz.prec[i]
    ind.L <- max(ind-1,1); ind.R <- min(ind,n)
    #Lower boundary#
    go <- 1
    while(go) {
      if (ind.L <= 1) {
        go <- 0
      } else {
        if (abs(m.hat/mz.theory[ind.L-1]-1-mu.ppm)>delta) {
          go <- 0
        } else {
          ind.L <- ind.L - 1
        }
      }
    }
    #Upper boundary#
    go <- 1
    while(go) {
      if (ind.R >= n) {
        go <- 0
      } else {
        if (abs(m.hat/mz.theory[ind.R+1]-1-mu.ppm)>delta) {
          go <- 0
        } else {
          ind.R <- ind.R + 1
        }
      }
    }
    ind.final <- (ind.L:ind.R)[which(abs(m.hat/mz.theory[ind.L:ind.R]-1-mu.ppm)<=delta)]
    if (length(ind.final)==0) {return(rep(NA,2))}
    return(ind.final[c(1,length(ind.final))])
  })))
}

####Ion matching####
Delta.Mass <- function(mass.find, masses, sorted=F, charge.find=NULL, charge.masses=NULL, mu.ppm=0) {
  if (!sorted) {   #length(mass.find) < length(masses)
    order.mass <- order(masses)
    masses <- masses[order.mass]
    order.find <- order(mass.find)
    if (!is.null(charge.masses) && !is.null(charge.find)) {
      charge.masses <- charge.masses[order.mass]
      charge.find <- charge.find[order.find]
    }
    rank.find <- rank(mass.find)
    mass.find <- mass.find[order.find]
  }
  masses <- c(0, masses, Inf)
  int <- findInterval(x = mass.find, vec = masses)
  delta.1 <- (masses[int]/mass.find-1-mu.ppm)*10^6; delta.2 <- (masses[int+1]/mass.find-1-mu.ppm)*10^6
  delta.1[int==1] <- Inf
  ind.min <- int; delta <- delta.1
  ind.change <- abs(delta.2)<abs(delta.1)
  if (sum(ind.change)>0) {ind.min[ind.change]<-ind.min[ind.change]+1; delta[ind.change]<-delta.2[ind.change]}
  #ind.min <- apply(cbind(abs(delta.1),abs(delta.2),int),1,function(x){ifelse(which.min(x[1:2])==1,x[3],x[3]+1)})
  #delta <- apply(cbind(delta.1,delta.2),1,function(x){x[which.min(abs(x))]})
  if (!is.null(charge.masses) && !is.null(charge.find)) {
    tmp <- charge.find==charge.masses[ind.min-1]; tmp[is.na(tmp)] <- T
    delta[!tmp] <- 1e16
  }
  if (!sorted) {return(list(delta=delta[rank.find],ind=order.mass[ind.min[rank.find]-1]))}
  return(list(delta=delta,ind=ind.min-1))
}
