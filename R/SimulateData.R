library(mzR)
library(rlist)
###Simulate data###
#Add signal to di-signaled real data.

SimulateData.new <- function(ProcessedDataPath, NoiseInt.path="/Users/Chris/Desktop/UniversityofPittsburgh/Projects/YiShi/output/TrainingData/TrainingParameters/NoiseInt/", rt.path="/Users/Chris/Desktop/UniversityofPittsburgh/Projects/YiShi/output/TrainingData/TrainingParameters/RetentionTime/RetentionTime_new.txt", byInt.path="/Users/Chris/Desktop/UniversityofPittsburgh/Projects/YiShi/output/TrainingData/TrainingParameters/BYInt/", byMA.path="/Users/Chris/Desktop/UniversityofPittsburgh/Projects/YiShi/output/TrainingData/TrainingParameters/BYMassAccuracy/BYMassAccuracy_new.txt", byObs.path="/Users/Chris/Desktop/UniversityofPittsburgh/Projects/YiShi/output/TrainingData/TrainingParameters/BYObserved/", append.file="", df.by=Inf) {
  scale.type <- "Quant90"
  mzml.path <- ""
  tmp <- load(ProcessedDataPath); rm("Sim.N0")
  database.PSM <- sapply(data.F,function(x){paste0(x$Pep,".",x$z)}); ind.keep <- !duplicated(database.PSM); data.F <- data.F[ind.keep]; database.PSM <- database.PSM[ind.keep]
  if (grepl(pattern = "\\.gz$", x = mzml.path, ignore.case = F) && !is.null(file.spec)) {
    mzml.path.gz <- mzml.path
    tmp.path <- paste0("highph/",gsub(pattern = "\\.mzML$", replacement = "", x = gsub(pattern = "^.*/", replacement = "", x = file.spec), ignore.case = T),".mzML")
    mzml.path <- paste0("../data/",tmp.path)
    if (!file.exists(mzml.path)) {
      untar(tarfile = mzml.path.gz, files = tmp.path, list = F, exdir = "../data")
    }
  }
  missing.params <- NULL; file.spec <- NULL; data.F2 <- NULL; True.Results.F2 <- NULL; charge.bin <- list(2,3,4:10); return.MA <- F; use.all <- F; scale.data <- F; simulate.noise <- F
  
  if (!is.null(missing.params)) {
    pi0 <- missing.params$pi0; mz.shift <- missing.params$mz.shift
  } else {
    pi0 <- NULL; mz.shift <- NULL
  }
  
  #Noise Int#
  NoiseInt.path <- gsub(pattern = "/$", replacement = "", x = NoiseInt.path)
  param.Noise.Int <- list()
  #for (c in 1:length(charge.bin)) {
  #  Noise.tmp <- readLines(con = paste0(NoiseInt.path,"/","NoiseInt.",c,append.file,".txt"))
  #  param.Noise.Int[[c]] <- list()
  #  param.Noise.Int[[c]]$beta <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^beta", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  #  param.Noise.Int[[c]]$norm2 <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^norm2", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  #  param.Noise.Int[[c]]$alpha <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^alpha", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  #  param.Noise.Int[[c]]$bound <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^min", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  #  if (simulate.noise) {
  #   param.Noise.Int[[c]]$hist <- Hist.Noise(beta = param.Noise.Int[[c]]$beta, norm2 = param.Noise.Int[[c]]$norm2, alpha = param.Noise.Int[[c]]$alpha, min = param.Noise.Int[[c]]$bound[1], param.Noise.Int[[c]]$bound[2])
  #  }
    #}
  
  #RT#
  RT.tmp <- readLines(con = rt.path)
  param.RT <- list()
  param.RT$beta <- as.numeric(strsplit(x = strsplit(x = RT.tmp[1], split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.RT$lambda <- as.numeric(strsplit(x = strsplit(x = RT.tmp[2], split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.RT$lambda <- param.RT$lambda/sum(param.RT$lambda)
  param.RT$sigma <- sqrt(as.numeric(strsplit(x = strsplit(x = RT.tmp[3], split = ":\\s*")[[1]][2], split = ",")[[1]]))
  
  #BY Int#
  byInt.path <- gsub(pattern = "/$", replacement = "", x = byInt.path)
  param.BYInt <- list()
  for (c in 1:length(charge.bin)) {
    Int.tmp <- readLines(con = paste0(byInt.path,"/","BYInt.",c,append.file,".txt"))
    param.BYInt[[c]] <- list()
    param.BYInt[[c]]$beta <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^beta", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$beta) <- c("beta", "se")
    param.BYInt[[c]]$mu <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^mu,", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$mu) <- c("mu", "se", "rho")
    #if (scale.type=="Quant90") {param.BYInt[[c]]$mu[3] <- 0}
    param.BYInt[[c]]$mu.noise <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^mu\\.noise", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$mu.noise) <- c("mu.noise", "se.noise")
    param.BYInt[[c]]$var <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^d", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$var) <- c("d", "phi")
  }
  
  #BY Observed#
  byObs.path <- gsub(pattern = "/$", replacement = "", x = byObs.path)
  param.BYObs <- list()
  for (c in 1:length(charge.bin)) {
    Obs.tmp <- readLines(con = paste0(byObs.path,"/","BYObs.",c,append.file,".txt"))
    param.BYObs[[c]] <- list()
    param.BYObs[[c]]$ad <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^alpha", x = Obs.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYObs[[c]]$ad) <- c("alpha", "delta")
    param.BYObs[[c]]$Omega <- matrix(as.numeric(strsplit(x = strsplit(x = grep(pattern = "^Omega", x = Obs.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]]), nrow=2)
  }
  
  #MA Observed#
  MA.path <- gsub(pattern = "/$", replacement = "", x = byMA.path)
  MA.tmp <- readLines(con = MA.path)
  param.MA <- list()
  param.MA$beta <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^beta", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.MA$sigma <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^sigma", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.MA$bound <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^min", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  
  #Get real data#
  out <- list()
  database.search <- Sim.0$database.search$Results
  Sequest <- Sim.0$Sequest
  #Sequest <- data.frame(read.csv(Sequest.path, stringsAsFactors = F, check.names = F))
  min.score.use <- 15
  #ind.train <- which(Sequest$PSM.Ambiguity!="Rejected"); tab.sequ <- table(Sequest$First.Scan[ind.train]); ind.train <- ind.train[!(Sequest$First.Scan[ind.train]%in%as.numeric(names(tab.sequ[tab.sequ>1])))]
  #Sequest <- Sequest[-ind.train,]
  Sequest$Annotated.Sequence <- gsub(pattern = "c", replacement = "C", x = Sequest$Annotated.Sequence, ignore.case = F)
  Sequest$length <- sapply(Sequest$Annotated.Sequence,function(x){length(strsplit(x,"")[[1]])})
  Sequest$PSM <- paste0(Sequest$Annotated.Sequence,".",Sequest$Charge)
  Sequest <- Sequest[Sequest$PSM%in%database.PSM,]
  if (!is.null(Sim.0)) {
    Map.Peps <- Map.IDs(Sequest.scan = Sequest$First.Scan, Sequest.pep = Sequest$PSM, Sim0.scan = Sim.0$Headers.sim$scan, Results = database.search, min.score = min.score.use, use.all = use.all)
  } else {
    Map.Peps <- Map.IDs(Sequest.scan = Sequest$First.Scan, Sequest.pep = Sequest$PSM, Sim0.scan = NULL, Results = NULL, min.score = min.score.use, use.all = use.all)
  }
  if (!is.null(Map.Peps$All.results)) {
    ttt <- c(unlist(Map.Peps$All.results),Map.Peps$AllPeps); ttt <- unique(ttt[!is.na(ttt)])
    data.F <- data.F[database.PSM%in%ttt]; database.PSM <- database.PSM[database.PSM%in%ttt]; rm(ttt)
  } else {
    data.F <- data.F[database.PSM%in%Map.Peps$AllPeps]; database.PSM <- database.PSM[database.PSM%in%Map.Peps$AllPeps]
  }
  All.scans <- Map.Peps$Scans
  if (is.null(Sim.0)) {
    Spectra.prosit <- lapply(data.F,function(x){tmp <- Clean.prosit(Prosit = x$mat, add = T); tmp <- tmp[tmp[,2]>=1e-8,]; return(tmp[order(tmp[,1]),])})
  } else {
    Spectra.prosit <- lapply(data.F,function(x){tmp <- x$mat[x$mat[,2]>=1e-8,]; return(tmp[order(tmp[,1]),])})
  }
  irt.prosit <- sapply(data.F,function(x){x$irt})
  if (is.null(Sim.0)) {
    mZ <- mzR::openMSfile(filename = mzml.path)
    Spectra <- peaks(mZ)
    Headers <- header(mZ); rm(mZ)
    Spectra[All.scans] <- lapply(All.scans, function(s){
      (Deisotope.workhorse(Spectrum = Spectra[[s]], Prec.mz = Headers$precursorMZ[s], Iso.mz = Headers$isolationWindowTargetMZ[s], z = Headers$precursorCharge[s], PIW.left = Headers$isolationWindowLowerOffset[s], PIW.right = Headers$isolationWindowUpperOffset[s], ppm = 15, ppm.filter = 30, max.addition = 2, remove.precursor = F))[,1:3]
    })
    out$Spectra <- Spectra; out$Headers <- Headers; if (is.null(out$Headers$scan)) {out$Headers$scan <- 1:nrow(Headers)}
    rm(Spectra,Headers)
  } else {
    out$Spectra <- Sim.0$Spectra; out$Headers <- Sim.0$Headers; if (is.null(out$Headers$scan)) {out$Headers$scan <- 1:nrow(Headers)}
  }
  
  #Simulate data#
  out$Pep <- rep(NA, length(All.scans))
  out$Charge <- rep(NA, length(All.scans))
  out$Spectra.sim <- out$Spectra[All.scans]; out$Headers.sim <- out$Headers[All.scans,]
  out$BYObs <- matrix(NA, nrow=length(All.scans), ncol=2); colnames(out$BYObs) <- c("Alpha", "Delta")
  out$BYInt <- matrix(NA, nrow=length(All.scans), ncol=2); colnames(out$BYInt) <- c("mu", "beta")
  out$NoiseInt <- rep(NA, length(All.scans))
  out$Sigma <- rep(NA, length(All.scans))
  out$Scale <- rep(NA, length(All.scans))
  if (!is.null(pi0)) {
    out$Hg <- rbinom(n = length(All.scans), size = 1, prob = 1-pi0)
  } else {
    out$Hg <- NA
  }
  if (return.MA) {out$MA <- vector(mode = "list", length = length(All.scans))}
  for (j in 1:length(All.scans)) {
    s <- All.scans[j]
    if (!is.na(out$Hg[1])) {
      mz.null <- Random.null(mz = Headers$precursorMZ[s], mz.all = Headers$precursorMZ)
      next
    }
    ind.sequest.s <- which(Sequest$First.Scan==s)
    ind.prosit.s <- match(Map.Peps$Peps.Sequ[[j]],database.PSM)
    if (is.null(Map.Peps$All.results)) {
      tmp.s <- sample(x = Map.Peps$Peps.Sequ[[j]], size = 1, replace = F)
    } else {
      ttt <- c(Map.Peps$All.results[[j]],Map.Peps$AllPeps[[j]])
      tmp.s <- sample(x = unique(ttt[!is.na(ttt)]), size = 1, replace = F); rm(ttt)
      ind.prosit.s <- unique(c(ind.prosit.s, which(database.PSM==tmp.s)))
    }
    out$Charge[j] <- as.numeric(gsub(pattern = "^.*\\.([0-9]+)$", replacement = "\\1", x = tmp.s))
    out$Pep[j] <- gsub(pattern = "\\.[0-9]+$", replacement = "", x = tmp.s)
    c <- which(sapply(charge.bin,function(x){out$Charge[j]%in%x}))
    Log.int.prosit <- log(Spectra.prosit[[which(database.PSM==tmp.s)]][,2])
    mz.prosit <- Spectra.prosit[[which(database.PSM==tmp.s)]][,1]
    tmp.s.copy <- tmp.s
    
    if (is.null(Map.Peps$Peps.mine)) {
      tmp.s <- Remove.signal(S = out$Spectra[[s]], header = out$Headers[s,], mz.vec = rlist::list.rbind(Spectra.prosit[ind.prosit.s]), deisotope = F)
    } else {
      if (is.na(Map.Peps$Peps.mine[[j]][1])) {
        tmp.s <- Remove.signal(S = out$Spectra[[s]], header = out$Headers[s,], mz.vec = rlist::list.rbind(Spectra.prosit[ind.prosit.s]), deisotope = F)
      } else {
        tmp.s <- Remove.signal(S = out$Spectra[[s]], header = out$Headers[s,], mz.vec = rbind(rlist::list.rbind(Spectra.prosit[ind.prosit.s]),rlist::list.rbind(Spectra.prosit[match(Map.Peps$Peps.mine[[j]],database.PSM)])), deisotope = F)
      }
    }
    Noise.s <- rbind(tmp.s$E)[,1:2]
    if (scale.type == "base-peak") {scale.s <- tmp.s$basepeak}
    if (scale.type == "Quant90") {scale.s <- quantile(c(rbind(tmp.s$E)[,2],rbind(tmp.s$O)[,2]),0.9)}
    out$Scale[j] <- scale.s
    
    #BY Observed#
    Obs.params <- param.BYObs[[c]]$ad + as.vector(SQRT.mat(solve(param.BYObs[[c]]$Omega))%*%rnorm(2))
    out$BYObs[j,] <- Obs.params
    ind.obs <- rbinom(n = length(Log.int.prosit), size = 1, prob = 1/(1+exp(-( Obs.params[2]+Obs.params[1]*Log.int.prosit )))) == 1
    
    #Noise#
    mu.noise.s <- mean(log(Noise.s[,2]/scale.s))
    out$NoiseInt[j] <- mu.noise.s
    if (simulate.noise) {
      mu.all <- Sim.Mu(mu.beta = param.BYInt[[c]]$mu[1], mu.noise = param.BYInt[[c]]$mu.noise[1], sigma.beta = param.BYInt[[c]]$mu[2], sigma.noise = param.BYInt[[c]]$mu.noise[2], rho = param.BYInt[[c]]$mu[3])
      ind.noise.s <- apply(rmultinom(n = length(Noise.s[,2]), size = 1, prob = param.Noise.Int[[c]]$hist$probs), 2, function(x){which(x==1)})
      random.noise.s <- (param.Noise.Int[[c]]$hist$breaks[ind.noise.s+1]-param.Noise.Int[[c]]$hist$breaks[ind.noise.s])*runif(length(Noise.s[,2])) + param.Noise.Int[[c]]$hist$breaks[ind.noise.s]
      Noise.s[,2] <- exp(mu.all[2] + random.noise.s)
      #scale.s <- max(scale.s, Noise.s[,2])
      Int.mu <- mu.all[1]
    } #else {
    #Noise.s[,2] <- Noise.s[,2]/scale.s
    #}
    
    #BY Int#
    Int.beta <- param.BYInt[[c]]$beta[1] + param.BYInt[[c]]$beta[2]*rnorm(1)
    if (!simulate.noise) {Int.mu <- param.BYInt[[c]]$mu[1] + param.BYInt[[c]]$mu[3]*param.BYInt[[c]]$mu[2]/param.BYInt[[c]]$mu.noise[2]*(mu.noise.s - param.BYInt[[c]]$mu.noise[1]) + param.BYInt[[c]]$mu[2]*sqrt(1-param.BYInt[[c]]$mu[3]^2)*rnorm(1)}
    out$BYInt[j,] <- c(Int.mu, Int.beta)
    out$Sigma[j] <- 1/(param.BYInt[[c]]$var[2]/param.BYInt[[c]]$var[1]*rchisq(1,param.BYInt[[c]]$var[1]))
    if (scale.data) {
      Int.obs <- scale.s*exp( Int.beta*(Log.int.prosit[ind.obs]-mean(Log.int.prosit[ind.obs])) + Int.mu + sqrt(out$Sigma[j])*Sim.Tdist(n = sum(ind.obs), df = df.by) )
    } else {
      Int.obs <- exp( Int.beta*(Log.int.prosit[ind.obs]-mean(Log.int.prosit[ind.obs])) + Int.mu + sqrt(out$Sigma[j])*Sim.Tdist(n = sum(ind.obs), df = df.by) )
    }
    
    #BY mz#
    mz.obs <- Simulate.mz(mz.mean = mz.prosit[ind.obs], x = log(Int.obs)+ifelse(scale.data,0,log(scale.s)), theta = param.MA$beta, sigma = sort(param.MA$sigma), df = Inf)
    if (return.MA) {out$MA[[j]] <- cbind(log(Int.obs)+ifelse(scale.data,0,log(scale.s)),10^6*(mz.obs/mz.prosit[ind.obs]-1))}
    
    if (scale.data) {
      out$Spectra.sim[[j]] <- rbind(Noise.s, cbind(mz.obs,Int.obs))
    } else {
      Noise.s[,2] <- Noise.s[,2]/scale.s
      out$Spectra.sim[[j]] <- rbind(Noise.s, cbind(mz.obs,Int.obs))
    }
    tmp.sim <- out$Spectra.sim[[j]][order(out$Spectra.sim[[j]][,1]),]
    out$Spectra.sim[[j]] <- tmp.sim[tmp.sim[,2]>=1e-8,]
    
    #RT#
    out$Headers.sim$retentionTime[j] <- Sim.RT(irt = irt.prosit[which(database.PSM==tmp.s.copy)], beta = param.RT$beta, sigma = param.RT$sigma, lambda = param.RT$lambda)
  }
  if (!is.null(data.F2)) {
    out$data.F2 <- data.F2[which(!sapply(data.F2,function(x){paste0(x$Pep,".",x$z)}) %in% database.PSM)]
    if (!is.null(True.Results.F2)) {
      tmp <- unlist(lapply(True.Results.F2,function(x){
        if (is.na(x$FM[1])) {return(NA)}
        x$labels <- !x$labels
        if (sum(x$labels)==0) {return(NA)}
        return(x$Peptides[x$labels])
      }))
      tmp <- tmp[!is.na(tmp)]; tmp <- unique(tmp)
      out$data.F2 <- out$data.F2[sapply(out$data.F2,function(x){x$Pep})%in%tmp]; rm(tmp)
    }
    out$data.F2 <- lapply(out$data.F2,function(x){x$Pep<-paste0(".",x$Pep); x$mat <- Clean.prosit(Prosit = x$mat, add = T); return(x)})
  }
  return(out)
}

Logit.RT <- function(x, max=90, type=c("reg","inv")) {
  type <- match.arg(type, c("reg","inv"))
  if (type=="reg") {
    return( log(x/max) - log(1-x/max) )
  }
  return( max*exp(x)/(1+exp(x)) )
}

SQRT.mat <- function(X) {
  s <- eigen(x = X, symmetric = T)
  return( s$vectors%*%diag(sqrt(s$values))%*%t(s$vectors) )
}

Remove.signal <- function(S, header, mz.vec=NULL, ppm=20, deisotope=T, mu.ppm=0) {
  mu.ppm <- mu.ppm/10^6
  Prec.mz <- header$precursorMZ
  if (NROW(mz.vec) > 1) {
    mz.vec <- mz.vec[,1]
  }
  mz.vec <- unique(mz.vec)
  if (deisotope) {
    if (NROW(S) > 20) {
      S.out <- Deisotope.workhorse(Spectrum = S, Prec.mz = Prec.mz, Iso.mz = header$isolationWindowTargetMZ, z = header$precursorCharge, PIW.left = header$isolationWindowLowerOffset, PIW.right = header$isolationWindowUpperOffset, ppm = 15, max.addition = 2, ppm.filter = 30, remove.precursor = F)
      S <- S.out[,1:2]
    }
  }
  delta <- ppm/10^6
  ind.prec <- which(abs(S[,1]/Prec.mz-1)<=delta)
  if (length(ind.prec) > 0) {
    S <- S[-ind.prec,]
  }
  ind.remove <- sapply(mz.vec, function(m) {
    ind <- which(abs(S[,1]/m-1-mu.ppm)<=delta)
    return( ifelse(length(ind)==0, NA, ind[which.max(S[ind,2])]) )
  })
  ind.remove <- ind.remove[!is.na(ind.remove)]
  if (length(ind.remove)==0) {return(list(E=S,O=NULL))}
  return(list(E=S[-ind.remove,],O=rbind(S[ind.remove,]),basepeak=max(S[,2])))
}

Map.IDs <- function(Sequest.scan, Sequest.pep, Sim0.scan, Results, min.score=15, use.all=F) {
  Scans <- sort(unique(Sequest.scan))
  out.sequ <- lapply(Scans,function(scan){sequ.pep <- unique(Sequest.pep[which(Sequest.scan==scan)])})
  if (!is.null(Sim0.scan)) {
    out.mine <- lapply(1:length(Scans),function(j){
      scan <- Scans[j]
      my.pep <- unlist(lapply(which(Sim0.scan==scan),function(i){
        x <- Results[[i]]$FM
        if (is.na(x[1])){return(NA)}
        ind.s <- which(rowSums(rbind(x[,1:4])) >= min.score)
        if (length(ind.s)==0){return(NA)}
        return(rownames(x)[ind.s])
      }))
      my.pep <- my.pep[!is.na(my.pep)]; if (length(my.pep)==0) {return(NA)}
      my.pep <- my.pep[!my.pep%in%out.sequ[[j]]]; if (length(my.pep)==0) {return(NA)}
      return(unique(my.pep))
    })
    all <- c(unlist(out.sequ),unlist(out.mine))
  } else {
    out.mine <- NULL; all <- unlist(out.sequ)
  }
  if (use.all) {
    All.results <- lapply(1:length(Scans),function(j){
      scan <- Scans[j]
      my.pep <- unlist(lapply(which(Sim0.scan==scan),function(i){
        x <- Results[[i]]$FM
        if (is.na(x[1])){return(NA)}
        return(rownames(x))
      }))
      my.pep <- my.pep[!is.na(my.pep)]; if (length(my.pep)==0) {return(NA)}
      my.pep <- my.pep[!my.pep%in%out.sequ[[j]]]; if (length(my.pep)==0) {return(NA)}
      return(unique(my.pep))
    })
  } else {
    All.results <- NULL
  }
  return(list(Scans=Scans, Peps.Sequ=out.sequ, Peps.mine=out.mine, AllPeps=unique(all[!is.na(all)]), All.results=All.results))
}

Sim.RT <- function(irt, beta, sigma, lambda) {
  inv <- sum(beta*c(1,irt,irt^2)) + ifelse(rbinom(n = 1, size = 1, prob = lambda[1]) == 1, sigma[1]*rnorm(1), sigma[2]*rnorm(1))
  return(60*Logit.RT(x = inv, type = "inv"))
}

Hist.Noise <- function(beta, norm2, alpha, min, max, n.breaks=1e3) {
  breaks <- seq(min, max, length=n.breaks)
  probs <- sapply(1:(length(breaks)-1), function(i){
    integrate(f = int.f, lower = breaks[i], upper = breaks[i+1], beta=beta, norm2=norm2, alpha=alpha)$value
  })
  probs <- probs/sum(probs)
  return(list(breaks=breaks,probs=probs))
}

int.f <- function(x, beta, norm2, alpha) {
  if (length(x) == 1) {
    return( exp(sum( poly(x = x, coefs = list(norm2=norm2,alpha=alpha), degree = length(beta)-1)*beta[2:length(beta)] + beta[1] )) )
  }
  exp(as.vector(poly(x = x, coefs = list(norm2=norm2,alpha=alpha), degree = length(beta)-1)%*%beta[2:length(beta)]) + beta[1])
}

Sim.Mu <- function(mu.beta, mu.noise, sigma.beta, sigma.noise, rho) {
  V <- matrix(c(sigma.beta^2,sigma.beta*sigma.noise*rho,sigma.beta*sigma.noise*rho,sigma.noise^2), nrow=2)
  return( c(mu.beta,mu.noise) + as.vector(SQRT.mat(V)%*%rnorm(2)) )
}

Get.True.Indices <- function(Sim.Results, min.match=3, p.forward=NULL, z.forward=NULL, mz.forward=NULL, Sim, ppm.precursor=10, weight.H0=1, rownames.pep=F, my.ind=1:5, prosit.ind=10, sim.score.ind=13, min.match.ind=7, n.ions.theory=6) {
  Headers <- Sim$Headers.sim
  if (rownames.pep) {
    out.true <- t(sapply(1:length(Sim.Results$Results), function(i){
      x <- Sim.Results$Results[[i]]$FM
      if (is.na(x[1])) {return(c(NA,Headers$precursorCharge[i],0))}
      tmp <- which(gsub(pattern="\\.[0-9]+$",replacement="",x=rownames(x))==Sim$Pep[i])
      return(c(ifelse(length(tmp)==0,NA,tmp),Headers$precursorCharge[i],length(rownames(x))))
    }))
  } else {
    out.true <- t(sapply(1:nrow(Headers), function(scan){
      ind <- which(abs(Headers$precursorMZ[scan]/mz.forward-1) <= ppm.precursor/10^6 & z.forward==Headers$precursorCharge[scan])
      tmp <- which(p.forward[ind] == Sim$Pep[scan])
      c(ifelse(length(tmp)==0,NA,tmp),Headers$precursorCharge[scan])
    }))
  }
  out.rest <- t(sapply(Sim.Results$Results,function(x){
    if (is.na(x$FM[1])) {return(rep(NA,6))}
    ind.use <- which(x$FM[,min.match.ind] >= min.match)
    if (length(ind.use)==0) {return(rep(NA,6))}
    if (weight.H0 > 0) {
      tmp <- rowSums(rbind(x$FM[ind.use,my.ind])) - log(weight.H0)
      out.mine <- which.max(tmp)
      out.prosit <- which.max(rbind(x$FM[ind.use,prosit.ind])); out.prosit <- ifelse(length(out.prosit)==0,NA,out.prosit)
      out.simscore <- which.max(rbind(x$FM[ind.use,sim.score.ind])); out.simscore <- ifelse(length(out.simscore)==0,NA,out.simscore)
      out.ionfrac <- which.max(x$FM[ind.use,min.match.ind]/x$FM[ind.use,n.ions.theory]); out.ionfrac <- ifelse(length(out.ionfrac)==0|is.nan(out.ionfrac),NA,out.ionfrac)
      return(c(1/(sum(exp(c(tmp,0)-max(tmp)))),max(tmp),ind.use[c(out.mine, out.prosit, out.simscore, out.ionfrac)]))
    } else {
      tmp <- rowSums(rbind(x$FM[ind.use,my.ind]))
      out.mine <- which.max(tmp)
      out.prosit <- which.max(rbind(x$FM[ind.use,prosit.ind])); out.prosit <- ifelse(length(out.prosit)==0,NA,out.prosit)
      out.simscore <- which.max(rbind(x$FM[ind.use,sim.score.ind])); out.simscore <- ifelse(length(out.simscore)==0,NA,out.simscore)
      out.ionfrac <- which.max(x$FM[ind.use,min.match.ind]/x$FM[ind.use,n.ions.theory]); out.ionfrac <- ifelse(length(out.ionfrac)==0|is.nan(out.ionfrac),NA,out.ionfrac)
      return(c(1/(sum(exp(c(tmp)-max(tmp)))),max(tmp),ind.use[c(out.mine, out.prosit, out.simscore, out.ionfrac)]))
    }
  }))
  if (ncol(out.true)==2) {
    out <- cbind(out.true[,1],out.rest,out.true[,2])
    colnames(out) <- c("True", "Pprob", "Score", "Mine", "PROSIT", "SimScore", "IonFrac", "Charge")
  } else {
    out <- cbind(out.true[,1],out.rest,out.true[,2:3])
    colnames(out) <- c("True", "Pprob", "Score", "Mine", "PROSIT", "SimScore", "IonFrac", "Charge", "Npeps")
  }
  return(out)
}

Get.Ion.Frac <- function(Search.Results, Results, ind.theory=6, ind.match=7) {
  i.max <- min(which(sapply(Search.Results,is.null))[1]-1,which(sapply(Results,is.null))[1]-1)
  if ("IonFrac"%in%colnames(Results[[1]]$out)) {return(Results)}
  Results <- lapply(1:i.max,function(i){
    copy <- Results[[i]]
    out <- copy$out
    ion.ind <- sapply(Search.Results[[i]]$Results,function(z){
      z <- z$FM
      if (is.na(z[1])) {return(NA)}
      tmp <- z[,ind.match]/z[,ind.theory]; ind.tmp <- which(!is.nan(tmp))
      if (length(ind.tmp)==0) {return(NA)}
      return(ind.tmp[which.max(tmp[ind.tmp])])
    })
    nn <- colnames(out)
    out <- cbind(out,ion.ind)
    colnames(out) <- c(nn,"IonFrac")
    copy$out <- out
    return(copy)
  })
  return(Results)
}

Calibrate.PProb <- function(out=NULL, Sim.Results=NULL, Sim.Results.DB2=NULL, p.forward=NULL, z.forward=NULL, mz.forward=NULL, Sim=NULL, bins=c(0,0.4,0.6,0.8,0.9,1), weight.H0=1, rownames.pep=F, min.score=-Inf, plot.it=T, cols=c("black","blue","violet"), cex=2, charge.bin=list(2,3,4:10), prob.min=0, correct.weight=F) {
  if (is.null(out)) {out <- Get.True.Indices(Sim.Results = Sim.Results, p.forward = p.forward, z.forward = z.forward, mz.forward = mz.forward, Sim = Sim, weight.H0 = weight.H0, rownames.pep = rownames.pep)}
  if (!is.null(Sim.Results.DB2)) {out2 <- Get.True.Indices(Sim.Results = Sim.Results, p.forward = p.forward, z.forward = z.forward, mz.forward = mz.forward, Sim = Sim, weight.H0 = weight.H0, rownames.pep = rownames.pep)}
  if (length(min.score)==1) {min.score <- rep(min.score, length(charge.bin))}
  mids <- bins[2:length(bins)]/2 + bins[1:(length(bins)-1)]/2
  out.true <- matrix(NA, nrow=length(charge.bin), ncol=length(mids))
  colnames(out.true) <- round(mids,digits=2); rownames(out.true) <- sapply(charge.bin,function(x){x[1]})
  out.theory <- out.true; out.N <- out.true
  if (correct.weight) {out[,2] <- 1/(1/out[,2] - exp(-out[,3]))}
  if (!is.null(out)) {col.charge <- ifelse(ncol(out)==7,7,ncol(out)-1)}
  for (c in 1:length(charge.bin)) {
    ind.c <- out[,col.charge]%in%charge.bin[[c]] & !is.na(out[,1]) & !is.na(out[,4]) & out[,3]>=min.score[c]
    tmp <- sapply(1:length(mids),function(i){
      tmp.ind <- ind.c & out[,2]>bins[i] & out[,2]<=bins[i+1]
      return(c(sum(tmp.ind),sum(tmp.ind&out[,1]==out[,4])/sum(tmp.ind),mean(out[tmp.ind,2])))
    })
    out.true[c,] <- tmp[2,]
    out.theory[c,] <- tmp[3,]
    out.N[c,] <- tmp[1,]
  }
  if (plot.it) {
    if (length(prob.min)==1) {prob.min <- rep(prob.min,length(charge.bin))}
    for (c in 1:length(charge.bin)) {
      #cex.point <- cex.max*sqrt(out.N[c,])/sqrt(max(out.N))
      plot(out.theory[c,], out.true[c,], col=cols[c], xlim=c(prob.min[c],1), ylim=c(prob.min[c],1), xlab="Predicted", ylab="Observed", pch=4, cex=cex)
      abline(a=0, b=1, col="red")
    }
  }
  return(list(true=out.true,theory=out.theory,N=out.N,out=out))
}

Results.Other <- function(Results, Sim.Results) {
  n.return <- 6
  return(lapply(1:length(Results),function(i){
    r.i <- Results[[i]]$out
    s.i <- Sim.Results[[i]]$Results
    out <- t(sapply(1:length(s.i),function(j){
      s.ij <- s.i[[j]]$FM
      if (is.na(s.ij[1])) {return(rep(NA,n.return))}
      x.prosit <- ifelse(is.na(r.i[j,5]),NA,s.ij[r.i[j,5],7])
      x.sim <- ifelse(is.na(r.i[j,6]),NA,s.ij[r.i[j,6],10])
      return(c(r.i[j,1],r.i[j,5],r.i[j,6],r.i[j,length(r.i[j,])],x.prosit,x.sim))
    }))
    colnames(out) <- c("True.ind", "Prosit.ind", "Sim.ind", "IonFrac.ind", "Prosit", "Sim")
    return(out)
  }))
}

Combine.Charge <- function(Results.other, Results.all) {
  return(lapply(1:length(Results.all),function(i){
    out <- cbind(Results.other[[i]][,1:5], Results.all[[i]]$out[,ncol(Results.all[[i]]$out)])
    colnames(out) <- c("True.ind", "Prosit.ind", "Sim.ind", "Prosit", "Sim", "Charge")
    return(out)
  }))
}

###Plot %incorrect###
BoxPlot.Incorrect <- function(Data, min.y=NULL, col.plot=c("black", "blue", "violet"), diff.axis=1.5, axis.labels=c("Ours","Prosit"), ind.use=NULL, ylab="% incorrect matches") {
  if (is.null(ind.use)) {ind.use <- 1:length(Data)}
  Data <- Data[ind.use]; axis.labels <- axis.labels[ind.use]
  if (nrow(Data[[1]])>ncol(Data[[1]])) {Data <- lapply(Data,function(x){t(x)})}
  tmp <- unlist(lapply(Data,function(x){c(x)})); max.y <- 100*max(tmp,na.rm=T)
  if (is.null(min.y)) {min.y <- max(100*min(tmp,na.rm=T)-5,0)}
  tmp <- rlist::list.rbind(Data)
  plot.list <- lapply(1:nrow(tmp), function(i){100*tmp[i,]})
  
  tmp <- 1:length(col.plot)
  n.tmp <- length(tmp)
  colors <- rep(col.plot,times=length(Data))
  at.plot <- rep(tmp, length(Data)) + rep(0:(length(Data)-1), each=length(tmp))*(length(tmp)+diff.axis)
  
  at.1 <- seq(mean(tmp),mean(tmp)+(length(Data)-1)*length(tmp),by=length(tmp))+(0:(length(Data)-1))*diff.axis
  boxplot(x=plot.list, at = at.plot, boxcol=colors, outcol=colors, whiskcol=colors, axes=F, medcol=colors, staplecol=colors, xlab="", ylab="", ylim=c(min.y,max.y))
  axis(side = 1, at = at.1, labels = axis.labels, las = 2 )
  axis(side = 2, at = pretty(seq(min.y,max.y,length=10)))
  title(ylab = ylab)
  return(0)
}

Simulate.mz <- function(mz.mean, x, theta, sigma, df=Inf) {
  n <- length(x)
  pi.large <- 1/(1+exp(-(cbind(rep(1,n),poly(x=x,degree=length(theta)-1,raw=T))%*%theta)))
  ind.large <- rbinom(n = n, size = 1, prob = pi.large)==1
  delta <- rep(0, n)
  if (sum(ind.large) > 0) {
    delta.tmp <- sigma[2]*Sim.Tdist(n = sum(ind.large), df = df)
    ind.redo <- abs(delta.tmp)>20
    while(sum(ind.redo)) {
      delta.tmp[ind.redo] <- sigma[2]*Sim.Tdist(n = sum(ind.redo), df = df)
      ind.redo <- abs(delta.tmp)>20
    }
    delta[ind.large] <- delta.tmp
  }
  if (sum(!ind.large) > 0) {
    delta.tmp <- sigma[1]*Sim.Tdist(n = sum(!ind.large), df = df)
    ind.redo <- abs(delta.tmp)>20
    while(sum(ind.redo)) {
      delta.tmp[ind.redo] <- sigma[1]*Sim.Tdist(n = sum(ind.redo), df = df)
      ind.redo <- abs(delta.tmp)>20
    }
    delta[!ind.large] <- delta.tmp
  }
  return(mz.mean*(1 + delta/10^6))
}

Sim.Tdist <- function(n, df=Inf) {
  if (is.infinite(df)) {return(rnorm(n))}
  return(rt(n = n, df = df)*sqrt((df-2)/df))
}

Combine.Pprob <- function(Sim.Results1, Sim.Results2=NULL, Sim, min.n=5) {
  out <- list()
  if (is.null(Sim.Results2)) {
    Sim.Results2 <- list()
    Sim.Results2$Results <- lapply(Sim.Results1$Results,function(x){
      if (is.na(x$FM[1])) {return(x)}
      tmp <- list()
      ind.use <- grep(pattern = "^\\.", x = rownames(x$FM), ignore.case = F)
      if (length(ind.use)==0) {tmp$FM <- NA; return(tmp)}
      tmp$FM <- rbind(x$FM[ind.use,])
      return(tmp)
    })
    Sim.Results1$Results <- lapply(Sim.Results1$Results,function(x){
      if (is.na(x$FM[1])) {return(x)}
      tmp <- list()
      ind.use <- which(!grepl(pattern = "^\\.", x = rownames(x$FM), ignore.case = F))
      if (length(ind.use)==0) {tmp$FM <- NA; return(tmp)}
      tmp$FM <- rbind(x$FM[ind.use,])
      return(tmp)
    })
  }
  out$out1 <- Get.True.Indices(Sim.Results = Sim.Results1, p.forward = NULL, z.forward = NULL, mz.forward = NULL, Sim = Sim, weight.H0 = 0, rownames.pep = T)
  out$out2 <- Get.True.Indices(Sim.Results = Sim.Results2, p.forward = NULL, z.forward = NULL, mz.forward = NULL, Sim = Sim, weight.H0 = 0, rownames.pep = T)
  ind.use <- which(!is.na(rowSums(out$out1)) & !is.na(rowSums(out$out2[,-1])) & out$out1[,8]>=min.n & out$out2[,8]>=min.n)
  cnames <- c("Pprob", "Score1", "Score2", "N1", "N2", "Incorrect", "MatchDB2", "Charge")
  out$out <- matrix(NA, nrow=nrow(out$out1), ncol=length(cnames))
  colnames(out$out) <- cnames
  out$out[ind.use,c(2:length(cnames))] <- cbind(out$out1[ind.use,3], out$out2[ind.use,3], out$out1[ind.use,8], out$out2[ind.use,8], as.numeric(out$out1[ind.use,1]!=out$out1[ind.use,4]|out$out2[ind.use,3]>=out$out1[ind.use,3]), as.numeric(out$out2[ind.use,3]>=out$out1[ind.use,3]), out$out1[ind.use,7])
  ind.1 <- intersect(ind.use,which(out$out1[,3]>out$out2[,3])); ind.2 <- intersect(ind.use,which(out$out1[,3]<=out$out2[,3]))
  out$out[ind.1,1] <- 1/( 1/out$out1[ind.1,2] + exp(out$out2[ind.1,3]-out$out1[ind.1,3])/out$out2[ind.1,2] )
  out$out[ind.2,1] <- 1/( 1/out$out2[ind.2,2] + exp(out$out1[ind.2,3]-out$out2[ind.2,3])/out$out1[ind.2,2] )
  return(out)
  
}

Combine.Results.DB <- function(Results.list, bins=c(0,seq(0.5,1,by=0.1)), min.n=5, charge.bin=list(2,3,4:10), rem.mid=F) {
  out <- list()
  if (!is.list(charge.bin)) {charge.bin <- list(charge.bin)}
  #out$N <- matrix(0,nrow=length(charge.bin),ncol=length(bins)-1); N.2 <- out$N
  out$Obs <- array(0,dim=c(length(charge.bin),length(bins)-1,length(Results.list)))
  out$N <- out$Obs
  out$Theory <- out$Obs
  for (c in 1:length(charge.bin)) {
    z <- charge.bin[[c]]
    for (j in 1:length(Results.list)) {
      tt.j <- Results.list[[j]]$out
      tmp.cj <- sapply(1:(length(bins)-1),function(i){
        if (rem.mid) {
          ind.cj <- which(tt.j[,8]%in%z & tt.j[,1]>bins[i]&tt.j[,1]<=bins[i+1] & tt.j[,4]>=min.n & tt.j[,5]>=min.n & (tt.j[,1]>0.75 | tt.j[,1]<0.25))
        } else {
          ind.cj <- which(tt.j[,8]%in%z & tt.j[,1]>bins[i]&tt.j[,1]<=bins[i+1] & tt.j[,4]>=min.n & tt.j[,5]>=min.n)
        }
        if (length(ind.cj)==0) {return(rep(0,4))}
        return(c( sum(tt.j[ind.cj,6]), sum(tt.j[ind.cj,7]), median(tt.j[ind.cj,5]/(tt.j[ind.cj,4]+tt.j[ind.cj,5])), mean(tt.j[ind.cj,1]) ))
      })
      #out$N[c,] <- out$N[c,] + tmp.cj[1,]; N.2[c,] <- N.2[c,] + tmp.cj[2,]
      out$N[c,,j] <- tmp.cj[1,]
      out$Obs[c,,j] <- tmp.cj[2,]/tmp.cj[1,]
      out$Theory[c,,j] <- tmp.cj[3,]
    }
  }
  #out$Obs <- N.2/out$N
  #colnames(out$Obs) <- (bins[2:length(bins)] + bins[1:(length(bins)-1)])/2
  return(out)
}
