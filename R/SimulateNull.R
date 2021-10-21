#' Simulate data from the paper
#' 
#' Simulate MS/MS spectra subject to score-ordering and database incompleteness errors. The former contain simulated signal and real noise peaks, and the latter are real spectra whose precursors have been shifted +/-10m/z.
#' 
#' @param ProcessedDataPath Path to processed data.
#' @param n.null Number of spectra to simulate subject to database incompleteness errors.
#' @param t.df.by Log-signal intensity errors are drawn from a t distribution with \code{t.df.by} degrees of freedom.
#' @param list.train.paths A list of training parameter file paths: \code{NoiseInt.path: Path to noise intensity parameter file; \code{rt.path}: Path to retention time parameter file; \code{byInt.path}: Path to signal intensity parameter file; \code{byMA.path}: Path to signal mass accuracy parameter file; \code{byObs.path}: Path to signal peak presence parameter file; \code{append.file}: A character string that is appended to each paramter file
#' @return A list. \item{Spectra.sim}{A list of simulated spectra. It is assumed each MS/MS peak's charge is unknown in these simulations.} \item{Headers.sim}{Headers for the simulated spectra} \item{Pep}{Generating peptide} \item{Charge}{Precursor charge} \item{Scale}{The intensity scale for each simulated spectrum}
#' @export
SimulateDataPaper <- function(ProcessedDataPath, n.null, t.df.by=4, list.train.paths) {
  Sim.Pep <- SimulateData.new(ProcessedDataPath = ProcessedDataPath, NoiseInt.path = list.train.paths$NoiseInt.path, rt.path = list.train.paths$rt.path, byInt.path = list.train.paths$byInt.path, byMA.path = list.train.paths$byMA.path, byObs.path = list.train.paths$byObs.path, append.file = list.train.paths$append.file, df.by = t.df.by)
  tmp <- load(ProcessedDataPath); rm(list=tmp[!tmp%in%c("Sim.N0")])
  Sim.Null <- SimulateData.Null(mzml.path = "", missing.params = list(n.null=n.null,mz.shift=10,max.map=5), Sim.N0 = Sim.N0)
  out <- list()
  out$Spectra.sim <- c(Sim.Pep$Spectra.sim,Sim.Null$Spectra.sim)
  n.total <- length(out$Spectra.sim); n.pep <- length(Sim.Pep$Spectra.sim)
  out$Headers.sim <- rbind(Sim.Pep$Headers.sim,Sim.Null$Headers.sim)
  out$Pep <- rep(NA,n.total); out$Pep[1:n.pep] <- Sim.Pep$Pep
  out$Scale <- c(Sim.Pep$Scale,Sim.Null$Scale)
  out$Charge <- rep(NA,n.total); out$Charge[1:n.pep] <- Sim.Pep$Charge
  return(out)
}

SimulateData.Null <- function(mzml.path, missing.params=list(n.null=1600,mz.shift=10,max.map=5), other.path=NULL, Sim.N0=NULL, data.F=NULL, delta.prec.ppm=10, files.out=NULL) {
  if (is.character(data.F)) {
    data.F.copy <- data.F
    cat("Reading Sim.0 file...")
    tmp <- load(other.path)
    Headers.0 <- Sim.0$Headers.sim
    cat("done\n")
    rm(list=tmp); data.F <- data.F.copy; rm(data.F.copy)
    cat("Reading database file...")
    tmp <- load(data.F)
    data.F <- out.F; rm(list=tmp)
    cat("done\n")
    go.Sim0 <- T
  } else {
    go.Sim0 <- F
  }
  if (!is.null(data.F)) {database.PSM <- sapply(data.F,function(x){paste0(x$Pep,".",x$z)}); ind.keep <- !duplicated(database.PSM); data.F <- data.F[ind.keep]; rm(database.PSM)}
  if (grepl(pattern = "\\.gz$", x = mzml.path, ignore.case = F) && !is.null(file.spec)) {
    mzml.path.gz <- mzml.path
    tmp.path <- paste0("highph/",gsub(pattern = "\\.mzML$", replacement = "", x = gsub(pattern = "^.*/", replacement = "", x = file.spec), ignore.case = T),".mzML")
    mzml.path <- paste0("../data/",tmp.path)
    if (!file.exists(mzml.path)) {
      untar(tarfile = mzml.path.gz, files = tmp.path, list = F, exdir = "../data")
    }
  }
  
  n.null <- missing.params$n.null; mz.shift <- missing.params$mz.shift; max.map <- missing.params$max.map
  out <- list()
  if (is.null(Sim.N0)) {
    cat("Cleaning spectra...")
    mZ <- mzR::openMSfile(filename = mzml.path)
    Spectra <- peaks(mZ)
    Headers <- header(mZ); rm(mZ)
    if (go.Sim0) {
      Scans.consider <- Prune.spectra(mz.prec = Headers$precursorMZ, ind.exclude = (1:nrow(Headers))%in%Headers.0$scan, MSlevel = Headers$msLevel, charge.prec = Headers$precursorCharge, mz.theory = sapply(data.F, function(x){x$mz}), delta.ppm = delta.prec.ppm, mu.ppm = 0, charge.theory = sapply(data.F, function(x){x$z}))
    } else {
      Scans.consider <- Prune.spectra(mz.prec = Headers$precursorMZ, MSlevel = Headers$msLevel, charge.prec = Headers$precursorCharge, mz.theory = sapply(data.F, function(x){x$mz}), delta.ppm = delta.prec.ppm, mu.ppm = 0, charge.theory = sapply(data.F, function(x){x$z}))
    }
    tmp <- Get.Random.Scans(n = n.null*5, mz.prec = Headers$precursorMZ[Scans.consider], charge.prec = Headers$precursorCharge[Scans.consider], mz.shift = mz.shift, max.map = max.map)
    All.scans <- Scans.consider[tmp$ind.scan]
    Spectra[All.scans] <- lapply(All.scans, function(s){
      (Deisotope.workhorse(Spectrum = Spectra[[s]], Prec.mz = Headers$precursorMZ[s], Iso.mz = Headers$isolationWindowTargetMZ[s], z = Headers$precursorCharge[s], PIW.left = Headers$isolationWindowLowerOffset[s], PIW.right = Headers$isolationWindowUpperOffset[s], ppm = 15, ppm.filter = 30, max.addition = 2, remove.precursor = F))[,1:3]
    })
    out$Spectra <- Spectra; out$Headers <- Headers; if (is.null(out$Headers$scan)) {out$Headers$scan <- 1:nrow(Headers)}
    out$possible.mz <- tmp$mz; out$All.Scans <- All.scans
    rm(Spectra,Headers)
    cat("done\n")
  } else {
    out$Spectra <- Sim.N0$Spectra; out$Headers <- Sim.N0$Headers
    out$possible.mz <- Sim.N0$possible.mz; if (is.null(out$Headers$scan)) {out$Headers$scan <- 1:nrow(Headers)}
    out$All.Scans <- Sim.N0$All.Scans
  }
  cat("Randomly choosing spectra...")
  ind.spectra <- sort(sample(x = 1:length(out$All.Scans), size = n.null, replace = F))
  out$Headers.sim <- out$Headers[out$All.Scans[ind.spectra],]
  out$Spectra.sim <- out$Spectra[out$All.Scans[ind.spectra]]
  out$Scale <- sapply(out$Spectra.sim,function(x){quantile(x[,2],0.9)})
  for (i in 1:length(out$Spectra.sim)) {out$Spectra.sim[[i]][,2] <- out$Spectra.sim[[i]][,2]/out$Scale[i]}
  out$Headers.sim$precursorMZ <- sapply(ind.spectra,function(j){
    sample(x = out$possible.mz[[j]], size = 1, replace = F)
  })
  cat("done\n")
  if (is.null(data.F)) {return(out)}
  saveRDS(object = out, file = files.out[[1]])
  cat("Creating files...")
  mz.dataf <- c(Headers.0$precursorMZ, unlist(out$possible.mz))
  charge.dataf <- c(Headers.0$precursorCharge, unlist(lapply(1:length(out$possible.mz),function(i){
    rep(out$Headers$precursorCharge[out$All.Scans[i]],length(out$possible.mz[[i]]))
  })))
  rm(out)
  tmp.F <- Clean.Database.sim(Prune.database(data.F = data.F, all.mz = mz.dataf, all.charge = charge.dataf, Sim.N0 = NULL, delta.ppm = delta.prec.ppm))
  rm(data.F)
  saveRDS(object = tmp.F, file = files.out[[2]]); l.f <- length(tmp.F); rm(tmp.F)
  tmp.R <- Create.Decoy.Database(file.decoy.vec = files.out[[3]][1:2], file.Sim.0 = NULL, Sim.N0 = NULL, all.mz = mz.dataf, all.charge = charge.dataf, delta.prec.ppm = delta.prec.ppm)
  saveRDS(object = tmp.R, file = files.out[[3]][3]); l.r <- length(tmp.R); rm(tmp.R)
  cat("done\n")
  return(c(l.f,l.r))
}

Create.Decoy.Database <- function(file.decoy.vec, file.Sim.0=NULL, Sim.N0=NULL, all.mz=NULL, all.charge=NULL, delta.prec.ppm=10) {
  if (!is.null(Sim.N0)) {
    cat("Reading Sim.0 file...")
    tmp <- load(file.Sim.0)
    mz.0 <- Sim.0$Headers.sim$precursorMZ
    charge.0 <- Sim.0$Headers.sim$precursorCharge
    rm(list=tmp)
    cat("done\n")
    mz.N0 <- unlist(Sim.N0$possible.mz)
    charge.N0 <- unlist(lapply(1:length(Sim.N0$possible.mz),function(i){
      rep(Sim.N0$Headers$precursorCharge[Sim.N0$All.Scans[i]],length(Sim.N0$possible.mz[[i]]))
    }))
    rm(Sim.N0)
    all.mz <- c(mz.0,mz.N0)
    all.charge <- c(charge.0,charge.N0)
  }
  
  cat("Creating decoy database...")
  tmp <- load(file.decoy.vec[1])
  data.R <- out.R.0; rm(list=tmp)
  tmp <- load(file.decoy.vec[2])
  data.R <- c(data.R,out.R.1); rm(list=tmp)
  cat("done\n")
  
  cat("Cleaning database...")
  out <- Clean.Database.sim(Prune.database(data.F = data.R, Sim.N0 = NULL, all.mz = all.mz, all.charge = all.charge, delta.ppm = delta.prec.ppm))
  rm(data.R)
  cat("done\n")
  return(out)
}

Prune.data.which <- function(Data, mz.prec, charge.prec, delta=10) {
  delta <- delta/10^6
  mz.data <- sapply(Data,function(x){x$mz})
  z.data <- sapply(Data,function(x){x$z})
  out <- vector(mode = "list", length = 9)
  for (z in 2:10) {
    ind.z <- which(z.data==z)
    if (length(ind.z)==0) {next}
    tmp.z <- which(charge.prec==z)
    if (length(tmp.z)==0) {next}
    out[[z-1]] <- unique(unlist(lapply(tmp.z,function(i){
      tmp <- which(abs(mz.prec[i]/mz.data[ind.z]-1)<=delta)
      if (length(tmp)==0) {return(NULL)}
      return(ind.z[tmp])
    })))
  }
  return(Clean.Database.sim(data = Data[unique(unlist(out))]))
}

Prune.database.list <- function(Data, Sim.0=NULL, Sim.N0, naive=F) {
  if (is.null(Sim.0)) {
    mz.prec <- unlist(Sim.N0$possible.mz)
    charge.prec <- unlist(lapply(1:length(Sim.N0$possible.mz),function(i){
      rep(Sim.N0$Headers$precursorCharge[Sim.N0$All.Scans[i]],length(Sim.N0$possible.mz[[i]]))
    }))
  } else {
    mz.prec <- c(Sim.0$Headers.sim$precursorMZ, unlist(Sim.N0$possible.mz))
    charge.prec <- c(Sim.0$Headers.sim$precursorCharge, unlist(lapply(1:length(Sim.N0$possible.mz),function(i){
      rep(Sim.N0$Headers$precursorCharge[Sim.N0$All.Scans[i]],length(Sim.N0$possible.mz[[i]]))
    })))
  }

  if (!naive) {return(Clean.Database.sim(data = Prune.database(data.F = Data, all.mz = mz.prec, all.charge = charge.prec)))}
  mz.data <- sapply(Data,function(x){x$mz})
  z.data <- sapply(Data,function(x){x$z})
  out <- vector(mode = "list", length = 9)
  for (z in 2:10) {
    ind.z <- which(z.data==z)
    if (length(ind.z)==0) {next}
    tmp.z <- which(charge.prec==z)
    if (length(tmp.z)==0) {next}
    out[[z-1]] <- unique(unlist(lapply(tmp.z,function(i){
      tmp <- which(abs(mz.prec[i]/mz.data[ind.z]-1)<=10/10^6)
      if (length(tmp)==0) {return(NULL)}
      return(ind.z[tmp])
    })))
  }
  return(Clean.Database.sim(data = Data[unique(unlist(out))]))
}

Get.Random.Scans <- function(n, mz.prec, charge.prec, mz.shift=10, max.map=5) {
  mz.choose <- mz.prec
  charge.choose <- charge.prec
  ind.sample <- sort(sample(x = 1:length(mz.prec), size = n, replace = F))
  out <- vector(mode = "list", length = n)
  for (i in 1:n) {
    mz.i <- mz.prec[ind.sample[i]]
    charge.i <- charge.prec[ind.sample[i]]
    ind.charge <- which(charge.choose==charge.i)
    out.rand <- Random.null(mz = mz.i, mz.all = mz.choose[ind.charge], mz.delta = mz.shift, mz.add = 1, max.map = max.map)
    if (is.na(out.rand[1])) {out[[i]] <- NA; next}
    out[[i]] <- mz.choose[ind.charge[out.rand]]
    mz.choose <- mz.choose[-ind.charge[out.rand]]; charge.choose <- charge.choose[-ind.charge[out.rand]]
  }
  ind.use <- sapply(out,function(x){!is.na(x[1])})
  return(list(ind.scan=ind.sample[ind.use],mz=out[ind.use]))
}

Random.null <- function(mz, mz.all, mz.delta=10, mz.add=1, max.map=5) {
  ind.choose <- c( which(mz.all<=mz-(mz.delta-mz.add) & mz.all>=mz-(mz.delta+mz.add)), which(mz.all>=mz+(mz.delta-mz.add) & mz.all<=mz+(mz.delta+mz.add)) )
  if (length(ind.choose) > 0) {
    return(sample(x = ind.choose, size = min(max.map,length(ind.choose)), replace = F))
  }
  if (mz.delta>30) {return(NA)}
  return(Random.null(mz = mz, mz.all = mz.all, mz.delta = mz.delta+2*mz.add, mz.add = mz.add, max.map = max.map))
}

Prune.spectra <- function(mz.prec, charge.prec, mz.theory, charge.theory, MSlevel, ind.exclude=NULL, delta.ppm=10, mu.ppm=0) {
  out <- vector(mode = "list", length = 9)
  if (!is.null(ind.exclude)) {
    if (is.logical(ind.exclude[1])) {ind.exclude <- which(ind.exclude)}
    ind.use <- (1:length(mz.prec))[-ind.exclude]
    mz.prec <- mz.prec[ind.use]
    charge.prec <- charge.prec[ind.use]
    MSlevel <- MSlevel[ind.use]
  }
  ind.level2 <- MSlevel==2
  for (z in 2:10) {
    ind.z <- which(charge.theory==z)
    if (length(ind.z)==0) {next}
    tmp.z <- which(charge.prec==z & ind.level2)
    if (length(tmp.z)==0) {next}
    Map.z <- Prune.database.sub(mz.prec = mz.prec[tmp.z], mz.theory = sort(mz.theory[ind.z]), delta.ppm = delta.ppm, mu.ppm = mu.ppm)
    out[[z-1]] <- tmp.z[which(!is.na(Map.z[,1]))]
  }
  out <- unique(unlist(out))
  if (is.null(ind.exclude)) {return(out)}
  return(ind.use[out])
}

Prune.spectra.sub <- function(mz.prec, mz.theory, delta.ppm=10, mu.ppm=0) {
  delta <- delta.ppm/10^6
  mz.theory <- sort(mz.theory)
  vec <- c(-Inf,mz.theory,Inf)
  n <- length(vec) - 2
  Int <- findInterval(x = mz.prec, vec = vec)
  out <- t(sapply(1:length(Int),function(i){
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
  }))
  return(which(rowSums(is.na(out))==0))
}

Clean.Database.sim <- function(data) {
  data <- data[!duplicated(sapply(data,function(x){paste0(x$Pep,".",x$z)}))]
  tmp <- order(sapply(data,function(x){x$mz}))
  return(data[tmp])
}

Prune.database <- function(data.F, Sim.N0=NULL, all.mz=NULL, all.charge=NULL, delta.ppm=10) {
  if (!is.null(Sim.N0)) {
    all.mz <- unlist(Sim.N0$possible.mz)
    all.charge <- unlist(lapply(1:length(Sim.N0$possible.mz),function(i){
      rep(Sim.N0$Headers$precursorCharge[Sim.N0$All.Scans[i]],length(Sim.N0$possible.mz[[i]]))
    }))
  }
  z.F <- sapply(data.F, function(x){x$z})
  mz.F <- sapply(data.F, function(x){x$mz})
  tmp <- order(mz.F); mz.F <- mz.F[tmp]; z.F <- z.F[tmp]; data.F <- data.F[tmp]
  Ind.F.Theory <- sapply(1:10,function(zz){which(z.F==zz)})
  Index.mat <- matrix(NA, nrow=length(all.mz), ncol=2)
  out <- vector(mode = "list", length = length(Ind.F.Theory))
  for (zz in 1:length(Ind.F.Theory)) {
    Ind.obs.zz <- which(all.charge == zz)
    if (!(length(Ind.F.Theory[[zz]]) > 0 && length(Ind.obs.zz) > 0)) {next}
    out.zz <- Prune.database.sub(mz.prec = all.mz[Ind.obs.zz], mz.theory = mz.F[Ind.F.Theory[[zz]]], delta.ppm = delta.ppm, mu.ppm = 0)
    ind.use.zz <- which(!is.na(out.zz[,1])); if (length(ind.use.zz)==0) {next}
    out[[zz]] <- unique(unlist(lapply(ind.use.zz,function(j){Ind.F.Theory[[zz]][out.zz[j,1]:out.zz[j,2]]})))
  }
  return(data.F[sort(unique(unlist(out)))])
}

####Precursor matching####
Prune.database.sub <- function(mz.prec, mz.theory, delta.ppm=10, mu.ppm=0) {
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