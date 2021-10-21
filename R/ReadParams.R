#####Read training files#####
Read.BYobs <- function(path.dir, append.file=NA) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  byObs.path <- gsub(pattern = "/$", replacement = "", x = path.dir)
  param.BYObs <- list()
  for (c in 1:3) {
    Obs.tmp <- readLines(con = paste0(byObs.path,"/","BYObs.",c,append.file,".txt"))
    param.BYObs[[c]] <- list()
    param.BYObs[[c]]$ad <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^alpha", x = Obs.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYObs[[c]]$ad) <- c("alpha", "delta")
    param.BYObs[[c]]$Omega <- matrix(as.numeric(strsplit(x = strsplit(x = grep(pattern = "^Omega", x = Obs.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]]), nrow=2)
  }
  return(param.BYObs)
}

Read.BYint <- function(path.dir, append.file=NA) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  byInt.path <- gsub(pattern = "/$", replacement = "", x = path.dir)
  param.BYInt <- list()
  for (c in 1:3) {
    Int.tmp <- readLines(con = paste0(byInt.path,"/","BYInt.",c,append.file,".txt"))
    param.BYInt[[c]] <- list()
    param.BYInt[[c]]$beta <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^beta", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$beta) <- c("beta", "se")
    param.BYInt[[c]]$mu <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^mu,", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$mu) <- c("mu", "se", "rho")
    param.BYInt[[c]]$mu.noise <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^mu\\.noise", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$mu.noise) <- c("mu.noise", "se.noise")
    param.BYInt[[c]]$var <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^d", x = Int.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    names(param.BYInt[[c]]$var) <- c("d", "phi")
  }
  return(param.BYInt)
}

Read.BYMA <- function(path.dir, append.file=NA) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  byMA.path <- gsub(pattern = "/$", replacement = "", x = path.dir)
  if (grepl(pattern = "\\.txt$", x = byMA.path, perl = T)) {
    MA.tmp <- readLines(con = byMA.path)
  } else {
    MA.tmp <- readLines(con = paste0(byMA.path,"/","BYMassAccuracy",append.file,".txt"))
  }
  param.MA <- list()
  param.MA$beta <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^beta", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.MA$sigma <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^sigma", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.MA$bound <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^min", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  return(param.MA)
}

Read.NoiseMA <- function(path.dir, append.file=NA) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  Noise.MA.path <- gsub(pattern = "/$", replacement = "", x = path.dir)
  MA.tmp <- readLines(con = paste0(Noise.MA.path,"/","NoiseMassAccuracy",append.file,".txt"))
  param.MA.noise <- list()
  param.MA.noise$beta <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^beta", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.MA.noise$sigma <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^sigma", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.MA.noise$bound <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^min", x = MA.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
  return(param.MA.noise)
}

Read.RT <- function(path.dir, append.file=NA) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  rt.path <- gsub(pattern = "/$", replacement = "", x = path.dir)
  if (grepl(pattern = "\\.txt$", x = rt.path, perl = T)) {
    RT.tmp <- readLines(con = rt.path)
  } else {
    RT.tmp <- readLines(con = paste0(rt.path,"/","RetentionTime",append.file,".txt"))
  }
  param.RT <- list()
  param.RT$beta <- as.numeric(strsplit(x = strsplit(x = RT.tmp[1], split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.RT$lambda <- as.numeric(strsplit(x = strsplit(x = RT.tmp[2], split = ":\\s*")[[1]][2], split = ",")[[1]])
  param.RT$lambda <- param.RT$lambda/sum(param.RT$lambda)
  param.RT$sigma <- sqrt(as.numeric(strsplit(x = strsplit(x = RT.tmp[3], split = ":\\s*")[[1]][2], split = ",")[[1]]))
  return(param.RT)
}

Read.NoiseInt <- function(path.dir, append.file=NA, include.k=T) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  NoiseInt.path <- gsub(pattern = "/$", replacement = "", x = path.dir)
  param.Noise.Int <- list()
  for (c in 1:3) {
    Noise.tmp <- readLines(con = paste0(NoiseInt.path,"/","NoiseInt.",c,append.file,".txt"))
    param.Noise.Int[[c]] <- list()
    param.Noise.Int[[c]]$beta <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^beta", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    param.Noise.Int[[c]]$norm2 <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^norm2", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    param.Noise.Int[[c]]$alpha <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^alpha", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    param.Noise.Int[[c]]$bound <- as.numeric(strsplit(x = strsplit(x = grep(pattern = "^min", x = Noise.tmp, ignore.case = F, value = T), split = ":\\s*")[[1]][2], split = ",")[[1]])
    param.Noise.Int[[c]]$k <- ifelse(include.k,integrate(f = Int.Noise.int, lower = param.Noise.Int[[c]]$bound[1], upper = param.Noise.Int[[c]]$bound[2], beta=param.Noise.Int[[c]]$beta, coefs=list(norm2=param.Noise.Int[[c]]$norm2,alpha=param.Noise.Int[[c]]$alpha))$value,NA)
  }
  return(param.Noise.Int)
}

Int.Noise.int <- function(x, beta, coefs) {
  X <- poly(x = x, degree = length(beta)-1, coefs = coefs)
  if (length(x)==1) {return(exp(sum(X*beta[2:length(beta)] + beta[1])))}
  return( exp(as.vector(X%*%beta[2:length(beta)]) + beta[1]) )
}

Read.NoiseObs <- function(path.dir, append.file=NA) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  path.spline <- ifelse(length(grep(pattern = "/$", x = path.dir))>0,path.dir,paste0(path.dir,"/"))
  fit.2 <- ReadSplineFile(file1 = paste0(path.spline,"Fit.2",append.file,".txt"))
  fit.3 <- ReadSplineFile(file1 = paste0(path.spline,"Fit.3.1",append.file,".txt"), file2 = paste0(path.spline,"Fit.3.2",append.file,".txt"))
  fit.4 <- ReadSplineFile(file1 = paste0(path.spline,"Fit.4.1",append.file,".txt"), file2 = paste0(path.spline,"Fit.4.2",append.file,".txt"))
  Spline.Fit <- list(fit.2$out, fit.3$out, fit.4$out)
  Spline.cutoff <- c(fit.2$cut, fit.3$cut, fit.4$cut)
  return(list(Spline.Fit=Spline.Fit,Spline.cutoff=Spline.cutoff))
}

Read.AbsentIons <- function(dir.path, append.file=NA, parametrization=c("normal","binom")) {
  if (is.na(append.file)) {append.file <- path.dir[2]}
  path.dir <- path.dir[1]
  parametrization <- match.arg(parametrization, c("normal","binom"))
  out <- vector(mode = "list", length = 3)
  for (c in 1:3) {
    tmp.c <- readLines(con = paste0(dir.path,"/","AbsIons.",c,append.file,".txt"))
    p <- as.numeric(strsplit(x = grep(pattern = "^p", x = tmp.c, ignore.case = F, value = T), split = ":\\s*")[[1]][2])
    N <- as.numeric(strsplit(x = grep(pattern = "^N", x = tmp.c, ignore.case = F, value = T), split = ":\\s*")[[1]][2])
    if (parametrization == "normal") {
      out[[c]] <- c(p*N,(1-p)*N) 
    } else {
      out[[c]] <- c(p,N)
    }
  }
  return(out)
}

#######Read spline data#######
ReadSplineFile <- function(file1, file2=NULL) {
  X1 <- readLines(file1)
  out.1 <- list()
  for (i in 1:length(X1)) {
    t.i <- strsplit(x = X1[i], split = ":\\s+", perl = T)[[1]]
    out.1[[t.i[1]]] <- as.numeric(strsplit(t.i[2],",")[[1]])
  }
  if (is.null(file2)) {return(list(out=out.1,cut=NA))}
  out.2 <- list()
  X2 <- readLines(file2)
  for (i in 1:length(X2)) {
    t.i <- strsplit(x = X2[i], split = ":\\s+", perl = T)[[1]]
    out.2[[t.i[1]]] <- as.numeric(strsplit(t.i[2],",")[[1]])
  }
  return(list(out=list(out.1,out.2),cut=out.2$Bound[1]))  #as.numeric(strsplit(strsplit(X2[length(X2)],":\\s+")[[1]][2],",")[[1]][1])
}
