Clean.prosit <- function(Prosit, ppm.filter=40, add=F, remove.0s=F, pep=NULL, z=NULL) {
  if (!is.null(pep) && !is.null(z)) {
    tmp <- Get.BY.Ions(peptide = pep, z = z)
    Prosit[,1] <- as.numeric(tmp)[match(gsub(pattern = "/.*$", replacement = "", x = rownames(Prosit)),names(tmp))]
  }
  if (remove.0s) {Prosit <- Prosit[Prosit[,2]>1e-8,]}
  Prosit <- Prosit[order(Prosit[,1]),]
  m.vec <- Prosit[,1]
  I.vec <- Prosit[,2]
  clean.rest <- F
  if (ncol(Prosit)>2) {Rest <- cbind(Prosit[,3:ncol(Prosit)]); clean.rest <- T}
  n <- length(m.vec)
  names.prosit <- rownames(Prosit)
  
  #Filter peaks#
  if (!add) {
    Delta.1 <- (m.vec[2:n]/m.vec[1:(n-1)] - 1)*10^6
    ind.1 <- which(Delta.1<=ppm.filter & I.vec[1:(n-1)]>=I.vec[2:n])   #Keep the first peak
    if (length(ind.1) > 0) {m.vec <- m.vec[-(ind.1+1)]; I.vec <- I.vec[-(ind.1+1)]; if (clean.rest){Rest<-cbind(Rest[-(ind.1+1),])}; names.prosit <- names.prosit[-(ind.1+1)]; n <- length(m.vec)}
    Delta.1 <- (1-m.vec[1:(n-1)]/m.vec[2:n])*10^6
    ind.1 <- which(Delta.1<=ppm.filter & I.vec[1:(n-1)]<=I.vec[2:n])   #Keep the second peak
    if (length(ind.1) > 0) {m.vec <- m.vec[-ind.1]; I.vec <- I.vec[-ind.1]; if (clean.rest){Rest<-cbind(Rest[-ind.1,])}; names.prosit <- names.prosit[-ind.1]; n <- length(m.vec)}
    out <- cbind(m.vec,I.vec)
    colnames(out) <- colnames(Prosit); rownames(out) <- names.prosit
    return(out)
  }
  Delta.1 <- (m.vec[2:n]/m.vec[1:(n-1)] - 1)*10^6
  ind.1 <- which(Delta.1<=ppm.filter & I.vec[1:(n-1)]>=I.vec[2:n])   #Keep the first peak
  if (length(ind.1) > 0) {
    m.vec <- m.vec[-(ind.1+1)]
    I.vec[ind.1] <- I.vec[ind.1+1] + I.vec[ind.1]; I.vec <- I.vec[-(ind.1+1)]
    if (clean.rest) {Rest <- cbind(Rest[-(ind.1+1),])}
    names.prosit <- names.prosit[-(ind.1+1)]; n <- length(m.vec)
  }
  Delta.1 <- (1-m.vec[1:(n-1)]/m.vec[2:n])*10^6
  ind.1 <- which(Delta.1<=ppm.filter & I.vec[1:(n-1)]<=I.vec[2:n])   #Keep the second peak
  if (length(ind.1) > 0) {
    m.vec <- m.vec[-ind.1]
    I.vec[ind.1+1] <- I.vec[ind.1+1] + I.vec[ind.1]; I.vec <- I.vec[-ind.1]
    if (clean.rest) {Rest <- cbind(Rest[-ind.1,])}
    names.prosit <- names.prosit[-ind.1]; n <- length(m.vec)
  }
  if (clean.rest) {
    out <- cbind(m.vec,I.vec/max(I.vec),Rest)
  } else {
    out <- cbind(m.vec,I.vec/max(I.vec))
  }
  colnames(out) <- colnames(Prosit); rownames(out) <- names.prosit
  return(out)
}