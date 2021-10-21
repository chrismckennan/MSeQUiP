library(mzR)
library(parallel)

#' Deisotope MS/MS spectra
#'
#' Deisotope MS/MS spectra using a modified chi squared statistic scoring function.
#'
#' @param mzml.file.path A character string containing the mzML or mzXML file path
#' @param n_cores Number of cores. Defaults to the total number of available cores minus 1.
#' @param ppm Fragment isotope mass error, in the ppm scale. Detaults to 20.
#' @param ppm.filter The software will collapse any peaks <=\code{ppm.filter}ppm apart into 1 peak. Defaults to 40.
#' @param remove.precursor Should the precursor be removed. Defaults to \code{F}.
#' @param write.mzml Should an mzML file be returned. Defaults to \code{F}.
#' 
#' @return If \code{write.mzml} is set to \code{T}, returns an mzML file in the same directory as \code{mzml.file.path}. Otherwise, it returns a list with the following entries: \item{Headers}{A data frame of headers for each MS and MS/MS spectrum.} \item{Spectra}{A list containing the MS and MS/MS spectra. The outputed MS/MS spectra contain 3 columns: m/z, intensity, and inferred charge.}
#' @export
#####Deisotope mzML file#####
Deisotope <- function(mzml.file.path, n_cores=NULL, ppm=20, ppm.filter=40, ppm.precursor=20, remove.precursor=F, max.addition=2, min.ions=20, min.chi=0.26, write.mzml=F) {
  mZ <- openMSfile(filename = mzml.file.path)
  Spectra <- peaks(mZ)
  Headers <- header(mZ)
  ind.msms <- which(Headers$msLevel==2)
  cl <- makeCluster(ifelse(is.null(n_cores),max(detectCores()-1,1),n_cores))
  clusterEvalQ(cl = cl, expr = {
    library(MSeQUiP)
  })
  clusterExport(cl = cl, varlist = c("ppm", "ppm.filter", "min.chi", "max.addition", "min.ions", "ppm.precursor", "remove.precursor"), envir = environment())
  out <- parLapply(cl = cl, X = lapply(ind.msms, function(i){return( list(spectrum=Spectra[[i]],header=Headers[i,]) )}),
                   function(x) {
                     if (NROW(x$spectrum) < min.ions) {
                       if (NCOL(x$spectrum) == 1) {return( rbind(c(x$spectrum,NA,NA)) )}
                       return(cbind(x$spectrum,matrix(NA,nrow=nrow(x$spectrum),ncol=2)))
                     }
                     tmp <- Deisotope.workhorse(Spectrum = x$spectrum, Prec.mz = x$header$precursorMZ, Iso.mz = x$header$isolationWindowTargetMZ, z = x$header$precursorCharge, PIW.left = x$header$isolationWindowLowerOffset, PIW.right = x$header$isolationWindowUpperOffset, ppm = ppm, ppm.filter = ppm.filter, min.chi = min.chi, max.addition = max.addition, remove.precursor = remove.precursor)
                     return(tmp[,1:3])
                   })
  stopCluster(cl = cl)
  Spectra[ind.msms] <- out
  if (write.mzml) {
    Spectra[ind.msms] <- lapply(Spectra[ind.msms], function(x){rbind(x[,1:2])})
    mzR::writeMSData(object = Spectra, file = gsub(pattern = "(\\.mzX*ML)$", replacement = "_deisotoped\\1", x = mzml.file.path, ignore.case = T, perl = T), header=Headers, outformat="mzml")
    return(0)
  } else {
    return(list(Spectra=Spectra, Headers=Headers))
  }
}

#####De-isotope MS/MS data#####
Deisotope.workhorse <- function(Spectrum, Prec.mz, Iso.mz, z=2, PIW.left=0.8, PIW.right=0.8, ppm=15, ppm.precursor=20, remove.precursor=T, ppm.filter=40, min.chi=0.26, Step=1.00284187613033, m.H=1.007825032, max.addition=1e2) {
  Prec.MH <- Prec.mz*z + (z-1)*m.H
  range.s <- Find.Range(mz=Prec.mz, z=z, ppm=ppm)
  Spectrum <- Spectrum[Spectrum[,1]>=range.s[1] & Spectrum[,1]<=range.s[2],]
  Spectrum <- Spectrum[order(Spectrum[,1]),]
  m.vec <- Spectrum[,1]; I.vec <- Spectrum[,2]
  n <- length(m.vec)
  peak.max <- min(floor((Iso.mz-Prec.mz+PIW.right)*z),4)
  peak.min <- max(floor((Iso.mz-Prec.mz-PIW.left)*z),0)
  if (peak.min>=peak.max) {return(cbind(Spectrum,matrix(NA,nrow=n,ncol=2)))}
  delta <- ppm/10^6
  
  #Filter peaks#
  Delta.1 <- (m.vec[2:n]/m.vec[1:(n-1)] - 1)*10^6
  ind.1 <- which(Delta.1<=ppm.filter & I.vec[1:(n-1)]>=I.vec[2:n])   #Keep the first peak
  if (length(ind.1) > 0) {m.vec <- m.vec[-(ind.1+1)]; I.vec <- I.vec[-(ind.1+1)]; n <- length(m.vec)}
  Delta.1 <- (1-m.vec[1:(n-1)]/m.vec[2:n])*10^6
  ind.1 <- which(Delta.1<=ppm.filter & I.vec[1:(n-1)]<=I.vec[2:n])   #Keep the second peak
  if (length(ind.1) > 0) {m.vec <- m.vec[-ind.1]; I.vec <- I.vec[-ind.1]; n <- length(m.vec)}
  
  if (remove.precursor) {
    tmp.out <- Prec.ions(S = cbind(m.vec,I.vec), Prec.mz = Prec.mz, z = z, ppm = ppm.precursor, peak.min = peak.min, peak.max = peak.max, min.chi = min.chi, Step = Step)
    m.vec <- tmp.out[,1]; I.vec <- tmp.out[,2]
    n <- length(m.vec)
  }
  
  #Filter peaks#
  ind.z <- t(sapply(1:(n-1), function(i) {
    m.i <- m.vec[i]
    ind.tmp <- (i+1):n
    return(sapply(1:z, function(charge) {
      if (m.i*charge+Step*(1-delta)>=range.s[2]) {return(NA)}
      ind.m <- ind.tmp[abs(m.vec[ind.tmp]/(m.i+Step/charge)-1)<=delta]
      return(ifelse(length(ind.m)>0,ind.m[which.max(I.vec[ind.m])],NA))
    }))
  }))
  ind.tmp <- which(apply(ind.z,1,function(x) {!all(is.na(x))}))
  if (length(ind.tmp) == 0) {return(cbind(m.vec,I.vec,matrix(NA,nrow=n,ncol=2)))}
  ind.z <- cbind(ind.tmp, rbind(ind.z[ind.tmp,]))
  
  #Test isotope peaks#
  ind.remove <- rep(NA, n); count.remove <- 1   #Indices to remove
  Iso.out <- matrix(NA, nrow=3, ncol=n); count.out <- 1  #3 x #isotope_peaks matrix giving monoisotopic peak index, total intensity and inferred charge
  Stats <- rep(NA, n)
  for (i in 1:nrow(ind.z)) {
    if (ind.z[i,1] %in% ind.remove) {next}
    ind.tmp <- (ind.z[i,1]+1):n
    ind.charge <- ind.z[i,2:(z+1)]  #length z vector
    charge.tmp <- which(!is.na(ind.charge))
    out.z <- lapply( charge.tmp, function(c) {
      if ( ind.z[i,c+1] %in% ind.remove ) {return(list(stat=Inf))}
      out <- list()
      ind.obs <- ind.z[i,c(1,c+1)]
      iso.theory <- Fragment.Profiles(M = m.vec[ind.obs[1]]*c - (c-1)*m.H, Prec = Prec.MH, peak.min = peak.min, peak.max = peak.max)
      iso.theory <- iso.theory[iso.theory > 1e-3]
      n.theory <- length(iso.theory)
      if (n.theory<=1) {return(list(stat=Inf))}
      
      count.tmp <- 2
      test <- ind.tmp[abs(m.vec[ind.tmp]/(m.vec[ind.z[i,1]]+Step*count.tmp/c)-1)<=delta]
      while( length(test) > 0 && length(ind.obs) <= n.theory+max.addition ) {
        ind.obs <- c(ind.obs,test[which.max(I.vec[test])])
        count.tmp <- count.tmp + 1
        test <- ind.tmp[abs(m.vec[ind.tmp]/(m.vec[ind.z[i,1]]+Step*count.tmp/c)-1)<=delta]
      }
      n.obs <- length(ind.obs)
      if (n.obs <= n.theory) {
        out$I.obs <- c(I.vec[ind.obs],rep(0,n.theory-n.obs))
        iso.obs <- out$I.obs/sum(out$I.obs)
        out$stat <- sum((iso.theory-iso.obs)^2/iso.theory)/(length(iso.theory)-1)
        out$ind.obs <- ind.obs
        return(out)
      }
      test <- sapply(1:(n.obs-n.theory+1),function(j) {
        iso.theory.tmp <- Fragment.Profiles(M = m.vec[ind.obs[j]]*c - (c-1)*m.H, Prec = Prec.MH, peak.min = peak.min, peak.max = peak.max)
        iso.theory.tmp <- iso.theory.tmp[iso.theory.tmp > 1e-3]; n.theory.tmp <- length(iso.theory.tmp)
        I.tmp <- I.vec[ind.obs[j:(n.theory.tmp+j-1)]]; I.tmp[is.na(I.tmp)] <- 0
        return( c(sum(I.tmp), sum((I.tmp/sum(I.tmp)-iso.theory.tmp)^2/iso.theory.tmp)/(n.theory.tmp-1), j+sum(I.tmp>1e-8)-1) )
      })
      j.return <- ifelse( min(test[2,])>min.chi, which.max(test[1,]), which(test[2,]<=min.chi)[which.max(test[1,test[2,]<=min.chi])] )
      out$ind.obs <- ind.obs[j.return:test[3,j.return]]
      out$I.obs <- I.vec[out$ind.obs]
      out$stat <- test[2,j.return]
      return(out)
    } )
    ind.possible <- which(sapply(out.z,function(x){x$stat}) <= min.chi)
    if (length(ind.possible) == 0) {next}
    if (length(ind.possible) == 1) {
      ind.obs <- out.z[[ind.possible]]$ind.obs
      ind.remove[count.remove:(count.remove-2+length(ind.obs))] <- ind.obs[2:length(ind.obs)]
      count.remove <- count.remove + length(ind.obs) - 1
      Iso.out[,count.out] <- c( ind.obs[1], sum(out.z[[ind.possible]]$I.obs), charge.tmp[ind.possible] )
      Stats[count.out] <- out.z[[ind.possible]]$stat
      count.out <- count.out + 1
    } else {
      ind.possible <- ind.possible[which.max( sapply(out.z[ind.possible],function(x){sum(x$I.obs)}) )]
      ind.obs <- out.z[[ind.possible]]$ind.obs
      ind.remove[count.remove:(count.remove-2+length(ind.obs))] <- ind.obs[2:length(ind.obs)]
      count.remove <- count.remove + length(ind.obs) - 1
      Iso.out[,count.out] <- c( ind.obs[1], sum(out.z[[ind.possible]]$I.obs), charge.tmp[ind.possible] )
      Stats[count.out] <- out.z[[ind.possible]]$stat
      count.out <- count.out + 1
    }
  }
  if (all(is.na(Iso.out[1,]))) {return( cbind(m.vec, I.vec, matrix(NA,nrow=n,ncol=2)) )}
  ind.remove <- ind.remove[!is.na(ind.remove)]
  Iso.out <- t(Iso.out[,!is.na(Iso.out[1,])])
  out <- cbind(m.vec, I.vec, matrix(NA,nrow=n,ncol=2)); out[Iso.out[,1],2] <- Iso.out[,2]
  out[Iso.out[,1],3] <- Iso.out[,3]; out[Iso.out[,1],4] <- Stats[!is.na(Stats)]
  out <- out[-ind.remove,]
  return( out[order(out[,1]),] )
}

##Estimate isotopic distribution given M+H and Prec+H##

Fragment.Profiles <- function(M, Prec, peak.min, peak.max) {   #Prec is M+H; M is M+H; +/-1 dalton does not matter; peak.max is 1-4
  d.max <- min(peak.max,4)+1
  d.min <- max(0,peak.min)+1
  P.M <- Mono.profiles(M)
  P.F <- Mono.profiles(Prec-M)
  tmp.F <- sapply(1:d.max, function(i){sum(P.F[max(d.min-i+1,1):(d.max+1-i)])})
  return( P.M[1:d.max]*tmp.F/sum(P.M[1:d.max]*tmp.F) )
}

Mono.profiles <- function(M, edit=T) {   #M is M+H
  out <- rep(0, 5)
  max.mass <- 3433.70753144348
  if (M<=58.02928875348) {out[1] <- 1; return(out)}
  if (M >= max.mass) {
    out <- c(0.132563514071985, 0.263898120083401, 0.242163515705798, 0.170695853663258, 0.0958189871160826)
    if (edit) {out <- out/sum(out)}
    return(out)
  }
  
  ##Monoisotopic##
  out[1] <- sum(c(1,poly(x = M, degree = 3, coefs = list(norm2=c(1,2528370,823074875441.7515869140625,434293780108884608,278793704858184690171904),alpha=c(869.6954795377945401924080215394496917724609375, 1376.013738615819193000788800418376922607421875, 1603.0337106112428955384530127048492431640625))))*c(0.6374016288729791313016903586685657501220703125, -284.0188014808364869168144650757312774658203125, 53.9315841840325589373605907894670963287353515625, -7.34704648996829945417630369774997234344482421875))
  
  ##Peak 1##
  out[2] <- sum(c(1,poly(x = M, degree = 6, coefs = list(norm2=c(1, 2528370, 823074875441.7515869140625, 434293780108884608, 278793704858184690171904, 179911221996490290932559642624, 113867962513018724434482124715720704, 79147930640191131495199470311773462069248),alpha=c(869.6954795377945401924080215394496917724609375, 1376.013738615819193000788800418376922607421875, 1603.0337106112428955384530127048492431640625, 1687.32364298276252156938426196575164794921875, 1642.3182210755066989804618060588836669921875, 1711.18289866309714852832257747650146484375))))*c(0.23956887373800242468924182048795046284794807434082, 119.1639373660456868719847989268600940704345703125, -65.964280192262549462611787021160125732421875, 15.7292604220928726732608993188478052616119384765625, -2.60291802087605450566343279206193983554840087890625, 0.2018888104888031109762636106097488664090633392334, 0.1430664647751463125491255823362735100090503692627))
  
  ##Peak 2##
  if (M <= 132.0483101212) {
    if (edit) {out <- out/sum(out)}
    return(out)
  }
  out[3] <- sum(c(1,poly(x = M, degree = 5, coefs = list(norm2=c(1, 2334125, 711067755352.1357421875, 352472492426501888, 212623713217233209524224, 127746566949076572034205483008, 77618306749099971285762527182454784),alpha=c(930.34736787362135146395303308963775634765625, 1440.23651739157548945513553917407989501953125, 1665.5682120117935482994653284549713134765625, 1739.842921892001186279230751097202301025390625, 1686.61193720223627678933553397655487060546875))))*c(0.09325157029458637503172013794028316624462604522705, 89.2142700354662139261563424952328205108642578125, -8.024732356579807657226410810835659503936767578125, -5.39156299246661863122653812752105295658111572265625, 2.5654275079934460990216393838636577129364013671875, -0.84936260885983716484304295590845867991447448730469))
  
  ##Peak 3##
  if (M <= 545.25461449429) {
    if (edit) {out <- out/sum(out)}
    return(out)
  }
  out[4] <- sum(c(1,poly(x = M, degree = 5, coefs = list(norm2=c(1, 1521165, 363196789164.4503173828125, 153076018967164576, 82003848818709267218432, 42331543472448438922368778240, 23150872463547553903818060574228480),alpha=c(1206.300587584608592806034721434116363525390625, 1660.987886540066028828732669353485107421875, 1819.68764118964872977812774479389190673828125, 1821.4168081170455479878000915050506591796875, 1800.581916116849924947018735110759735107421875))))*c(0.04161754702996537519377184821678383741527795791626, 35.0127226021805739719638950191438198089599609375, 3.81194469511035105568907965789549052715301513671875, -2.63010528613787508689370042702648788690567016601562, 1.11783675158196871279869810678064823150634765625, -0.48106472090092189386467680378700606524944305419922))
  
  ##Peak 4##
  if (M <= 1006.33987750452) {
    if (edit) {out <- out/sum(out)}
    return(out)
  }
  out[5] <- sum(c(1,poly(x = M, degree = 4, coefs = list(norm2=c(1, 606966, 94222014571.230560302734375, 27691160690034396, 9046997781152573423616, 3229142531314495864549933056),alpha=c(1680.617306065823186145280487835407257080078125, 2040.96787865556461838423274457454681396484375, 2152.4223626263510595890693366527557373046875, 2078.84789659548869167338125407695770263671875))))*c(0.0234880963131236970842241618129264679737389087677, 9.506401309100954932773674954660236835479736328125, 2.01545080654551389542916695063468068838119506835938, -0.41186931055645065180570441043528262525796890258789, 0.01802332823100914055425292303880269173532724380493))
  
  if (edit) {out <- out/sum(out)}
  return(out)
}

Find.Range <- function(mz, z, ppm=20, m.G=57.0214637233, m.H=1.007825032, m.O=15.9949146221, z.min.glycine=3) {
  return( c( (m.G/z.min.glycine+m.H)*(1-ppm/10^6), ((mz*z-(z-1)*m.H)-m.G-2*m.H-m.O)*(1+ppm/10^6)) )
}

Prec.ions <- function(S, Prec.mz, z, ppm=20, peak.min, peak.max, min.chi=0.26, Step=1.00284187613033, m.H=1.007825032, m.NH3=17.0265491015, m.H2O=18.0105646863) {
  delta <- ppm/10^6
  ind.include <- which(S[,1] >= (Prec.mz - m.H2O/z)*(1-delta) & S[,1] <= (Prec.mz + m.H2O/z + peak.max*Step/z)*(1+delta))
  if (length(ind.include) == 0) {return(S)}
  m.vec <- S[ind.include,1]; I.vec <- S[ind.include,2]
  MZ.search <- c(Prec.mz, Prec.mz-m.H2O/z, Prec.mz-m.NH3/z)
  ind.remove <- matrix(NA, ncol=peak.max-peak.min+1, nrow=length(MZ.search))
  for (i in 1:length(MZ.search)) {
    mz.i <- MZ.search[i]
    ind.i <- which(abs(m.vec/mz.i-1)<=delta)
    if (length(ind.i)==0) {next}
    iso.theory <- Precursor.Isotope(mz.i*z-z*m.H)[1:(peak.max+1)]
    ind.chi <- iso.theory > 1e-2
    iso.theory <- iso.theory[ind.chi]; iso.theory <- iso.theory/sum(iso.theory)
    n.theory <- length(iso.theory)
    if (n.theory <= 1) {next}
    count.tmp <- 1
    test <- which(abs(m.vec/(m.vec[ind.i]+Step*count.tmp/z)-1)<=delta)
    ind.obs <- c(ind.i)
    while( length(test) > 0 && length(ind.obs) < n.theory ) {
      ind.obs <- c(ind.obs,test[which.max(I.vec[test])])
      count.tmp <- count.tmp + 1
      test <- which(abs(m.vec/(m.vec[ind.i]+Step*count.tmp/z)-1)<=delta)
    }
    n.obs <- length(ind.obs)
    if (n.obs==1) {ind.remove[i,1] <- ind.obs}
    iso.obs <- c(I.vec[ind.obs],rep(0,n.theory-n.obs)); iso.obs <- iso.obs/sum(iso.obs)
    if (sum((iso.theory-iso.obs)^2/iso.theory)/(n.theory-1) <= min.chi) {
      ind.remove[i,1:n.obs] <- ind.obs
    } else {
      ind.remove[i,1] <- ind.obs[1]
    }
  }
  ind.remove <- c(ind.remove); ind.remove <- ind.remove[!is.na(ind.remove)]
  if (length(ind.remove) == 0) {return(S)}
  return( S[-ind.include[ind.remove],] )
}

######Estimate precursor monoisotopic profile######
Precursor.Isotope <- function(M) {
  out <- rep(0, 6)
  if (M <= 589.34352612789) {
    out[1] <- 0.720264861672686
    out[2] <- 0.211931377562435
    out[3] <- 0.0440513205968971
    out <- out/sum(out)
    return(out)
  }
  
  #0 sulfers
  ##Peak 0##
  out[1] <- ifelse(M <= 3625.6824190989, sum(c(1,poly(x = M, degree = 3, coefs = list(norm2=c(1, 100000, 32331455988.514492034912109375, 16192890739808146, 9344450970232491606016),alpha=c(1511.13287023050133939250372350215911865234375, 1892.430386485211101899039931595325469970703125, 2031.19165279261278556077741086483001708984375))))*c(0.4412487567029614865532494150102138519287109375, -40.6407069092517332364877802319824695587158203125, 8.4079722572383115419825116987340152263641357421875, -1.52310811644549670695880649873288348317146301269531)), 0.11684873769142)
  
  ##Peak 1##
  out[2] <- ifelse(M <= 3505.92431840724, sum(c(1,poly(x = M, degree = 6, coefs = list(norm2=c(1, 49100, 16655640222.5086994171142578125, 7557586901856920, 5238973844206462697472, 2224358339126750864321019904, 1147380301804887825829699187113984, 537493881711065640087465430758501384192),alpha=c(1386.461262834273611588287167251110076904296875, 1922.9359978192869675694964826107025146484375, 1942.7087743206384402583353221416473388671875, 2148.961167157563977525569498538970947265625, 1802.59318813256322755478322505950927734375, 2127.962271738490017014555633068084716796875))))*c(0.31423739750447043928005541602033190429210662841797, 4.75196906514023709178218268789350986480712890625, -5.792980444820482688328411313705146312713623046875, 1.77708924656061895142045159445842728018760681152344, -0.52646834176096124480181970284320414066314697265625, 0.30102681686665472282626865307975094765424728393555, -0.18390885524724409627594923222204670310020446777344)), 0.255030296445825)
  
  ##Peak 2##
  out[3] <- ifelse(M <= 3505.92431840724, sum(c(1,poly(x = M, degree = 6, coefs = list(norm2=c(1, 49100, 16655640222.5086994171142578125, 7557586901856920, 5238973844206462697472, 2224358339126750864321019904, 1147380301804887825829699187113984, 537493881711065640087465430758501384192),alpha=c(1386.461262834273611588287167251110076904296875, 1922.9359978192869675694964826107025146484375, 1942.7087743206384402583353221416473388671875, 2148.961167157563977525569498538970947265625, 1802.59318813256322755478322505950927734375, 2127.962271738490017014555633068084716796875))))*c(0.13576193824353033345886387905920855700969696044922, 12.352779839927354288420247030444443225860595703125, -2.0381444615478532256247490295208990573883056640625, -0.31332043139200721482140465923293959349393844604492, 0.06816758533656194773975300904567120596766471862793, 0.11603171442985597827668442505455459468066692352295, -0.08054499297153419612449454234592849388718605041504)), 0.25660248135987)

  ##Peak 3##
  if (M > 665.30205523921) {
    out[4] <- ifelse(M <= 3505.92431840724, sum(c(1,poly(x = M, degree = 4, coefs = list(norm2=c(1, 45049, 14135085426.8064746856689453125, 5442412002512283, 3470782085741346291712, 1172251672770180691984908288),alpha=c(1454.39105066049023662344552576541900634765625, 2005.4407777584556242800317704677581787109375, 2093.68815796585067801061086356639862060546875, 2167.0437557501127230352722108364105224609375))))*c(0.05082876266955050487128531244707119185477495193481, 7.449025424317216703684607637114822864532470703125, 0.37866263765482860437217027538281399756669998168945, -0.38307552908073244646303123772668186575174331665039, 0.02237824916905332625627345066732232226058840751648)), 0.174878725107753)
  }

  ##Peak 4##
  if (M > 1278.59208923994) {
    out[5] <- ifelse(M <= 3505.92431840724, sum(c(1,poly(x = M, degree = 4, coefs = list(norm2=c(1, 20882, 4760552653.87915897369384765625, 1318691001612939.75, 321068448378206617600, 88837074071699416917475328),alpha=c(1925.95404238942546726320870220661163330078125, 2213.412738688932222430594265460968017578125, 2640.27963355222482277895323932170867919921875, 2256.63641385679648010409437119960784912109375))))*c(0.02759230357948707496040796627312374766916036605835, 2.40144052101060667681053928390610963106155395507812, 0.26752525730208587262382025073748081922531127929688, -0.03742231525016473381306525425316067412495613098145, 0.00134532369591727264147851261810728829004801809788)), 0.0937215220581988)
  }

  ##Peak 5##
  if (M > 2122.09597260249) {
    out[6] <- ifelse(M <= 3505.92431840724, sum(c(1,poly(x = M, degree = 2, coefs = list(norm2=c(1, 7753, 560436671.416139125823974609375, 63431241722938.84375),alpha=c(2455.12919764773732822504825890064239501953125, 2925.579381191853826749138534069061279296875))))*c(0.01534045596756576125552218314851415925659239292145, 0.51072066226535040023293277045013383030891418457031, 0.04576330562136381174820343176179449073970317840576)), 0.0410984507445831)
  }
  out.s0 <- out/sum(out)
  
  #>0 sulfers#
  if (M < 635.33124736653) {
    return(out.s0)
  }
  bins <- c(600, 1000, 1250, 1500, 1750, 2000, 2100, 2150, 2175, 2200, 2300, 2400, 2600, 2800, 3000, 1e6)
  Probs <- c(0.3652041, 0.3736631, 0.5006995, 0.6992348, 0.8297292, 0.7018779, 0.6924167, 0.6302521, 0.5361446, 0.3150956, 0.3292810, 0.4603650, 0.6435157, 0.7746479, 0.5868925)
  
  ##Peak 0##
  out <- rep(0, 6)
  if (M > 589.34352612789) {
    out[1] <- ifelse (M <= 3625.6824190989, sum(c(1,poly(x = M, degree = 3, coefs = list(norm2=c(1, 100000, 32331455988.514492034912109375, 16192890739808146, 9344450970232491606016),alpha=c(1511.13287023050133939250372350215911865234375, 1892.430386485211101899039931595325469970703125, 2031.19165279261278556077741086483001708984375))))*c(0.4412487567029614865532494150102138519287109375, -40.6407069092517332364877802319824695587158203125, 8.4079722572383115419825116987340152263641357421875, -1.52310811644549670695880649873288348317146301269531)), 0.11684873769142)
  }
  
  ##Peak 1##
  if (M > 635.33124736653) {
    out[2] <- ifelse (M <= 3625.6824190989, sum(c(1,poly(x = M, degree = 6, coefs = list(norm2=c(1, 50900, 14176480253.5966167449951171875, 7438483872119910, 3299095246461914841088, 1702801182693151977299247104, 912412455639437311645069137149952, 511976745077614399020144071803584643072),alpha=c(1631.395658504661923871026374399662017822265625, 1949.805524416352227490278892219066619873046875, 2091.53329450518822341109625995159149169921875, 2069.56273829989504520199261605739593505859375, 2110.529500393020498449914157390594482421875, 2079.80333515572056057862937450408935546875))))*c(0.31465695697432272126548014057334512472152709960938, 1.75491217514560671197898500395240262150764465332031, -4.38082962554083277240124516538344323635101318359375, 1.13100404127716647195711630047298967838287353515625, -0.2872562055197562247421672054770169779658317565918, 0.10405218209644753746268008853803621605038642883301, -0.07190122947234960093521038970720837824046611785889)), 0.230773200191172)
  }
  
  ##Peak 2##
  if (M > 635.33124736653) {
    out[3] <- ifelse (M <= 3625.6824190989, sum(c(1,poly(x = M, degree = 6, coefs = list(norm2=c(1, 50900, 14176480253.5966167449951171875, 7438483872119910, 3299095246461914841088, 1702801182693151977299247104, 912412455639437311645069137149952, 511976745077614399020144071803584643072),alpha=c(1631.395658504661923871026374399662017822265625, 1949.805524416352227490278892219066619873046875, 2091.53329450518822341109625995159149169921875, 2069.56273829989504520199261605739593505859375, 2110.529500393020498449914157390594482421875, 2079.80333515572056057862937450408935546875))))*c(0.1730102925470380081929278048846754245460033416748, 9.119463178602746467049655620940029621124267578125, -2.06482265814179921648019444546662271022796630859375, -0.2027600500473546729640617058976204134523868560791, 0.17759389861743452021514144689717795699834823608398, -0.01053010246016774890953637822121891076676547527313, -0.09955997108886370516867003743755049072206020355225)), 0.234942431753195)
  }
  
  ##Peak 3##
  if (M > 635.33124736653) {
    out[4] <- ifelse (M <= 3625.6824190989, sum(c(1,poly(x = M, degree = 4, coefs = list(norm2=c(1, 50900, 14176480253.5966167449951171875, 7438483872119910, 3299095246461914841088, 1702801182693151977299247104),alpha=c(1631.395658504661923871026374399662017822265625, 1949.805524416352227490278892219066619873046875, 2091.53329450518822341109625995159149169921875, 2069.56273829989504520199261605739593505859375))))*c(0.07346959741014562339511684285753290168941020965576, 7.4198242496895883135721305734477937221527099609375, 0.03944639289699625878116151511676434893161058425903, -0.3076248862888153468020391301251947879791259765625, 0.10078248577518292283805578790634172037243843078613)), 0.190148593040845)
  }
  
  ##Peak 4##
  if (M > 865.33501259485) {
    out[5] <- ifelse (M <= 3625.6824190989, sum(c(1,poly(x = M, degree = 4, coefs = list(norm2=c(1, 42309, 8871770774.2035503387451171875, 3507382819918933, 1237472304918037004288, 540998779660240034233057280),alpha=c(1775.601365004076114928466267883777618408203125, 2206.77027101012981802341528236865997314453125, 2224.809020797727498575113713741302490234375, 2246.48290721955572735168971121311187744140625))))*c(0.03024084562367039888397535207786859245970845222473, 3.4086977466223373767206794582307338714599609375, 0.3908721560564931141179556561837671324610710144043, -0.07529162349381973706741177920775953680276870727539, 0.06380336242826320503063186606595991179347038269043)), 0.115156131101243)
  }
  
  ##Peak 5##
  if (M > 1364.59607239652) {
    out[6] <- ifelse (M <= 3625.6824190989, sum(c(1,poly(x = M, degree = 2, coefs = list(norm2=c(1, 14042, 2023600030.6966474056243896484375, 376789132951063.1875),alpha=c(2276.25640400846532429568469524383544921875, 2600.36759212138349539600312709808349609375))))*c(0.01790623140638929536305745671143085928633809089661, 0.88466485095815172723376917929272167384624481201172, 0.19004024122041190958398715338262263685464859008789)), 0.0565789202824795)
  }
  out.s1 <- out/sum(out)
  prob.s1 <- Probs[findInterval(x = M, vec = bins)]
  return( (1-prob.s1)*out.s0 + prob.s1*out.s1 )
}
