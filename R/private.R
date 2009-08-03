.prlossgain <- function(aber, sign, bound, nsam) {
    numga   <- (length(sign[sign == aber]))/nsam
    bound1  <- c(0, bound, 1)
    for (i in 1:(length(bound)+1)) {
        if (numga >= bound1[i] & numga <= bound1[i+1]) clas <- i
    }
    return(clas)
}

.freqlossgain <- function(aber, sign, nsam) {
    numga <- (length(sign[sign == aber]))/nsam
    return(numga)
}

.findbp <- function(reg, bppos) {
    st <- bppos[reg[1]]
    end <- bppos[reg[2]]
    c(st, end)
}

.findchr <- function(reg, chr) {
    chr[reg[1]]
}

.dist1 <- function(c1, c2) {
    d <- abs(c1-c2)
    apply(d, 1, sum)
}


.dist2 <- function(c1, c2) {
    c3 <- c2
    if (nrow(c1) >1) for (i in 1:(nrow(c1)-1)) c3<-rbind(c3, c2)
    d <- abs(c1-c3)
    apply(d, 1, sum)
}

#function insbr inserts breaks per region based on distance larger than c
#some problems with null-regions: therefore always insert c(0,0)
#Also determine mono-regions: distance to both neighbours > 3*c/4.

.insbr <- function(reg, c, ct) {
    datareg     <- ct[reg[1]:reg[2],-1]
    datareg2    <- rbind(datareg, datareg[nrow(datareg),])[-1,]
    totdist     <- CGHregions:::.dist1(datareg, datareg2)
    totdist2    <- append(totdist, 0, 0)[-length(totdist+1)]
    indc        <- ct[reg[1]:reg[2],1]
    totdisti    <- cbind(indc, totdist, totdist2)
    newbr       <- rbind(c(0,0,0), totdisti[totdist>c,])
    newbr1      <- as.vector(newbr[,1]) 
    mono        <- rbind(c(0,0,0), totdisti[totdist>floor(3*c/4)&totdist2>floor(3*c/4),])
    mono1       <- as.vector(mono[,1])[-1]
    mono2       <- c(mono1, mono1-1)
    sort(unique(c(newbr1,mono2)))
}

#jump size, this function computes the jump-sizes, it returns the 
#absolute clone index and its index within the region
.jump <- function(reg, ctdat) {
    datareg     <- CGHregions:::.td(reg, ctdat)
    ind         <- CGHregions:::.tdind(reg, ctdat)
    datareg2    <- rbind(datareg, datareg[nrow(datareg),])[-1,]
    totdist     <- CGHregions:::.dist1(datareg, datareg2)
    maxim       <- max(totdist)
    selectmax   <- cbind(ind, totdist)[totdist==maxim,1]
    selectmax
}



#concatenate regions 
.concat <- function(reg, c, ctdat, breakchr) {
    regionsst   <- c(1)
    regionsend  <- c()
    for (i in 1:(nrow(reg)-1)) {
        regionstart <- reg[i+1,1]
        c1 <- CGHregions:::.td(c(reg[i,2], reg[i,2]), ctdat)
        c2 <- CGHregions:::.td(c(reg[i+1,1], reg[i+1,1]), ctdat)
        if (CGHregions:::.dist1(c1, c2) > c | is.element(regionstart, breakchr)) {
            regionsst   <- c(regionsst, regionstart)
            regionsend  <- c(regionsend, reg[i,2])
        }
    }
    regionsend <- c(regionsend, reg[nrow(reg),2])
    cbind(regionsst, regionsend)
}

#take datarows from counts data
.td <- function(reg, ct) {
    ct[ct$ind<=reg[2] & ct$ind>= reg[1],,drop=FALSE][,-1]
}

.ntd <- function(reg, ct) {
    datareg <- CGHregions:::.td(reg, ct)
    return(ifelse(!is.null(dim(datareg)), nrow(datareg), 1))
}

.tdind <- function(reg, ct) {
    ct[ct$ind<=reg[2]&ct$ind>= reg[1],,drop=FALSE][,1]
}


#check basic condition, check it on first and last and 6 equally spaced clones maximally
.checkcond <- function(reg, c, ctdat) {
    ctr <- CGHregions:::.td(reg, ctdat)
    ncl <- nrow(ctr) 
    if (ncl > 10) {
        skip    <- floor(ncl/6)
        takerow <- c(seq(1, ncl, skip), ncl)
    }
    else {
        takerow <- 1:ncl
    }
    ctreg   <- ctr[takerow,,drop=FALSE]
    ctuni   <- unique(ctreg)
    luni    <- dim(ctuni)[1]
    i       <- 1
    cond    <- 0
    while (i < (luni-1) & cond == 0) {
        for (j in luni:(i+1)) {
            totdist <- CGHregions:::.dist1(ctuni[i,], ctuni[j,])      
            if (totdist>c) cond <- 1
        }
        i <- i+1 
    }
    cond
}

#function for unique confs within region
#.check <- function(reg,c) 
#         {
#         ncl <- reg[2]-reg[1] + 1 
#         if (ncl > 10)
#         {
#         skip <- floor(ncl/6)
#         takerow <- c(seq(1,ncl,skip),ncl)
#         }
#         else {takerow <- 1:ncl}
#         ctreg <- (counts[reg[1]:reg[2],-1])[takerow,]
#         ctuni <- unique(ctreg)
#ctuni
#}

#compute right-gradient
.rightgrad <- function(cl,ct,nro) #changed 28/10/08
       {
         #ct<-ctreg;nro<-nr;cl<-1;
         mincl <- nro-cl+1        
         if (mincl==1) rgrad <- 0 else
         {
         minim <- min(5,mincl)
         cllst <- ct[(cl+1):(cl+minim-1),,drop=FALSE]
         #distvec <- dist2(cllst,ct[cl,]) 
         ct0 <- ct[cl,]
         distvec2 <- apply(cllst,1,function(x) {dist(rbind(x,ct0),method="manhattan")})
         weightvec <- 1/(1:(minim-1))
         rgrad <- (distvec2 %*% weightvec)/sum(weightvec)       
         rgrad
         }
         }


.wh <- function(x, i) {
    which(x == i)
}

.gradients <- function(reg, index, ctdat) {
    ctreg       <- CGHregions:::.td(reg, ctdat)
    ctind       <- CGHregions:::.tdind(reg, ctdat)
    clones      <- sapply(index, CGHregions:::.wh, x=ctind)
    nr          <- nrow(ctreg)
    grads       <- sapply(clones, CGHregions:::.rightgrad, ct=ctreg, nro=nr)
    maxim       <- max(grads)
    selectmax   <- cbind(index, grads)[grads==maxim,1]
    selectmax
}

#this function is only applied if both the jump-size and the gradient 
#for two clones are equal. Currently not re-computed!!!
.dist2mid <-  function(reg, index, ctdat) {
    ncl     <- nrow(CGHregions:::.td(reg, ctdat))
    nmid    <- (ncl+1)/2
    ctreg   <- CGHregions:::.td(reg, ctdat)
    ctind   <- CGHregions:::.tdind(reg, ctdat)
    clones  <- sapply(index, CGHregions:::.wh, x=ctind)
    di2mi   <- abs((1:ncl) - nmid +0.25)
    cl      <- which.min(di2mi[clones])
    index[cl]        
}
             
#defines new breaks
.breek <- function(reg, ctdat) {
    jee <- CGHregions:::.jump(reg, ctdat)
    if (length(jee)==1) {
        br <- jee[1]
    } else {
        grads <- CGHregions:::.gradients(reg, jee, ctdat)   
        if (length(grads)==1) {
            br <- grads[1]
        } else {
            br <- CGHregions:::.dist2mid(reg, grads, ctdat)     
        }   
    }
    br
}   

#compute gradient
#.grad <- function(cl,ct,nro)
#       {
#         mincl <- nro-cl+1        
#         if (cl==1) 0
#         else 
#         {
#         minim <- min(5,cl)
#         cllst <- ct[(cl-minim+1):(cl-1),]
#         distvec <- CGHregions:::.dist2(cllst,ct[cl,])
#         weightvec <- (minim-1):1
#         lgrad <- (distvec %*% weightvec)/sum(weightvec)
#         if (mincl==1) rgrad <- 0 else
#         {
#         minim <- min(6,mincl)
#         cllst <- ct[(cl+1):(cl+minim-1),]
#         distvec <- CGHregions:::.dist2(cllst,ct[cl,]) 
#         weightvec <- 1:(minim-1)
#         rgrad <- (distvec %*% weightvec)/sum(weightvec)       
#         }
#         max(lgrad,rgrad)
#         }
#         }

#COMPUTE unique signatures within region and their frequency

.distmat <- function(uni) {
    el      <- nrow(uni)
    cmat    <- array(NA, c(el, el))
    for (i in (1:el)) {
        for (j in (1:el)) {
            ifelse(j<i,cmat[i,j]<-cmat[j,i],cmat[i,j] <- CGHregions:::.dist1(uni[i,], uni[j,]))
        }
    }
    cmat
}

.countrow <- function(unireg, datareg) {
    testje <- function(unireg, reg1) {
        min(reg1==unireg)
    }
    filt    <- apply(datareg, 1, testje, unireg=unireg)
    datarel <- cbind(filt,datareg)[filt==1,-1]
    ifelse(is.null(dim(datarel)), 1, nrow(datarel))
}

#this is the SLOW function!!!!!
#.whichsign <- function(reg,ctdat)
#{
#datareg <-CGHregions:::.td(reg,ctdat)
#uni <- unique(datareg)
#if (!is.null(dim(uni)))
#{
#afst <- CGHregions:::.distmat(uni)%*%apply(uni,1,CGHregions:::.countrow,datareg=datareg)
##afst <- distmat(uni)
##afst <- apply(uni,1,countrow,datareg=datareg)
#minim <- which.min(afst[,1])
#return(list(uni[minim,],afst[minim,1]))
##return(apply(uni,1,countrow,datareg=datareg))
#} 
#else return(list(uni,0))
#}

#this is the FAST function!
.whichsign2 <- function(reg, ctdat, levels) {
    datareg <-CGHregions:::.td(reg, ctdat)
    uni     <- unique(datareg)
    if (!is.null(dim(uni))) {
        cummat  <- apply(datareg, 2, CGHregions:::.countlevels, levels=levels)
        afst    <- apply(uni, 1, CGHregions:::.dm, cm=cummat, levels=levels)
        minim   <- which.min(afst)
        return(list(uni[minim,], afst[minim]))
    } 
    else return(list(uni, 0))
}


.countlevels <- function(datacol, levels) {
    sapply(levels, whl <- function(x,i) length(which(x==i)), x=datacol)
}

.dm <- function(uniseq, cm, levels) {
    ele <- length(levels)+1
    all <- rbind(cm, uniseq)
    concol <-function(colvec, el=ele, lev) {
        elem    <- colvec[el]
        dt      <- colvec[1:(el-1)]
        crossprod(dt,   abs(elem-lev))
    }
    sum(apply(all, 2, concol, lev=levels))
}

#.regionact <- function(reg) {
#    datareg <- countsnordel[reg[1]:reg[2],]
#    return(sum(apply(datareg, 1, sum))/nrow(datareg))
#}

.regionact2 <- function(reg, ctdat) {
    datareg <- CGHregions:::.td(reg, ctdat)
    countf <- function(dr) {
        length(dr[dr != 0])
    }
    return(sum(apply(datareg, 1, countf))/nrow(datareg))
}


#initilisation: inserts breaks at chromosome borders; 
#breakpoint is index of first clone in new region
.deterreg <- function(CGHdata, crit, ncolm, normstate, levels) {
#CGHdata <- CGHdataTry;crit<- critst
    chr     <- as.numeric(CGHdata[,1])
    numcl   <- length(chr)
    ind     <- 1:numcl
    chrsh   <- append(chr, 1, 0)
    chrsh   <- chrsh[-length(chrsh)]
    dif     <- chr-chrsh
    difind  <- cbind(ind, dif)
    counts  <- cbind(ind, CGHdata[,-(1:2)])
    counts  <- counts[,1:(ncolm+1)]
    nchr    <- length(unique(chr))

    if (nchr > 1) {
        if (nchr > 2) {
            breaks <- difind[dif != 0,][,1] 
        } else {
            breaks <- difind[dif != 0,][1]
        }
        regions     <- cbind(append(breaks, 1, 0),append(breaks-1, numcl, length(breaks)))
        allbreaks   <- as.vector(apply(regions, 1, .insbr, c=crit, ct=counts))
        newb        <- c()
        for (i in 1:length(allbreaks)) {
            tp      <- allbreaks[[i]]
            newb    <- c(newb, tp[tp != 0]+1)
        }
    }
    
    if (nchr <= 1) {
        breaks      <- c()
        regions     <- c(1, numcl)
        allbreaks   <- .insbr(regions, crit, counts)
        newb        <- c()
        for (i in 1:length(allbreaks)){
            tp      <- allbreaks[[i]]
            newb    <- c(newb, tp[tp != 0]+1)
        }
    }


    #merge with old breaks
    allb        <- sort(c(breaks, newb))
    regions     <- cbind(append(allb, 1, 0), append(allb-1, numcl, length(allb)))
    del         <- regions[(regions[,1]-regions[,2]) == 0,1]
    countsdel   <- if (length(del) > 0) counts[-del,] else counts
    if (nrow(regions) > 1) { #adapted 23/08
        regions2 <- regions[(regions[,1]-regions[,2])<0,,drop=FALSE]
        if (nrow(regions2) > 1) {
            regions3 <- .concat(regions2, crit, countsdel, breaks) 
        } else {
            regions3 <- regions2
        }
    } else {
        regions3 <- regions
    }
    
    allcond <- apply(regions3, 1, .checkcond, c=crit, ctdat=countsdel)

    if (max(allcond) == 0) {
        selreg <- regions3
    } else {#only apply gradient ruler if necessary
        #violreg     <- regions3[which.max(allcond),,drop=FALSE]
        #regions3    <- rbind(regions3, violreg)
        selreg      <- c()
        stop        <- 0

        #recursive application of gradient ruler
        while (stop==0) {
            allcond <- apply(regions3, 1, .checkcond, c=crit, ctdat=countsdel)
            #note '0' indicates that region satisfies criterion
            selreg0 <- cbind(regions3, allcond)
            selreg  <- rbind(selreg, selreg0[allcond==0,-3])
            newreg  <- selreg0[allcond==1,-3,drop=FALSE]

            if (!is.null(dim(newreg)) && dim(newreg)[1] != 0) {
                newbr   <- apply(newreg, 1, .breek, ctdat=countsdel)
                #insert newbr into newreg
                lbr     <- length(newbr)
                newreg2 <- c()
                for (i in (1:(lbr-1))) {
                    reg1 <- c(newreg[i,1], newbr[i])
                    reg2 <- c(newbr[i]+1, newreg[i,2])
                    newreg2 <- rbind(newreg2, c(newreg[i,1], newbr[i]), c(newbr[i]+1, newreg[i,2]))
                }
                regions3 <-newreg2
                #regions3 <- rbind(regions3, violreg)
            }
            if (is.null(dim(newreg)) || dim(newreg)[1] == 0) {
                stop <- 1
            }
        }
    }
    
    #DELETE MONO-REGIONS FROM SELREG
    
    selnew <- selreg[(selreg[,1]-selreg[,2])!=0,,drop=FALSE]
    
    #selects ACTIVE regions, assumes order loss, normal, gain.
    #seqnone <- seq(colnor,length(counts[1,-1]), by = nclass)
    #countsnordel <- counts[,-1][,-seqnone]


    sortedact       <- sort(apply(selnew, 1, .regionact2, ctdat=countsdel), index.return=TRUE, decreasing=TRUE)
    indicesActReg   <- sortedact$ix[1:ceiling(0.25*nrow(selnew))]
    ActReg25perc    <- selnew[indicesActReg,,drop=FALSE]

    if(is.null(nrow(ActReg25perc))) {
        ActReg25perc <- rbind(ActReg25perc, ActReg25perc)
    } # create two rows for rare case in which ActReg25perc consists of only 1 row

    #ASSUMES  missings are absent!!
    nsam    <- ncolm
    nclone  <- sum(apply(ActReg25perc, 1, function(y){y[2]-y[1]+1}))
    all     <- apply(ActReg25perc, 1, .whichsign2, ctdat=countsdel, levels=levels)
    avedist <- sum(as.vector(lapply(all, function(x) {x[[2]]}),mode="double"))/(nclone*nsam)
    return(list(avedist, selnew, countsdel))
}
