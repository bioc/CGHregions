CGHregions <- function(input, averror=0.01) {

    if (class(input) == "cghCall") {
        CGHdata <- input
        kolnam  <- sampleNames(CGHdata)
        ncolm   <- ncol(CGHdata)
        chromo  <- chromosomes(CGHdata)
        bppos   <- bpstart(CGHdata)
        numbr   <- nrow(CGHdata)
        CGHdata <- data.frame(chromosomes(CGHdata), bpstart(CGHdata), calls(CGHdata))
    } else if (class(input) == "character") {
        CGHdata <- read.table(input, sep="\t", header=T)
        kolnam  <- colnames(CGHdata)[5:ncol(CGHdata)]
        ncolm   <- ncol(CGHdata)-4
        chromo  <- CGHdata[,2]
        bppos   <- CGHdata[,3]
        numbr   <- nrow(CGHdata)
        CGHdata <- CGHdata[,-4]
        CGHdata <- CGHdata[,-1]
    } else if (class(input) == "data.frame") {
        CGHdata <- input
        kolnam  <- colnames(CGHdata)[5:ncol(CGHdata)]
        ncolm   <- ncol(CGHdata)-4
        chromo  <- CGHdata[,2]
        bppos   <- CGHdata[,3]
        numbr   <- nrow(CGHdata)
        CGHdata <- CGHdata[,-4]
        CGHdata <- CGHdata[,-1]       
    }
    
    critst      <- max(1, floor(ncolm/10))
    stepsize    <- max(1, round(ncolm/20))
    
    normstate   <- 0
    levels      <- c(-1, 0, 1, 2)
    thresh      <- averror
    
    if(numbr<=10000) {
        startnew    <- floor(numbr/3)
        minim       <- min(400, startnew)
        CGHdataTry  <- rbind(CGHdata[1:minim,], CGHdata[startnew:(startnew+minim),], CGHdata[(2*startnew):(2*startnew+minim),])
    } else {
        startnew    <- floor(numbr/6)
        minim       <- min(400, startnew)
        CGHdataTry  <- rbind(CGHdata[1:minim,], CGHdata[startnew:(startnew+minim),], CGHdata[(2*startnew):(2*startnew+minim),],
                            CGHdata[(3*startnew):(3*startnew+minim),], CGHdata[(4*startnew):(4*startnew+minim),], CGHdata[(5*startnew):(5*startnew+minim),])
    }
    
    pmt         <- proc.time()
    stoploop    <- 0
    cspr        <- CGHregions:::.deterreg(CGHdata=CGHdataTry, critst, ncolm, normstate, levels)
    critsatpr   <- cspr[[1]]
    print(c(critst,cspr[[1]],nrow(cspr[[2]])))
    ifelse(critsatpr <= thresh, cr <- critst+stepsize, cr <- max(0,critst-stepsize))
    while (stoploop == 0) {
        cs      <- CGHregions:::.deterreg(CGHdata=CGHdataTry, cr, ncolm, normstate, levels)
        critsat <- cs[[1]]
        print(c(cr, cs[[1]], nrow(cs[[2]])))
        if (critsat>= thresh) {
            if (critsatpr <= thresh) {
                stoploop    <- 1
                critfound   <- cr-stepsize + floor(stepsize*(thresh-critsatpr)/(critsat-critsatpr))
            } else {
                cr <- max(0,cr-stepsize)
            }
        } else {
            if (critsatpr >= thresh) {
                stoploop    <- 1
                critfound   <- cr + floor(stepsize*(thresh-critsat)/(critsatpr-critsat))
            } else {
                if (critsat == critsatpr) {
                    stoploop <- 1
                    critfound <- cr
                } else {
                    cr <- cr + stepsize
                }
            }     
        }
        critsatpr <- critsat
    }
    critfound
    proc.time()-pmt
    print("Tuning on small data set finished...started with entire data set")    
    
    pmt         <- proc.time()
    stoploop    <- 0
    cspr        <- CGHregions:::.deterreg(CGHdata=CGHdata, critfound, ncolm, normstate, levels)
    critsatpr   <- cspr[[1]]
    print(c(critfound, cspr[[1]], nrow(cspr[[2]]), nrow(CGHdata)-nrow(cspr[[3]])))
    ifelse(critsatpr <= thresh, cr <- critfound+1, cr <- critfound-1)
    while (stoploop == 0) {
        cs      <- CGHregions:::.deterreg(CGHdata=CGHdata, cr, ncolm, normstate, levels)
        critsat <- cs[[1]]
        print(c(cr, cs[[1]], nrow(cs[[2]])))
        if (critsat >= thresh) {
            if (critsatpr <= thresh) {
                stoploop        <- 1
                critfound       <- cr-1
                regionsfound    <- cspr[[2]]
                countnomono     <- cspr[[3]]
            } else {
                cr          <- cr-1
                critsatpr   <- critsat
                cspr        <- cs
            }
        } else {
            if (critsatpr >= thresh) {
                stoploop        <- 1
                critfound       <- cr
                regionsfound    <- cs[[2]]
                countnomono     <- cs[[3]]
            }
            else {
                cr          <- cr+1
                critsatpr   <- critsat
                cspr        <- cs
            }     
        }
    }
    
    critfound
    proc.time()-pmt
    print(paste("c = ",critfound,", nr of regions: ", nrow(regionsfound), sep=""))
    print("Finished with entire data set.")
    res     <- apply(regionsfound, 1, CGHregions:::.whichsign2, ctdat=countnomono, levels=levels)
    prof    <- t(sapply(res, function(x) {as.vector(x[[1]], mode="numeric")}))
    nclone  <- apply(regionsfound,1,CGHregions:::.ntd,ct=countnomono)
    aved    <- signif(as.vector(lapply(res, function(x) {x[[2]]}), mode="numeric")/nclone, digits=3)
    bp      <- t(apply(regionsfound, 1, CGHregions:::.findbp, bppos = bppos))
    chrreg  <- apply(regionsfound, 1, CGHregions:::.findchr, chr = chromo)
    towrite <- cbind(regionsfound, bp, chrreg, nclone, aved, prof)
    rownames(towrite) <- c()
    od      <- order(towrite[,1])
    towrite <- towrite[od, -(1:2)]
    kolnamnew <- c("bp start", "bp end", "chromosome", "nclone", "Ave Dist", kolnam)
    colnames(towrite) <- kolnamnew
    
    annotation  <- data.frame(Chromosome=towrite[,3], Start=towrite[,1], End=towrite[,2], Nclone=towrite[,4], AveDist=towrite[,5])
    metadata    <- data.frame(  labelDescription=c("Chromosomal position",
                                                    "Basepair position start",
                                                    "Basepair position end",
                                                    "Number of clones in region",
                                                    "Average distance"),
                                row.names=c("Chromosome",
                                            "Start",
                                            "End",
                                            "Nclone",
                                            "AveDist")
                                )
    
    dimLabels   <- c("featureNames", "featureColumns")
    annotation  <- new("AnnotatedDataFrame", data=annotation, dimLabels=dimLabels, varMetadata=metadata)
    result      <- new("cghRegions", regions=as.matrix(towrite[,6:ncol(towrite)]), featureData=annotation)
    result   
}
