cat("--------------------------------------------------\n")
cat("--- Breeding exercise version 0.2 (2022-10-30) ---\n")
cat("-----                                        -----\n")
cat("--- Questions/Issues? Contact cyang@sruc.ac.uk ---\n")
cat("--------------------------------------------------\n")

### function to calculate variance.
var2 <- function(x){
  
  return( sum((x-mean(x))^2)/length(x) )
  
}

### function to obtain the prediction equation.
gen.pred <- function(geno, pheno){
  
  # Cross Validation.
  CV <- sample(x=1:nrow(geno), size=nrow(geno), replace=FALSE)
  CV <- split(CV, sort(CV%%5))
  
  # rrBLUP.
  pred <- matrix(0, nrow=nrow(pheno), ncol=ncol(pheno))
  colnames(pred) <- colnames(pheno)
  rownames(pred) <- NULL
  mu <- rep(0, ncol(pheno))
  u <- matrix(0, nrow=ncol(geno), ncol=ncol(pheno))
  
  for(i in 1:ncol(pheno)){
    for(j in 1:5){
      temp <- mixed.solve(y=pheno[-CV[[j]], i], Z=geno[-CV[[j]], ])
      pred[CV[[j]], i] <- c(temp$beta) + c(geno[CV[[j]],]%*%temp$u)
      mu[i] <- mu[i] + c(temp$beta)
      u[,i] <- u[,i] + c(temp$u)
    }
  }
  mu <- mu/5
  u <- u/5
  
  # calculate prediction accuracy.
  PA <- data.frame(t(diag(cor(pheno[,1:2], pred[,1:2], use="pairwise.complete.obs"))))
  colnames(PA) <- c("size", "color")
  
  return(list(mu=mu, u=u, PA=PA))
  
}

### function to calculate the cost of breeding program.
cost.calc <- function(nf, nc, nF2, method, GS,
                      rHRT, rPYT, rAYT, rEYT,
                      sel.F2, sel.HRT, sel.PYT, sel.AYT,
                      sF2, sHRT, sPYT, sAYT){
  
  # cost.
  cost <- vector()
  
  # P to F1 (greenhouse = 100/month, grow = 1/plant, cross = 5/cross).
  cost <- c(cost, 1200 + 1*nf + 5*nc)
  
  # F2/DH to HRT.
  if(method=="SSD"){
    
    # F1 to F2 (greenhouse = 100/month, grow = 1/plant, self = 1/plant)
    temp.cost <- 1200 + 1*nc + 1*nc
    
    # F2 to F3 (greenhouse = 100/month, grow = 1/plant, self = 1/plant, PS = 5/plant, MAS = 10/plant).
    if(sel.F2=="PS"){
      temp.cost <- temp.cost + 1200 + 1*nc*nF2 + 1*nc*sF2 + 5*nc*nF2
    } else if(sel.F2=="MAS"){
      temp.cost <- temp.cost + 1200 + 1*nc*nF2 + 1*nc*sF2 + 10*nc*nF2
    } else if(sel.F2=="NS"){
      temp.cost <- temp.cost + 1200 + 1*nc*nF2 + 1*nc*sF2
    }
    
    # F3 to F4-F5-F6-HRT (greenhouse = 100/month, grow = 1/plant, self = 1/plant).
    temp.cost <- temp.cost + 4*300 + 4*1*nc*sF2 + 4*1*nc*sF2
    
  } else if(method=="DH"){
    
    # F1 to DH.
    temp.cost <- 40*nc*nF2
    
    # DH to HRT (greenhouse = 100/month, grow = 1/plant, self = 1/plant, PS = 5/plant, MAS = 10/plant).
    if(sel.F2=="PS"){
      temp.cost <- temp.cost + 1200 + 1*nc*nF2 + 1*nc*sF2 + 5*nc*nF2
    } else if(sel.F2=="MAS"){
      temp.cost <- temp.cost + 1200 + 1*nc*nF2 + 1*nc*sF2 + 10*nc*nF2
    } else if(sel.F2=="NS"){
      temp.cost <- temp.cost + 1200 + 1*nc*nF2 + 1*nc*sF2
    }
    
  }
  cost <- c(cost, temp.cost)
  
  # HRT to PYT (field = 40/plot, phenotyping = 10/plot, genotyping = 20/line).
  if(sel.HRT=="PS"){
    cost <- c(cost, 40*nc*sF2*rHRT + 10*nc*sF2*rHRT)
  } else if(sel.HRT=="GS"){
    cost <- c(cost, 40*nc*sF2*rHRT + 20*nc*sF2)
  } else if(sel.HRT=="NS"){
    cost <- c(cost, 40*nc*sF2*rHRT)
  }
  
  # PYT to AYT (field = 40/plot, phenotyping = 10/plot, genotyping = 20/line).
  if(sel.PYT=="PS"){
    cost <- c(cost, 40*nc*sHRT*rPYT + 10*nc*sHRT*rPYT)
  } else if(sel.PYT=="GS"){
    cost <- c(cost, 40*nc*sHRT*rPYT + 20*nc*sHRT)
  } else if(sel.PYT=="NS"){
    cost <- c(cost, 40*nc*sHRT*rPYT)
  }
  
  # AYT to EYT1 (field = 40/plot, phenotyping = 10/plot, genotyping = 20/line).
  if(sel.AYT=="PS"){
    cost <- c(cost, 40*nc*sPYT*rAYT + 10*nc*sPYT*rAYT)
  } else if(sel.AYT=="GS"){
    cost <- c(cost, 40*nc*sPYT*rAYT + 20*nc*sPYT)
  } else if(sel.AYT=="NS"){
    cost <- c(cost, 40*nc*sPYT*rAYT)
  }
  
  # EYT1 to EYT2 (field = 40/plot, phenotyping = 10/plot).
  cost <- c(cost, 40*sAYT*rEYT + 10*sAYT*rEYT)
  
  # EYT2 to variety (field = 40/plot, phenotyping = 10/plot).
  cost <- c(cost, 40*sAYT*rEYT + 10*sAYT*rEYT)
  
  # return the cost.
  names(cost) <- c("F1", "F2/DH", "HRT", "PYT", "AYT", "EYT1", "EYT2")
  return(cost)
  
}

### function to simulate and evaluate multiple breeding programs.
BP.eval <- function(va, vae, vres, ca, ew, BP, nsim=10){
  
  # fixed parameters.
  chr <- c("chr1", "chr2", "chr3", "chr4", "chr5") # 5 chromosomes.
  genlen <- c(2.0, 1.8, 1.6, 1.4, 1.2) # genetic lengths.
  qtl.pct <- 0.01 # percentage of markers as QTLs.
  snp.pct <- 0.10 # percentage of markers that go into SNP chip.
  
  nf <- 1000 # 1000 founders.
  nF1 <- 1 # 1 progeny per cross.
  nF3 <- 1 # 1 progeny per F2 self.
  nF4 <- 1 # 1 progeny per F3 self.
  nF5 <- 1 # 1 progeny per F4 self.
  nF6 <- 1 # 1 progeny per F5 self.
  nHRT <- 1 # 1 progeny per F6 self.

  rF0 <- 2 # Number of field trial for founders (for GS training population).
  rF1 <- 1/10 # number of field trial for F1 (proxy for higher residual variance).
  rF2 <- 1/10 # number of field trial for F2 (proxy for higher residual variance).
  rF3 <- 1/10 # number of field trial for F3 (proxy for higher residual variance).
  rF4 <- 1/10 # number of field trial for F4 (proxy for higher residual variance).
  rF5 <- 1/10 # number of field trial for F5 (proxy for higher residual variance).
  rF6 <- 1/10 # number of field trial for F6 (proxy for higher residual variance).
  rHRT <- 1/5 # number of field trial for HRT (proxy for smaller plot size).
  
  # checks.
  if(!(is.numeric(va) & length(va)==3 & all(va >= 0))) stop("va is a non-negative numeric vector of length 3.")
  if(!(is.numeric(vae) & length(vae)==3 & all(vae >= 0))) stop("vae is a non-negative numeric vector of length 3.")
  if(!(is.numeric(vres) & length(vres)==3 & all(vres >= 0))) stop("vres is a non-negative numeric vector of length 3.")
  if(!(is.numeric(ca) & length(ca)==1 & all(ca >= -1 & ca <= 1))) stop("ca is a [-1, 1] numeric vector of length 1.")
  #if(!(is.numeric(rF0) & length(rF0) & all(rF0 > 0))) stop("rF0 is a positive numeric vector of length 1.")
  if(!(is.numeric(ew) & length(ew)==2)) stop("ew is a numeric vector of length 2.")
  if(!(is.numeric(nsim) & length(nsim)==1 & all(nsim >= 1))) stop("nsim is a numeric vector of length 1 and nsim > 1.")
  if(!(is.data.frame(BP) & ncol(BP)==15)) stop("BP is a data.frame with 15 columns.")
  colnames(BP) <- c("program", "nc", "nF2", "method",
                    "sel.F2", "sel.HRT", "sel.PYT", "sel.AYT",
                    "sF2", "sHRT", "sPYT", "sAYT",
                    "rPYT", "rAYT", "rEYT")
  if(any(BP$nc < 100 | BP$nc > 10000)) stop("100 <= BP[,2] <= 10000.")
  if(any(BP$nF2 < 10 | BP$nF2 > 1000)) stop("10 <= BP[,3] <= 1000.")
  if(!all(BP$method%in%c("SSD", "DH"))) stop("BP[,4] = SSD/DH.")
  if(!all(BP$sel.F2%in%c("NS", "PS", "MAS"))) stop("BP[,5] = NS/PS/MAS.")
  if(!all(BP$sel.HRT%in%c("NS", "PS", "GS"))) stop("BP[,6] = NS/PS/GS.")
  if(!all(BP$sel.PYT%in%c("NS", "PS", "GS"))) stop("BP[,7] = NS/PS/GS.")
  if(!all(BP$sel.AYT%in%c("NS", "PS", "GS"))) stop("BP[,8] = NS/PS/GS.")
  if(any(BP$nF2 < BP$sF2 | BP$sF2 < BP$sHRT | BP$sHRT < BP$sPYT | BP$sPYT*BP$nc < BP$sAYT)) stop("BP[,3] >= BP[,9] >= BP[,10] >= BP[,11], BP[,2]*BP[,11] >= BP[,12].")
  if(any(BP$rPYT < 1 | BP$rAYT < 1 | BP$rEYT < 1)) stop("BP[,13], BP[,14], BP[,15] >= 1.")
  
  # check if GS is required.
  GS <- rep(FALSE, nrow(BP))
  for(i in 1:nrow(BP)){
    temp <- sum(c(BP$sel.HRT[i]=="GS", BP$sel.PYT[i]=="GS", BP$sel.AYT[i]=="GS"))
    if(temp==1) GS[i] <- TRUE else if(temp > 1) stop("Only one of HRT/PYT/AYT can have GS.")
  }
  
  # set up the founding population.
  founder <- runMacs(nInd=nf,
                     nChr=length(chr),
                     segSites=genlen*10000,
                     inbred=TRUE,
                     species="GENERIC",
                     ploidy=2L,
                     manualGenLen=genlen)
  
  # set up two correlated traits, and a third trait with penalty on first trait.
  SP <- SimParam$new(founder)
  SP$addTraitAG(nQtlPerChr=genlen*10000/(1/qtl.pct),
                mean=c(0,0),
                var=c(va[1], va[2]),
                varGxE=c(vae[1], vae[2]),
                varEnv=0,
                corA=matrix(c(1,ca,ca,1), nrow=2, ncol=2),
                name=c("size", "color"))
  
  # extract the first trait, QTL markers and maf.
  trait.size <- SP$traits[[1]]
  qtl <- pullQtlGeno(pop=founder, simParam=SP) - 1
  maf <- colSums(qtl+1)/(2*nrow(qtl))
  maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  
  # sample 10 QTL markers for third trait (and avoid maf < 0.10).
  idx3 <- which(maf >= 0.10)
  if(length(idx3) < 10){
    stop("Insufficient QTL markers for third trait.")
  } else {
    idx3 <- sort(sample(x=idx3, size=10, replace=FALSE))
  }
  
  # manually add a third trait (negatively correlated with first trait).
  trait.shape <- trait.size
  trait.shape@addEff <- -1*sign(trait.shape@addEff)
  trait.shape@addEff[-idx3] <- 0
  trait.shape@addEff <- trait.shape@addEff/sqrt(var2(c(qtl%*%trait.shape@addEff))/va[3])
  trait.shape@intercept <- 0 - mean(c(qtl%*%trait.shape@addEff))
  trait.shape@gxeEff <- rnorm(n=length(trait.shape@gxeEff), mean=0, sd=1)
  trait.shape@gxeEff[-idx3] <- 0
  trait.shape@gxeEff <- trait.shape@gxeEff/sqrt(var2(c(qtl%*%trait.shape@gxeEff))/vae[3])
  trait.shape@intercept <- 0 - mean(c(qtl%*%trait.shape@gxeEff))
  trait.shape@name <- "shape"
  SP$manAddTrait(lociMap=trait.shape)
  
  # set up residual effects.
  SP$setVarE(varE=vres)
  
  # set up marker genotype array.
  SP$addSnpChip(nSnpPerChr=genlen*10000/(1/snp.pct), name="BestEverSNPChip")
  
  # create founders.
  founder <- newPop(founder, simParam=SP)
  founder <- setPheno(pop=founder, reps=rF0, simParam=SP)
  
  # progress update.
  cat("Simulation for founders and traits: done.\n")
  
  # calculate genetic and phenotypic covariances of the 3 traits in the founders.
  covG <- varG(founder)
  covP <- varP(founder)

  # calculate Smith-Hazel index from Economic Weights (just approx, CovG/CovP change over time!).
  sh <- rbind(smithHazel(ew, covG[1:2,1:2], covP[1:2,1:2]), 0)
  rownames(sh)[3] <- "shape"
  
  # collect results.
  out <- replicate(n=nrow(BP), list(list(PYT=vector(), AYT=vector(), EYT=vector(), cost=vector())))
  names(out) <- paste("BP_", BP$program, sep="")
  
  # get the prediction model fitted with the founders.
  if(any(GS)){
    cat("GS training ... ")
    temp <- gen.pred(geno=pullSnpGeno(pop=founder, snpChip=1, simParam=SP) - 1, pheno=pheno(founder))
    mu <- temp$mu
    u <- temp$u
    out <- c(out, PA=list(temp$PA))
    cat("done.\n")
  } else {
    mu <- NULL
    u <- NULL
    out <- c(out, PA=NA)
  }
  
  # loop through breeding program-i and nsim-j.
  for(i in 1:nrow(BP)){
    
    # progress update.
    cat(paste("Breeding Program ", i, ": sim ", sep=""))
    
    for(j in 1:nsim){
      
      # F1.
      F1 <- randCross(pop=founder,
                      nCrosses=BP$nc[i],
                      nProgeny=nF1,
                      simParam=SP)
      F1 <- setPheno(pop=F1, reps=rF1, simParam=SP)
      
      # F2/DH to HRT.
      if(BP$method[i]=="SSD"){
        
        # F2.
        F2 <- self(pop=F1, nProgeny=BP$nF2[i], simParam=SP)
        F2 <- setPheno(pop=F2, reps=rF2, simParam=SP)
        
        # F3.
        if(BP$sel.F2[i]=="PS"){
          
          # Phenotypic Selection.
          F3 <- selectWithinFam(pop=F2, nInd=BP$sF2[i], trait=3, use="pheno", simParam=SP)
          
        } else if(BP$sel.F2[i]=="MAS"){
          
          # Marker Assisted Selection.
          F3 <- selectWithinFam(pop=F2, nInd=BP$sF2[i], trait=3, use="gv", simParam=SP)
          
        } else if(BP$sel.F2[i]=="NS"){
          
          # Random selection.
          F3 <- selectWithinFam(pop=F2, nInd=BP$sF2[i], trait=3, use="rand", simParam=SP)
          
        }
        F3 <- self(pop=F3, nProgeny=nF3, simParam=SP)
        F3 <- setPheno(pop=F3, reps=rF3, simParam=SP)
        
        # F4 (SSD).
        F4 <- self(pop=F3, nProgeny=nF4, simParam=SP)
        F4 <- setPheno(pop=F4, reps=rF4, simParam=SP)
        
        # F5 (SSD).
        F5 <- self(pop=F4, nProgeny=nF5, simParam=SP)
        F5 <- setPheno(pop=F5, reps=rF5, simParam=SP)
        
        # F6 (SSD).
        F6 <- self(pop=F5, nProgeny=nF6, simParam=SP)
        F6 <- setPheno(pop=F6, reps=rF6, simParam=SP)
        
        # HRT (SSD).
        HRT <- self(pop=F6, nProgeny=nHRT, simParam=SP)
        HRT <- setPheno(pop=HRT, reps=rHRT, simParam=SP)
        
      } else if(BP$method[i]=="DH"){
        
        # DH.
        DH <- makeDH(pop=F1, nDH=BP$nF2[i], simParam=SP)
        DH <- setPheno(pop=DH, reps=rF2, simParam=SP)
        
        # HRT (select from DH to go into HRT).
        if(BP$sel.F2[i]=="PS"){
          
          # Phenotypic Selection.
          HRT <- selectWithinFam(pop=DH, nInd=BP$sF2[i], trait=3, use="pheno", simParam=SP)
          
        } else if(BP$sel.F2[i]=="MAS"){
          
          # Marker Assisted Selection.
          HRT <- selectWithinFam(pop=DH, nInd=BP$sF2[i], trait=3, use="gv", simParam=SP)
          
        } else if(BP$sel.F2[i]=="NS"){
          
          # Random selection.
          HRT <- selectWithinFam(pop=DH, nInd=BP$sF2[i], trait=3, use="rand", simParam=SP)
          
        }
        
      }
      
      # PYT (select from HRT).
      if(BP$sel.HRT[i]=="PS"){
        
        # Phenotypic Selection.
        PYT <- selectWithinFam(pop=HRT, nInd=BP$sHRT[i], trait=selIndex, b=sh, scale=TRUE, use="pheno", simParam=SP)
        
      } else if(BP$sel.HRT[i]=="GS"){
        
        # Genomic Selection.
        geno.HRT <- pullSnpGeno(pop=HRT, snpChip=1, simParam=SP) - 1
        HRT@ebv <- matrix(mu, nrow=nrow(geno.HRT), ncol=ncol(u), byrow=TRUE) + geno.HRT%*%u
        rownames(HRT@ebv) <- NULL
        colnames(HRT@ebv) <- c("size", "color", "shape")
        PYT <- selectWithinFam(pop=HRT, nInd=BP$sHRT[i], trait=selIndex, b=sh, scale=TRUE, use="ebv", simParam=SP)
        
      } else if(BP$sel.HRT[i]=="NS"){
        
        # Random selection.
        PYT <- selectWithinFam(pop=HRT, nInd=BP$sHRT[i], trait=1, use="rand", simParam=SP)
        
      }
      PYT <- setPheno(pop=PYT, reps=BP$rPYT[i], simParam=SP)
      
      # AYT (select from PYT).
      if(BP$sel.PYT[i]=="PS"){
        
        # Phenotypic Selection.
        AYT <- selectWithinFam(pop=PYT, nInd=BP$sPYT[i], trait=selIndex, b=sh, scale=TRUE, use="pheno", simParam=SP)
        
      } else if(BP$sel.PYT[i]=="GS"){
        
        # Genomic Selection.
        geno.PYT <- pullSnpGeno(pop=PYT, snpChip=1, simParam=SP) - 1
        PYT@ebv <- matrix(mu, nrow=nrow(geno.PYT), ncol=ncol(u), byrow=TRUE) + geno.PYT%*%u
        rownames(PYT@ebv) <- NULL
        colnames(PYT@ebv) <- c("size", "color", "shape")
        AYT <- selectWithinFam(pop=PYT, nInd=BP$sPYT[i], trait=selIndex, b=sh, scale=TRUE, use="ebv", simParam=SP)
        
      } else if(BP$sel.PYT[i]=="NS"){
        
        # Random selection.
        AYT <- selectWithinFam(pop=PYT, nInd=BP$sPYT[i], trait=1, use="rand", simParam=SP)
        
      }
      AYT <- setPheno(pop=AYT, reps=BP$rAYT[i], simParam=SP)
      
      # EYT (select from AYT).
      if(BP$sel.AYT[i]=="PS"){
        
        # Phenotypic Selection.
        EYT <- selectInd(pop=AYT, nInd=BP$sAYT[i], trait=selIndex, b=sh, scale=TRUE, use="pheno", simParam=SP)
        
      } else if(BP$sel.AYT[i]=="GS"){
        
        # Genomic Selection.
        geno.AYT <- pullSnpGeno(pop=AYT, snpChip=1, simParam=SP) - 1
        AYT@ebv <- matrix(mu, nrow=nrow(geno.AYT), ncol=ncol(u), byrow=TRUE) + geno.AYT%*%u
        rownames(AYT@ebv) <- NULL
        colnames(AYT@ebv) <- c("size", "color", "shape")
        EYT <- selectInd(pop=AYT, nInd=BP$sAYT[i], trait=selIndex, b=sh, scale=TRUE, use="ebv", simParam=SP)
        
      } else if(BP$sel.AYT[i]=="NS"){
        
        # Random selection.
        EYT <- selectInd(pop=AYT, nInd=BP$sAYT[i], trait=1, use="rand", simParam=SP)
        
      }
      EYT <- setPheno(pop=EYT, reps=2*BP$rEYT[i], simParam=SP)
      
      # collect the results.
      out[[i]]$PYT <- rbind(out[[i]]$PYT, data.frame(sim=j, pheno(PYT)))
      out[[i]]$AYT <- rbind(out[[i]]$AYT, data.frame(sim=j, pheno(AYT)))
      out[[i]]$EYT <- rbind(out[[i]]$EYT, data.frame(sim=j, pheno(EYT)))
      
      # progress update.
      cat(paste(j, " ", sep=""))
      
    }
    
    # calculate the cost.
    out[[i]]$cost <- cost.calc(nf=nf, nc=BP$nc[i], nF2=BP$nF2[i], method=BP$method[i], GS=GS[i],
                               rHRT=rHRT, rPYT=BP$rPYT[i], rAYT=BP$rAYT[i], rEYT=BP$rEYT[i],
                               sel.F2=BP$sel.F2[i], sel.HRT=BP$sel.HRT[i], sel.PYT=BP$sel.PYT[i], sel.AYT=BP$sel.AYT[i],
                               sF2=BP$sF2[i], sHRT=BP$sHRT[i], sPYT=BP$sPYT[i], sAYT=BP$sAYT[i])
    
    # progress update.
    cat("done.\n")
    
  }
  
  # return the results.
  return(out)
  
}

### function to plot the results.
BP.plot <- function(out=out, threshold=c(3, 3, 3), nsim=10){
  
  # identify the number of breeding programs.
  nbp <- length(out) - 1
  
  # get the stages.
  stage <- c("PYT", "AYT", "EYT")
  
  # get the traits.
  trait <- c("size", "color", "shape")
  
  # loop through the breeding program and traits.
  df <- vector()
  for(i in 1:nbp){
    for(j in stage){
      for(k in trait){
        
        df <- rbind(df,
                    data.frame(program=i,
                               stage=j,
                               trait=k,
                               sim=1:nsim,
                               y0=tapply(out[[i]][[j]][,k], out[[i]][[j]]$sim, min),
                               y1=tapply(out[[i]][[j]][,k], out[[i]][[j]]$sim, mean),
                               y2=tapply(out[[i]][[j]][,k], out[[i]][[j]]$sim, max)))
        
      }
    }
  }
  
  # setup.
  df$program <- as.factor(df$program)
  df$stage <- as.factor(df$stage)
  df$stage <- factor(df$stage, levels=stage)
  df$trait <- as.factor(df$trait)
  df$trait <- factor(df$trait, levels=trait)
  df$sim <- as.factor(df$sim)
  
  anno <- data.frame(y=threshold, trait=trait)
  anno$trait <- as.factor(anno$trait)
  anno$trait <- factor(anno$trait, levels=trait)
  
  # plot.
  p <- ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#9999FF") +
    geom_segment(data=anno, aes(x=-Inf, xend=Inf, y=y, yend=y), color="#FF9999") +
    geom_linerange(data=df, aes(x=program, ymin=y0, ymax=y2, group=sim), position=position_dodge(width=0.7), color="#AAAAAA") +
    geom_point(data=df, aes(x=program, y=y1, group=sim), position=position_dodge(width=0.7), size=0.7, color="#555555") +
    facet_grid(rows=vars(trait), cols=vars(stage)) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(strip.background=element_blank()) +
    theme(strip.text.x=element_text(hjust=0)) +
    ylab("phenotypic value") +
    xlab("breeding program")
  
  temp <- format(Sys.time(), "%H%M%S")
  
  ggsave(plot=p,
         filename=paste(getwd(), "/breeding_exercise_", temp, ".png", sep=""),
         height=7,
         width=nbp*2*nsim/10,
         units="in",
         dpi=600)
  
  message(paste("The plot can be found at ", getwd(), "/breeding_exercise_", temp, ".png", sep=""))
  
}

### function to check the cost and if any passes the threshold in EYT.
BP.check <- function(out=out, threshold=c(3, 3, 3)){
  
  # identify the number of breeding programs.
  nbp <- length(out) - 1
  
  # get the cost.
  df <- vector()
  for(i in 1:nbp){
    temp <- out[[i]]$cost
    temp <- c(temp, sum(temp))
    names(temp)[length(temp)] <- "Total"
    df <- cbind(df, temp)
  }
  df <- data.frame(df)
  colnames(df) <- paste("Program ", 1:nbp, sep="")
  print(df)
  
  # extract only EYT.
  out2 <- out[1:nbp]
  for(i in 1:nbp) out2[[i]] <- out2[[i]]$EYT
  
  # remove any row that does not pass the threshold.
  for(i in 1:nbp) out2[[i]] <- out2[[i]][out2[[i]]$size > threshold[1] & out2[[i]]$color > threshold[2] & out2[[i]]$shape > threshold[3], ]
  
  # print some summary.
  message(paste("Out of ", nrow(out[[i]]$EYT),
                " EYT lines x simulations, the following pass the thresholds for size (", threshold[1],
                "), color (", threshold[2], ") and shape (", threshold[3], ").", sep=""))
  for(i in 1:nbp) message(paste("Program ", i, ": ", nrow(out2[[i]]), sep=""))
  
  return(out2)
  
}
