### Coefficient of coancestry
# load the packages.
library(kinship2)
library(ggplot2)
library(reshape2)

# load the pedigree.
ped <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/PERM/PERM_ped.csv", as.is=TRUE)

# manually plot the pedigree.
ped2 <- ped
ped2$x <- c(1.5,10.5,19.5,28.5, 1.5,7,12.5,18,23.5,29, 1:30)
ped2$y <- c(rep(3,4), rep(2,6), rep(1,30))
ped3 <- vector()
for(i in 5:40){
  ped3 <- rbind(ped3,
                data.frame(x1=ped2$x[i],
                           x2=ped2$x[ped2$P1[i]],
                           y1=ped2$y[i],
                           y2=ped2$y[ped2$P1[i]],
                           highlight=ifelse(i%in%c(5,6,11), TRUE, FALSE)),
                data.frame(x1=ped2$x[i],
                           x2=ped2$x[ped2$P2[i]],
                           y1=ped2$y[i],
                           y2=ped2$y[ped2$P2[i]],
                           highlight=ifelse(i%in%c(5,6,11), TRUE, FALSE)))
}
ggplot() +
  geom_segment(data=ped3, aes(x=x1, xend=x2, y=y1, yend=y2, color=highlight), linewidth=1) +
  geom_label(data=ped2, aes(x=x, y=y, label=ID), size=4) +
  theme_void() +
  theme(legend.position="none") +
  scale_color_manual(values=c("#999999", "#FF9999")) +
  coord_fixed(ratio=4)

# compute the coefficients of coancestry using R/kinship2.
coan <- kinship(id=ped$ID, dadid=ped$P1, momid=ped$P2)

# convert the matrix of coefficients for plotting.
coan4plot <- melt(data=coan)
coan4plot <- coan4plot[order(coan4plot$Var1, coan4plot$Var2), ]
colnames(coan4plot) <- c("ID1", "ID2", "coef")

# plot the coefficients of coancestry.
ggplot() +
  geom_tile(data=coan4plot, aes(x=ID1, y=ID2, fill=coef), color="#999999") +
  theme_bw() +
  scale_y_reverse(expand=c(0,0)) +
  scale_x_continuous(position="top", expand=c(0,0)) +
  scale_fill_gradient(low="#FFFFFF", high="#0000FF") +
  coord_fixed()

# calculate the coefficients of coancestry manually.
# create a table of all pairwise relationships
n <- nrow(ped)
df <- data.frame(ID1=sort(rep(1:n, n)), ID2=rep(1:n, n), manual=NA)
df <- df[!(df$ID2 > df$ID1), ]
rownames(df) <- NULL

# identify the individuals in each generation.
g0 <- c(1:4)
g1 <- c(5:10)
g2 <- c(11:40)

# compute the coefficients for g0 to itself.
for(i in g0) df$manual[df$ID1==i & df$ID2==i] <- 0.5

# compute the coefficients for g0 to g0.
for(i in g0){
  for(j in g0[g0 < i]){
    df$manual[df$ID1==i & df$ID2==j] <- 0
  }
}

# compute the coefficients for g1 to g0/g1.
for(i in g1){
  for(j in c(g0, g1)[c(g0, g1) < i]){
    a <- sort(c(ped$P1[ped$ID==i], j))
    b <- sort(c(ped$P2[ped$ID==i], j))
    df$manual[df$ID1==i & df$ID2==j] <- 0.5*(df$manual[df$ID1==a[2] & df$ID2==a[1]] + df$manual[df$ID1==b[2] & df$ID2==b[1]])
  }
}

# compute the coefficients for g1 to itself.
for(i in g1){
  a <- sort(c(ped$P1[ped$ID==i], ped$P2[ped$ID==i]))
  df$manual[df$ID1==i & df$ID2==i] <- 0.5*(1 + df$manual[df$ID1==a[2] & df$ID2==a[1]])
}

# compute the coefficients for g2 to g0/g1/g2.
for(i in g2){
  for(j in c(g0, g1, g2)[c(g0, g1, g2) < i]){
    a <- sort(c(ped$P1[ped$ID==i], j))
    b <- sort(c(ped$P2[ped$ID==i], j))
    df$manual[df$ID1==i & df$ID2==j] <- 0.5*(df$manual[df$ID1==a[2] & df$ID2==a[1]] + df$manual[df$ID1==b[2] & df$ID2==b[1]])
  }
}

# compute the coefficients for g2 to itself.
for(i in g2){
  a <- sort(c(ped$P1[ped$ID==i], ped$P2[ped$ID==i]))
  df$manual[df$ID1==i & df$ID2==i] <- 0.5*(1 + df$manual[df$ID1==a[2] & df$ID2==a[1]])
}

# compare the coefficients from R/kinship2 and manual calculation.
df$kinship2 <- coan4plot$coef[!(coan4plot$ID2 > coan4plot$ID1)]

ggplot() +
  geom_point(data=df, aes(x=manual, y=kinship2)) +
  theme_bw() +
  coord_fixed()


### Coefficient of fraternity
# load the packages.
library(AGHmatrix)

# compute the coefficients of fraternity using R/AGHmatrix.
frat <- Amatrix(data=ped, dominance=TRUE)

# convert the matrix of coefficients for plotting.
frat4plot <- melt(data=frat)
frat4plot <- frat4plot[order(frat4plot$Var1, frat4plot$Var2), ]
colnames(frat4plot) <- c("ID1", "ID2", "coef")

# plot the coefficients of coancestry.
ggplot() +
  geom_tile(data=frat4plot, aes(x=ID1, y=ID2, fill=coef), color="#999999") +
  theme_bw() +
  scale_y_reverse(expand=c(0,0)) +
  scale_x_continuous(position="top", expand=c(0,0)) +
  scale_fill_gradient(low="#FFFFFF", high="#0000FF") +
  coord_fixed()

# calculate the coefficients of fraternity manually.
# create a table of all pairwise relationships
df2 <- df[, c(1,2)]
df2$manual <- NA

# compute the coefficients for an individual to itself.
for(i in 1:n) df2$manual[df2$ID1==i & df2$ID2==i] <- 1

# compute the coefficients for g0 to g0/g1/g2.
for(i in 1:n){
  for(j in g0[g0 < i]){
    df2$manual[df2$ID1==i & df2$ID2==j] <- 0
  }
}

# compute the remaining coefficients.
for(i in c(g1,g2)){
  for(j in c(g1,g2)[c(g1,g2) < i]){
    a <- sort(c(ped$P1[i], ped$P1[j]))
    b <- sort(c(ped$P2[i], ped$P2[j]))
    d <- sort(c(ped$P1[i], ped$P2[j]))
    e <- sort(c(ped$P2[i], ped$P1[j]))
    
    if(is.na(df2$manual[df2$ID1==i & df2$ID2==j])){
      P1P1 <- df$manual[df$ID1==a[2] & df$ID2==a[1]]
      P2P2 <- df$manual[df$ID1==b[2] & df$ID2==b[1]]
      P1P2 <- df$manual[df$ID1==d[2] & df$ID2==d[1]]
      P2P1 <- df$manual[df$ID1==e[2] & df$ID2==e[1]]
      df2$manual[df2$ID1==i & df2$ID2==j] <- P1P1*P2P2 + P1P2*P2P1
    }
  }
}

# compare the coefficients from R/AGHmatrix and manual calculation.
df2$AGHmatrix <- frat4plot$coef[!(frat4plot$ID2 > frat4plot$ID1)]

ggplot() +
  geom_point(data=df2, aes(x=manual, y=AGHmatrix)) +
  theme_bw() +
  coord_fixed()


### Parentage inference
# load the package
library(qtl2)

# load the data
magic <- read_cross2("https://raw.github.com/cjyang-work/lecture/main/2024_genetic_relationship/MAGIC.zip")

# compute P(IBD)
gp <- calc_genoprob(magic)

# manually plot the P(IBD)
gp4plot <- vector()
for(i in 1:length(gp)){
  gp4plot <- rbind(gp4plot,
                   cbind(melt(gp[[i]]),
                         Chr=i,
                         Pos=sort(rep(magic$gmap[[i]], 4*10))))
}
colnames(gp4plot) <- c("RIL", "Founder", "SNP", "PIBD", "Chr", "Pos")
levels(gp4plot$Founder) <- c("B73", "Mo17", "W22", "FV2")

ggplot() +
  geom_line(data=gp4plot, aes(x=Pos, y=PIBD, color=Founder), linewidth=1) +
  facet_grid(rows=vars(RIL), cols=vars(Chr), scales="free_x", space="free_x") +
  theme_minimal() +
  theme(legend.position="bottom") +
  scale_color_manual(values=c("#CA0020", "#F4A582", "#0571B0", "#92C5DE")) +
  scale_y_continuous(breaks=c(0,0.5,1))

# threshold of 0.5
ip <- maxmarg(gp, minprob=0.5)

# manually plot the inferred parentage.
ip4plot <- data.frame(t(do.call(cbind, ip)))
temp <- sapply(1:length(magic$gmap), FUN=function(i) length(magic$gmap[[i]]))
ip4plot$Chr <- unlist(lapply(1:length(temp), FUN=function(i) rep(i, temp[i])))
ip4plot$Pos <- unlist(magic$gmap)
ip4plot <- melt(ip4plot, id.vars=c("Chr", "Pos"))
colnames(ip4plot)[3:4] <- c("RIL", "Founder")
ip4plot$Founder <- factor(x=ip4plot$Founder, labels=c("B73", "Mo17", "W22", "FV2"))

ggplot() +
  geom_point(data=ip4plot, aes(x=Pos, y=RIL, color=Founder), size=2) +
  facet_grid(cols=vars(Chr), scales="free_x", space="free_x") +
  theme_minimal() +
  theme(legend.position="bottom") +
  scale_color_manual(values=c("#CA0020", "#F4A582", "#0571B0", "#92C5DE")) +
  scale_y_discrete(limits=rev)


### Additive GRM
# load the packages
library(sommer)

# load the genomic marker data for PERM.
geno <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/PERM/PERM_geno.csv", as.is=TRUE, row.names=1)
geno <- as.matrix(geno)

# compute the A-GRM using sommer.
KA <- A.mat(geno-1)

# compute the allele frequency.
p <- colSums(geno)/(2*colSums(!is.na(geno)))

# center the marker data.
W <- geno - matrix(rep(2*p, nrow(geno)), nrow=nrow(geno), byrow=TRUE)

# obtain the denominator.
d <- sum(2*p*(1-p))

# compute the A-GRM.
KA2 <- W%*%t(W)/d

# check if KA and KA2 are identical.
KA4plot <- rbind(data.frame(sommer=KA[lower.tri(KA)], manual=KA2[lower.tri(KA2)], type="off-diagonal"),
                 data.frame(sommer=diag(KA), manual=diag(KA2), type="diagonal"))

ggplot() +
  geom_point(data=KA4plot, aes(x=sommer, y=manual, color=type)) +
  theme_bw() +
  coord_fixed()


### Dominance GRM
# compute the D-GRM using sommer.
KD <- D.mat(X=geno-1, nishio=FALSE)

# convert to dominance marker data.
genod <- geno
genod[genod==2] <- 0

# center the marker data.
Wd <- genod - matrix(rep(2*p*(1-p), nrow(genod)), nrow=nrow(genod), byrow=TRUE)

# obtain the denominator.
dd <- sum(2*p*(1-p)*( 1-2*p*(1-p) ))

# compute the D-GRM.
KD2 <- Wd%*%t(Wd)/dd

# check if KD and KD2 are identical.
KD4plot <- rbind(data.frame(sommer=KD[lower.tri(KD)], manual=KD2[lower.tri(KD2)], type="off-diagonal"),
                 data.frame(sommer=diag(KD), manual=diag(KD2), type="diagonal"))

ggplot() +
  geom_point(data=KD4plot, aes(x=sommer, y=manual, color=type)) +
  theme_bw() +
  coord_fixed()


### Compare A/D-GRM to coancestry/fraternity.
# load the packages.
library(kinship2)
library(AGHmatrix)

# load the pedigree for PERM.  
ped <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/PERM/PERM_ped.csv", as.is=TRUE)

# compute the pedigree-based A-GRM.
HA <- kinship(id=ped$ID, dadid=ped$P1, momid=ped$P2)
HD <- Amatrix(data=ped, dominance=TRUE)

# compare HA to KA, HD to KD.
compareAD <- rbind(data.frame(pedigree=HA[lower.tri(HA)],
                              marker=KA[lower.tri(KA)]/2,
                              type="off-diagonal",
                              coef="Additive/coancestry"),
                   data.frame(pedigree=diag(HA),
                              marker=diag(KA)/2,
                              type="diagonal",
                              coef="Additive/coancestry"),
                   data.frame(pedigree=HD[lower.tri(HD)],
                              marker=KD[lower.tri(KD)],
                              type="off-diagonal",
                              coef="Dominance/fraternity"),
                   data.frame(pedigree=diag(HD),
                              marker=diag(KD),
                              type="diagonal",
                              coef="Dominance/fraternity"))
temp <- data.frame(x1=c(0,0), x2=c(0.55,1), y1=c(0,0), y2=c(0.55,1), coef=c("Additive/coancestry", "Dominance/fraternity"))

ggplot() +
  geom_point(data=compareAD, aes(x=pedigree, y=marker, color=type)) +
  geom_segment(data=temp, aes(x=x1, xend=x2, y=y1, yend=y2), color="#999999", linewidth=1) +
  geom_smooth(data=compareAD, aes(x=pedigree, y=marker), method="lm", color="#000000") +
  facet_wrap(vars(coef), nrow=1) +
  theme_bw() +
  coord_fixed(xlim=c(-1,1.5))


### Computing AxA epistasis GRM
# subset the additive markers to use as an example.
X <- geno[, 1:50]

# center and scale the additive markers first - different approach from the previous way.
W <- scale(X, center=TRUE, scale=TRUE)/sqrt(ncol(X))

# compute the additive x additive epistasis markers.
WW <- lapply(1:(ncol(W)-1), FUN=function(i) W[,i]*W[,(i+1):ncol(W), drop=FALSE])
WW <- do.call(cbind, WW)

# compute the additive GRM.
MX <- W%*%t(W)

# method 1 for computing the additive x additive epistasis GRM.
MXX1 <- WW%*%t(WW)

# method 2 for computing the additive x additive epistasis GRM using trick.
MXX2 <- 0.5*(MX*MX - (W*W)%*%t(W*W))

# check if method 1 and 2 are identical.
MXX4plot <- rbind(data.frame(method.1=MXX1[lower.tri(MXX1)],
                             method.2=MXX2[lower.tri(MXX2)],
                             type="off-diagonal"),
                  data.frame(method.1=diag(MXX1),
                             method.2=diag(MXX2), type="diagonal"))

ggplot() +
  geom_point(data=MXX4plot, aes(x=method.1, y=method.2, color=type)) +
  theme_bw() +
  coord_fixed()

# compute KAA for the full data.
W <- scale(geno, center=TRUE, scale=TRUE)/sqrt(ncol(geno))
KAA <- 0.5*((W%*%t(W))*(W%*%t(W)) - (W*W)%*%t(W*W))

# compare KA to KD and KA to KAA.
K4plot <- rbind(data.frame(KA=KA[lower.tri(KA)],
                           KD=KD[lower.tri(KD)],
                           KAA=KAA[lower.tri(KAA)],
                           type="off-diagonal"),
                data.frame(KA=diag(KA),
                           KD=diag(KD),
                           KAA=diag(KAA),
                           type="diagonal"))
K4plot <- melt(K4plot, id.vars=c("KA", "type"))

ggplot() +
  geom_point(data=K4plot, aes(x=KA, y=value, color=type)) +
  facet_wrap(vars(variable), nrow=1) +
  theme_bw() +
  coord_fixed()


### Splitting and combining GRMs
# get the marker information.
marker <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/PERM/PERM_marker.csv", as.is=TRUE)

# loop through all 10 chromosomes.
d.list <- list()
KA.list <- list()
for(i in 1:10){
  # extract the marker for chr-i.
  temp.geno <- geno[, marker$Chr==i]
  
  # compute the allele frequencies.
  temp.p <- colSums(temp.geno)/(2*colSums(!is.na(temp.geno)))
  
  # center the marker data.
  temp.W <- temp.geno - matrix(2*temp.p, nrow=nrow(temp.geno), ncol=ncol(temp.geno), byrow=TRUE)
  
  # obtain the denominator.
  temp.d <- sum(2*temp.p*(1-temp.p))
  
  # compute the A-GRM.
  temp.KA <- temp.W%*%t(temp.W)/temp.d
  
  # collect the results.
  d.list[[i]] <- temp.d
  KA.list[[i]] <- temp.KA
}

# combine the GRMs.
KA.comb <- d.list[[1]]*KA.list[[1]]
for(i in 2:10) KA.comb <- KA.comb + d.list[[i]]*KA.list[[i]]
KA.comb <- KA.comb/sum(unlist(d.list))

# combine all but one chromosome.
KA.loco <- d.list[[2]]*KA.list[[2]]
for(i in 3:10) KA.loco <- KA.loco + d.list[[i]]*KA.list[[i]]
KA.loco <- KA.loco/sum(unlist(d.list)[-1])

# compare Original KA to Combined and LOCO.
KAcomb4plot <- rbind(data.frame(Original=KA[lower.tri(KA)],
                                Combined=KA.comb[lower.tri(KA.comb)],
                                LOCO=KA.loco[lower.tri(KA.loco)],
                                type="off-diagonal"),
                     data.frame(Original=diag(KA),
                                Combined=diag(KA.comb),
                                LOCO=diag(KA.loco),
                                type="diagonal"))
KAcomb4plot <- melt(KAcomb4plot, id.vars=c("Original", "type"))

ggplot() +
  geom_point(data=KAcomb4plot, aes(x=Original, y=value, color=type)) +
  facet_wrap(vars(variable), nrow=1) +
  theme_bw() +
  coord_fixed()


### Genetic distance

### Fst
# load packages.
library(hierfstat)

# load the marker data.
geno <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/DIV/DIV_geno.csv", as.is=TRUE, row.names=1)
geno <- as.matrix(geno)

# load the grouping data.
group <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/DIV/DIV_group.csv", as.is=TRUE)
table(group$Group)

# compute the Fst for each group.
fs.dosage(dos=geno, pop=group$Group)
#    CIMMYT EU_Inbred NA_Inbred    PVP Teosinte    All
#Fis 1.0000    1.0000    1.0000  1.000   1.0000 1.0000
#Fst 0.1304    0.1421   -0.0399 -0.043   0.2586 0.0896

# manually compute Fst.

# order geno by group.
group2 <- group[order(group$Group),]
geno2 <- geno[group2$ID,]

# identify the indices for each group.
idx <- lapply(unique(group2$Group), FUN=function(i) which(group2$Group==i))

# compute the matches as 1 - average Manhattan distance between individuals.
# this works for inbred data.
match.all <- 1 - as.matrix(dist(x=geno2/2, method="manhattan"))/ncol(geno2)

# set anything other than the lower triangle to 0 for convenience.
match.all[upper.tri(match.all, diag=TRUE)] <- 0

# compute match-within.
match.w <- vector()
for(i in 1:length(idx)){
  
  # extract the matrix containing all the within-i match.
  temp <- match.all[idx[[i]], idx[[i]]]
  
  # compute the mean of the match.
  match.w <- c(match.w, mean(temp[lower.tri(temp)]))
  
}

# compute match-between.
match.b <- vector()
for(i in 1:(length(idx)-1)){
  for(j in (i+1):length(idx)){
  
    # extract the matrix containing all the between-i match.
    match.b <- c(match.b, mean(c(match.all[idx[[j]], idx[[i]]])))

  }
}

# compute the mean of the match-between.
match.b <- mean(match.b)

# compute the Fst manually.
c((match.w - match.b)/(1 - match.b), (mean(match.w) - match.b)/(1 - match.b))
#0.13037383  0.14205002 -0.03989604 -0.04297350  0.25857899  0.08962666


### genome-wide Fst
# load package.
library(ggplot2)

# load the marker information.
marker <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/DIV/DIV_marker.csv", as.is=TRUE)

# compare the Fst between NA_Inbred and EU_Inbred.
# extract only the data for these two populations.
geno3 <- geno[group$Group%in%c("NA_Inbred", "EU_Inbred"), ]
group3 <- group[group$Group%in%c("NA_Inbred", "EU_Inbred"), ]

# compute genome-wide Fst.
gw.fst <- varcomp.glob(levels=group3[, "Group", drop=FALSE], loci=data.frame(geno3))
gw.fst <- gw.fst$loc[,1]/rowSums(gw.fst$loc)

# prepare the data for plotting.
gw4plot <- cbind(marker, Fst=gw.fst)

# plot the data.
ggplot() +
  geom_point(data=gw4plot, aes(x=Pos_bp, y=Fst), color="#555555") +
  facet_grid(cols=vars(Chr), switch="x", scales="free_x", space="free_x") +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

# null distribution for running Qst-Fst test.
ggplot() +
  geom_density(data=gw4plot, aes(x=Fst), fill="#DDDDDD") +
  theme_bw()


### NJ with Roger's and Nei's distance.
# load package.
library(adegenet)
library(poppr)
library(ape)

# convert the marker data into `adegenet::genind` format.
geno4 <- geno
geno4[geno4==0] <- "A/A"
geno4[geno4==1] <- "A/B" # we don't actually have heterozygotes.
geno4[geno4==2] <- "B/B"
geno4 <- df2genind(X=geno4, sep="/")

# compute Roger's distance.
dist.rog <- rogers.dist(x=geno4)

# compute Nei's distance.
dist.nei <- nei.dist(x=geno4)

# plot the NJ tree for Roger's distance.
plot(nj(dist.rog), cex=0.5)

# plot the NJ tree for Nei's distance.
plot(nj(dist.nei), cex=0.5)


### Population structure
### PCA
# load the packages
library(ggplot2)

# load the marker data.
geno <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/DIV/DIV_geno.csv", as.is=TRUE, row.names=1)
geno <- as.matrix(geno)

# load the grouping data.
group <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/DIV/DIV_group.csv", as.is=TRUE)

# compute the genetic relationship matrix.
geno.scaled <- scale(x=geno, center=TRUE, scale=TRUE)
temp <- geno.scaled%*%t(geno.scaled)/ncol(geno.scaled)

# compute the principal components.
PC <- eigen(temp)

# prepare the data to plot the percent variance explained by the PC.
pve <- data.frame(PC=1:nrow(geno), PVE=cumsum(100*PC$values/sum(PC$values)))

# plot the PVE.
ggplot() +
  geom_point(data=pve, aes(x=PC, y=PVE)) +
  theme_bw()

# prepare the data to create PCA plot.
pca <- data.frame(PC1=PC$vectors[,1],
                  PC2=PC$vectors[,2],
                  PC3=PC$vectors[,3],
                  PC4=PC$vectors[,4],
                  ID=group$ID,
                  Group=group$Group)

# plot the PCA (PC1 vs PC2, text label).
ggplot() +
  geom_text(data=pca, aes(x=PC1, y=PC2, label=ID), size=1) +
  theme_bw() +
  coord_fixed()

# plot the PCA (PC1 vs PC2, group-colored).
ggplot() +
  geom_point(data=pca, aes(x=PC1, y=PC2, color=Group), size=1) +
  theme_bw() +
  coord_fixed()

# plot the PCA (PC3 vs PC4, group-colored).
ggplot() +
  geom_point(data=pca, aes(x=PC3, y=PC4, color=Group), size=1) +
  theme_bw() +
  coord_fixed()



### Population admixture
# load the package.
library(LEA)

# set working directory
setwd("...")

# run the admixture analysis with K=2/3/4.
out <- snmf("BCF.geno", K=2:4, entropy=TRUE, repetitions=10, project="new")

# use to entropy to help us decide on the K.
plot(out, pch=16, col="#FF9999")

# K=2.
barchart(out,
         K=2,
         run=which.min(cross.entropy(out, K=2)),
         sort.by.Q=F,
         border=NA,
         space=0,
         col=c("#E69F00", "#56B4E9"),
         ylab="Ancestry proportions")

# K=3.
barchart(out,
         K=3,
         run=which.min(cross.entropy(out, K=3)),
         sort.by.Q=F,
         border=NA,
         space=0,
         col=c("#E69F00", "#56B4E9", "#009E73"),
         ylab="Ancestry proportions")

# K=4.
barchart(out,
         K=4,
         run=which.min(cross.entropy(out, K=4)),
         sort.by.Q=F,
         border=NA,
         space=0,
         col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
         ylab="Ancestry proportions")


### GWAS with/without population structure control.
# load packages.
library(sommer)

# load the marker information.
marker <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/DIV/DIV_marker.csv", as.is=TRUE)

# load the phenotype data.
pheno <- read.csv("https://raw.githubusercontent.com/cjyang-work/lecture/main/2024_genetic_relationship/DIV/DIV_pheno.csv", as.is=TRUE)

# calculate A-GRM.
KA <- A.mat(geno-1)

# GWAS with mixed model.
gw1 <- GWAS(fixed=trait~1,
            random=~vsr(ID, Gu=KA),
            rcov=~units,
            data=pheno,
            M=geno-1,
            gTerm="u:ID",
            dateWarning=FALSE)

# GWAS with fixed model.
gw2 <- sapply(1:ncol(geno), FUN=function(i) summary(glm(trait~geno[,i], data=pheno))$coefficients[2,4])

# prepare the data for plotting the mixed model results.
gw1plot <- cbind(marker[,c(2,4)], gw1$scores[,1])
colnames(gw1plot) <- c("Chrom", "Position", "p.val")
gw1plot$p.val[gw1plot$p.val==Inf] <- NA

# prepare the data for plotting the fixed model results.
gw2plot <- cbind(marker[,c(2,4)], -log10(gw2))
colnames(gw2plot) <- c("Chrom", "Position", "p.val")

# plot the mixed model results.
manhattan(gw1plot, pch=16, fdr.level=FALSE, show.fdr=FALSE)

# plot the fixed model results.
manhattan(gw2plot, pch=16, fdr.level=FALSE, show.fdr=FALSE)

# compare both.
ggplot() +
  annotate("segment", x=-Inf, xend=Inf, y=4, yend=4, color="#9999FF") +
  annotate("segment", x=4, xend=4, y=-Inf, yend=Inf, color="#9999FF") +
  geom_point(aes(x=gw1plot$p.val, y=gw2plot$p.val), na.rm=TRUE) +
  geom_abline(slope=1, intercept=0, color="#FF0000") +
  theme_bw() +
  xlab("With PS control") +
  ylab("Without PS control") +
  coord_fixed()

### eigenGWAS.
# compute the first principal component.
PC1 <- eigen(KA)$vectors[,1]

# add that to the phenotype data.
pheno <- cbind(pheno, PC1=PC1)

# eigenGWAS with mixed model.
gw3 <- GWAS(fixed=PC1~1,
            random=~vsr(ID, Gu=KA),
            rcov=~units,
            data=pheno,
            M=geno-1,
            gTerm="u:ID",
            dateWarning=FALSE)

# prepare the data for plotting the mixed model results.
gw3plot <- cbind(marker[,c(2,4)], gw3$scores[,1])
colnames(gw3plot) <- c("Chrom", "Position", "p.val")
gw3plot$p.val[gw3plot$p.val==Inf] <- NA

# plot the mixed model results.
manhattan(gw3plot, pch=16, fdr.level=FALSE, show.fdr=FALSE)


### GBLUP
# randomly sample 35 lines from EU_Inbred/NA_Inbred/PVP to predict 18 lines in Teosinte.
# identify the lines.
set.seed(252)
train1 <- sort(sample(group$ID[!(group$Group=="Teosinte")], 96))
train2 <- group$ID[group$Group=="PVP"]
test <- group$ID[group$Group=="Teosinte"]

# prepare the pheno data frame.
pheno <- pheno[, 1:2]
rownames(pheno) <- pheno$ID

# random 96 + Teosinte.
KA1 <- A.mat(geno[c(train1, test), ]-1)
pheno1 <- pheno[c(train1, test), ]
pheno1$trait[pheno1$ID%in%test] <- NA

# PVP + Teosinte.
KA2 <- A.mat(geno[c(train2, test), ]-1)
pheno2 <- pheno[c(train2, test), ]
pheno2$trait[pheno2$ID%in%test] <- NA

# GBLUP 1, get the predicted values.
mm1 <- mmer(fixed=trait~1,
            random=~vsr(ID, Gu=KA1),
            rcov=~units,
            data=pheno1,
            dateWarning=FALSE)
pred1 <- c(mm1$Beta$Estimate + mm1$U$`u:ID`$trait)[test]

# GBLUP 2, get the predicted values.
mm2 <- mmer(fixed=trait~1,
            random=~vsr(ID, Gu=KA2),
            rcov=~units,
            data=pheno2,
            dateWarning=FALSE)
pred2 <- c(mm2$Beta$Estimate + mm2$U$`u:ID`$trait)[test]

# compare the prediction accuracies.
cor(cbind(obs=pheno[test, "trait"], pred1, pred2))
