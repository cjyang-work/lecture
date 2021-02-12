### R scripts for "IBD, IBS, Genetic Distance, Population Structure"

## Part 1 - IBD practices
#install.packages("qtl2")
#install.packages("AlphaSimR")
#install.packages("ggplot2")

library(qtl2)
library(AlphaSimR)
library(ggplot2)

## Part 1.1 - IBD calculation in a simulated bi-parental RIL population.
set.seed(99999)

fhap <- matrix(c(rep(0,101),rep(1,101)), nrow=2, byrow=T)
pos <- seq(0,1,0.01)
names(pos) <- paste("SITE", 1:101, sep="_")

founder <- newMapPop(genMap=list(pos),
                     haplotypes=list(fhap),
                     inbred=T,
                     ploidy=2L)
SP <- SimParam$new(founder)
F0 <- newPop(founder, simParam=SP)

F1 <- makeCross(pop=F0, crossPlan=matrix(c(1,2),nrow=1), nProgeny=1, simParam=SP)
RIL <- self(pop=F1, nProgeny=50)
for(i in 1:3) RIL <- self(pop=RIL, nProgeny=1)

geno <- pullSegSiteGeno(pop=RIL)
geno[geno==1] <- 9
geno[geno==0] <- 1
geno[geno==9] <- 0

cross_info <- matrix(, nrow=nrow(geno), ncol=0)
class(cross_info) <- "integer"
rownames(cross_info) <- rownames(geno)

xdata <- list(crosstype="riself",
              geno=list(geno),
              gmap=list(pos*100),
              is_x_chr=F,
              is_female=rep(F,50),
              cross_info=cross_info,
              alleles=c("A", "B"))
class(xdata) <- "cross2"

xdata.gp <- calc_genoprob(xdata)

unname(geno[1,])

unname(geno[2,])

id <- 1
plot(x=pos*100, y=xdata.gp[[1]][id,1,], pch=16, xlab="genetic position (centiMorgan)", ylab="parent 1 probability")

id <- 2
plot(x=pos*100, y=xdata.gp[[1]][id,1,], pch=16, xlab="genetic position (centiMorgan)", ylab="parent 1 probability")


## Part 1.2 - IBD calculation in a sorghum MAGIC population.

xdata <- read_cross2("https://raw.github.com/cjyang-sruc/practice/main/data/sorghum_chr1.zip")

xdata.gp <- calc_genoprob(xdata)

#xdata$geno[[1]] <- xdata$geno[[1]][, seq(1,6933,70)]
#xdata$gmap[[1]] <- xdata$gmap[[1]][seq(1,6933,70)]
#xdata$founder_geno[[1]] <- xdata$founder_geno[[1]][, seq(1,6933,70)]
#xdata.gp <- calc_genoprob(xdata)

plot.gp <- function(gp, gmap, id, chr=1){

  n <- dim(gp[[chr]])[2]
  temp <- sort(rep(1:n, dim(gp[[chr]])[3]))
  
  dat <- gp[[chr]][id, ,]
  dat <- data.frame(probability=c(t(dat)), founder=temp, position=rep(gmap[[chr]], n))
  dat$founder <- as.factor(dat$founder)
    
  ggplot() +
    geom_line(data=dat, aes(x=position, y=probability, color=founder)) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(legend.position="bottom") +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#CCCCCC") +
    guides(color=guide_legend(nrow=2))
  
}

plot.gp(gp=xdata.gp, gmap=xdata$gmap, id=71)

xdata.fcall <- maxmarg(xdata.gp, minprob=0.5)

plot.fcall <- function(fcall, id, chr=1){
  
  dat <- data.frame(marker=1:ncol(fcall[[chr]]), founder=fcall[[chr]][id,])

  ggplot() +
    geom_point(data=dat, aes(x=marker, y=founder), size=1, na.rm=T) +
    theme(panel.background=element_blank(), panel.grid.minor=element_blank()) +
    theme(panel.grid.major.x=element_blank()) +
    theme(panel.grid.major.y=element_line(color="#CCCCCC", linetype=2)) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#CCCCCC") +
    scale_y_continuous(limits=c(1,19), breaks=c(1:19))
  
}

plot.fcall(fcall=xdata.fcall, id=71)


## Part 2 - IBS practices

#install.packages("rrBLUP")
library(rrBLUP)

map <- read.csv("https://raw.github.com/cjyang-sruc/practice/main/data/soybean_map.csv", as.is=T)
info <- read.csv("https://raw.github.com/cjyang-sruc/practice/main/data/soybean_info.csv", as.is=T)
geno <- read.csv("https://raw.github.com/cjyang-sruc/practice/main/data/soybean_geno.csv", as.is=T, row.names=1)

geno <- t(as.matrix(geno))

kin <- A.mat(geno-1)

p <- (colSums(geno==2) + 0.5*colSums(geno==1))/nrow(geno)
W <- geno - matrix(rep(2*p, nrow(geno)), nrow=nrow(geno), byrow=T)
d <- sum(2*p*(1-p))
kin2 <- W%*%t(W)/d

plot(x=kin[lower.tri(kin)],
     y=kin2[lower.tri(kin2)],
     pch=16,
     cex=0.5,
     col="#555555",
     xlab="kin",
     ylab="kin2",
     xlim=range(kin),
     ylim=range(kin2))

points(x=diag(kin),
       y=diag(kin2),
       pch=16,
       cex=0.5,
       col="#55FF55")

abline(a=0, b=1, col="#FF5555")

cor(kin[lower.tri(kin)], kin2[lower.tri(kin2)])

cor(diag(kin), diag(kin2))

geno.1 <- geno[, map$chr==1]
p.1 <- (colSums(geno.1==2) + 0.5*colSums(geno.1==1))/nrow(geno.1)
W.1 <- geno.1 - matrix(rep(2*p.1, nrow(geno.1)), nrow=nrow(geno.1), byrow=T)
d.1 <- sum(2*p.1*(1-p.1))
kin3.1 <- W.1%*%t(W.1)/d.1

geno.2 <- geno[, map$chr==2]
p.2 <- (colSums(geno.2==2) + 0.5*colSums(geno.2==1))/nrow(geno.2)
W.2 <- geno.2 - matrix(rep(2*p.2, nrow(geno.2)), nrow=nrow(geno.2), byrow=T)
d.2 <- sum(2*p.2*(1-p.2))
kin3.2 <- W.2%*%t(W.2)/d.2

kin3 <- kin3.1*d.1/d + kin3.2*d.2/d
cor(kin[lower.tri(kin)], kin2[lower.tri(kin3)])
cor(diag(kin), diag(kin3))

z.1 <- 10
z.2 <- 1
kin4 <- kin3.1*d.1*z.1/(z.1*d.1+z.2*d.2) + kin3.2*d.2*z.2/(z.1*d.1+z.2*d.2)

plot(x=kin[lower.tri(kin)],
     y=kin4[lower.tri(kin4)],
     pch=16,
     cex=0.5,
     col="#555555",
     xlab="kin",
     ylab="kin4",
     xlim=range(kin),
     ylim=range(kin4))

points(x=diag(kin),
       y=diag(kin4),
       pch=16,
       cex=0.5,
       col="#55FF55")

abline(a=0, b=1, col="#FF5555")

cor(kin[lower.tri(kin)], kin4[lower.tri(kin4)])

cor(diag(kin), diag(kin4))

cor(kin3.1[lower.tri(kin3.1)], kin4[lower.tri(kin4)])
cor(diag(kin3.1), diag(kin4))

cor(kin3.2[lower.tri(kin3.2)], kin4[lower.tri(kin4)])
cor(diag(kin3.2), diag(kin4))


## Part 3 - Genetic distance practices

#install.packages("hierfstat")
library(hierfstat)

library(rrBLUP)
map <- read.csv("https://raw.github.com/cjyang-sruc/practice/main/data/soybean_map.csv", as.is=T)
info <- read.csv("https://raw.github.com/cjyang-sruc/practice/main/data/soybean_info.csv", as.is=T)
geno <- read.csv("https://raw.github.com/cjyang-sruc/practice/main/data/soybean_geno.csv", as.is=T, row.names=1)
geno <- t(as.matrix(geno))
kin <- A.mat(geno-1)

fs.dosage(dos=geno, pop=info$type)

mean(diag(kin) - 1)

geno.C <- geno[info$type=="Cultivar", ]
geno.L <- geno[info$type=="Landrace", ]
geno.W <- geno[info$type=="Wild", ]

piw.C <- vector()
for(i in 1:(nrow(geno.C)-1)){
  for(j in (i+1):nrow(geno.C)){
    piw.C <- c(piw.C, sum(abs(geno.C[i,] - geno.C[j,])))
  }
}
piw.C <- mean(piw.C)

piw.L <- vector()
for(i in 1:(nrow(geno.L)-1)){
  for(j in (i+1):nrow(geno.L)){
    piw.L <- c(piw.L, sum(abs(geno.L[i,] - geno.L[j,])))
  }
}
piw.L <- mean(piw.L)

piw.W <- vector()
for(i in 1:(nrow(geno.W)-1)){
  for(j in (i+1):nrow(geno.W)){
    piw.W <- c(piw.W, sum(abs(geno.W[i,] - geno.W[j,])))
  }
}
piw.W <- mean(piw.W)

pib.CL <- vector()
for(i in 1:nrow(geno.C)){
  for(j in 1:nrow(geno.L)){
    pib.CL <- c(pib.CL, sum(abs(geno.C[i,] - geno.L[j,])))
  }
}
pib.CL <- mean(pib.CL)

pib.CW <- vector()
for(i in 1:nrow(geno.C)){
  for(j in 1:nrow(geno.W)){
    pib.CW <- c(pib.CW, sum(abs(geno.C[i,] - geno.W[j,])))
  }
}
pib.CW <- mean(pib.CW)

pib.LW <- vector()
for(i in 1:nrow(geno.L)){
  for(j in 1:nrow(geno.W)){
    pib.LW <- c(pib.LW, sum(abs(geno.L[i,] - geno.W[j,])))
  }
}
pib.LW <- mean(pib.LW)

pib <- mean(c(pib.CL, pib.CW, pib.LW))

(pib - piw.C)/pib
(pib - piw.L)/pib
(pib - piw.W)/pib
(3*pib - piw.C - piw.L - piw.W)/(3*pib)


## Part 4 - Population structure practices

#### Part 4.1 - Principal component analysis (PCA)

pca <- eigen(kin)

round(pca$values/sum(pca$values)*100, 2)

library(ggplot2)
dat <- data.frame(PC1=pca$vectors[,1], PC2=pca$vectors[,2], type=info$type, origin=info$origin)
ggplot() +
  geom_point(data=dat, aes(x=PC1, y=PC2, color=type)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#CCCCCC")

library(ggplot2)
dat <- data.frame(PC1=pca$vectors[,1], PC2=pca$vectors[,2], type=info$type, origin=info$origin)
dat$origin2 <- ifelse(dat$origin=="Korea", "Korea", "Elsewhere")
ggplot() +
  geom_point(data=dat, aes(x=PC1, y=PC2, color=type, pch=origin2)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#CCCCCC")


#### Part 4.2 - Population admixture analysis

#install.packages("BiocManager")
#library(BiocManager)
#install("LEA")
library(LEA)

LEA.info <- info[order(info$type), ]
LEA.geno <- geno[LEA.info$line,]
LEA.geno <- data.frame(t(LEA.geno))
LEA.geno <- do.call(paste, c(LEA.geno, sep=""))
write.table(LEA.geno, "LEA.geno", quote=F, row.names=F, col.names=F, sep="\t")

out <- snmf("LEA.geno", K=3:5, entropy=T, repetitions=10, project="new")

barchart(out,
         K=3,
         run=which.min(cross.entropy(out, K=3)),
         sort.by.Q=F,
         border=NA,
         space=0,
         col=c("#E69F00", "#56B4E9", "#009E73"),
         ylab="Ancestry proportions")
abline(v=c(47, 55), col="#000000", lwd=2)

outp <- barchart(out,
                 K=4,
                 run=which.min(cross.entropy(out, K=4)),
                 sort.by.Q=F,
                 border=NA,
                 space=0,
                 col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                 ylab="Ancestry proportions")
abline(v=c(47, 55), col="#000000", lwd=2)

outp <- barchart(out,
                 K=5,
                 run=which.min(cross.entropy(out, K=5)),
                 sort.by.Q=F,
                 border=NA,
                 space=0,
                 col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"),
                 ylab="Ancestry proportions")
abline(v=c(47, 55), col="#000000", lwd=2)
