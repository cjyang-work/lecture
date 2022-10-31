# load the required packages.
library(AlphaSimR)
library(rrBLUP)
library(ggplot2)

# set working directory.
setwd("")

# load the required functions.
#source("breeding_exercise_v2_modified_to_fit_old_AlphaSimR.R")
source("breeding_exercise_v2.R")

# read the file containing parameters for 3 default breeding programs.
BP <- read.csv("BP.csv", as.is=TRUE)

# run the simulation to test the 3 breeding programs.
out <- BP.eval(va=c(1,1,1),
               vae=c(0,0,0),
               vres=c(1,1,1),
               ca=0.5,
               ew=c(1,1),
               BP=BP,
               nsim=10)

# plot the results.
BP.plot(out=out, threshold=c(2, 2, 2), nsim=10)

# check the cost and identify if any line passes the threshold.
check <- BP.check(out=out, threshold=c(2, 2, 2))
check


