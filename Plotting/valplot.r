#! /usr/bin/Rscript
##
# Read multiple *.csv files and plot each column vs the 1st.
#
# Time-stamp: <2018-07-27 11:56:45 gepr>
#
argv <- commandArgs(TRUE)

if (length(argv) != 3) {
  print("Usage: valplot.r <prefix> <dat column> <exp column>")
  print("       e.g. valplot.r revals218x s_at Sucrose")
  print("       <dat column> can also be \"Average\" or \"MinMax\"")
  print("Note: The script expects both the <prefix> and reduced data")
  print("      directories to exist.")
  quit()
}

prefix <- argv[1]
refName <- argv[2]
datFile <- paste(prefix,"/","datModel-0000.csv",sep="")
expName <- argv[3]
expFile <- paste(prefix,"-reduced","/",prefix,"_avgOutput.csv",sep="")

########################################################################
##  Constants
##
MAGIC <- 12   # I don't know why we have to scale it????

########################################################################
##  Options
##
USETITLE <- TRUE
USELEGEND <- T
XTICLABELS <- TRUE
YTICLABELS <- TRUE
YLABVALUE <- paste(expName,"Fraction") # string or ""
XLABVALUE <- "Time" # string or ""
LOGVALUE <- "y" # "x", "y", or ""
RADIUS <- 0.56 # Default in case multiple ref profiles aren't used
PLOTNOMINAL <- FALSE
########################################################################

## extract the dose
dp <- scan(paste(prefix,"/cfg/delivery.properties",sep=""),what=character())
## get the percentage
dp.index <- grep(expName,dp) # find the right vCompound
dp.perc <- dp[dp.index+2] # get the dose percentage
dp.perc <- gsub("[^0-9.]", "", dp.perc) # remove any extra characters like ;
if (dp.perc == "") {
  print(paste("Couldn't find Dose percentage for", refName))
  q()
}
dp.perc <- as.numeric(dp.perc)
## get the total
dp.total.index <- grep("referenceDose",dp)
if (length(dp.total.index) > 1) {
  for (i in dp.total.index) {
    candidate <- dp[i+2]
    candidate <- gsub("[^0-9.]", "", candidate)
    if (candidate == "") next
    dp.total <- as.numeric(candidate)
    break
  }
}
if (!exists("dp.total")) {
  print("Couldn't find referencDose")
  q()
}
DOSE <- dp.perc*dp.total

## get output data
dat <- read.csv(datFile)
exp <- read.csv(expFile)

library(matrixStats) # for rowSds(), rowMax(), etc.
## get the radius
rawDat <- cbind(dat["Time"])
numobs <- nrow(exp)
## remove that annoying extra column
##dat <- dat[1:ncol(dat)-1]
for (colNdx in colnames(dat)) {
  ##print(colNdx)
  if (!is.na(charmatch("s_", colNdx))) {
    rawDat <- cbind(rawDat, dat[colNdx])
  }
}
rawDat <- as.matrix(rawDat)
rawDat[rawDat == 0.0] <- NA
means <- rowMeans(rawDat[,2:ncol(rawDat)], na.rm=TRUE)
stddevs <- rowSds(as.matrix(rawDat[,2:ncol(rawDat)]),na.rm=T)
cvs <- stddevs/means
RADIUS <- mean(cvs[!is.nan(cvs)])

if (refName != "Average" && refName != "MinMax") {
  ## else use the given time series from the dat dataset
  refNdx <- match(refName, colnames(dat))
  if (is.na(refNdx)) {
    print(paste("Warning!!! Couldn't find",refName,"in dat_results_0.csv."))
    print("      Continuing using Sucrose averages.")
    refName <- "Average"
  } else {
    print(paste("Comparing ", expName," to ", refName,".", sep=""))
    dat.final <- dat[refName]
  }
  ## upper and lower sliced from "initial" data frame
  upper <- as.matrix(((1.0 + RADIUS) * dat.final))
  lower <- as.matrix(((1.0 - RADIUS) * dat.final))
  ## data.frame

  colnames(upper) <- paste(refName,"+",format(RADIUS,digits=4),sep="")
  colnames(lower) <- paste(refName,"-",format(RADIUS,digits=4),sep="")
}

if (refName == "Average") {
  ## if refName == "Average" then get a mean
  print("Comparing exp Sucrose to data average ± CV.")

  dat.final <- as.matrix(means)
  ## define upper and lower bands from the dat series
  upper <- ((1.0 + RADIUS) * dat.final)
  lower <- ((1.0 - RADIUS) * dat.final)
  ## matrix
  colnames(upper) <- paste("μ+",format(RADIUS,digits=4),sep="")
  colnames(lower) <- paste("μ-",format(RADIUS,digits=4),sep="")
}

if (refName == "MinMax") {
  print("Comparing exp Sucrose to Min-Max band reference data.")
  upper <- as.matrix(rowMaxs(as.matrix(rawDat[,2:ncol(rawDat)]),na.rm=T))
  lower <- as.matrix(rowMins(as.matrix(rawDat[,2:ncol(rawDat)]),na.rm=T))
  middle <- upper - lower
  dat.final <- middle
  ## vectors
  colnames(upper) <- "data max"
  colnames(lower) <- "data min"
}

exp.final <- exp[expName] # don't smooth

## convert to dose fraction
exp.final <- exp.final/DOSE*MAGIC

##readline()

## handle different sizes
times <- as.matrix(dat["Time"])
if (nrow(dat.final) > nrow(exp.final)) {
  ## extend exp.final to the size of dat.final
  exp.final <- c(as.matrix(exp.final), vector("numeric",length=(nrow(dat.final)-nrow(exp.final))))
  times <- as.matrix(dat["Time"])
} else {
  if (nrow(exp.final) > nrow(dat.final)) {
    ## extend dat.final
    dat.final <- c(as.matrix(dat.final), vector("numeric",length=(nrow(exp.final)-nrow(dat.final))))
    times <- as.matrix(exp["Time"])
  }
}

##initial <- as.data.frame(cbind(times, as.numeric(as.character(dat.final)), exp.final))
initial <- cbind(times, dat.final, exp.final)
initial <- as.data.frame(initial)
colnames(initial) <- c("Time", "Dat", "Exp")

attach(initial)
dMaxI <- which.max(Dat)
rMaxI <- which.max(Exp)
maxI <- rMaxI
maxISeries <- Exp

if (maxI < dMaxI) {
  maxI <- dMaxI
  maxISeries <- Dat
}
prePeak <- cbind(Time[row(as.matrix(maxISeries)) <= maxI],
                 Dat[row(as.matrix(maxISeries)) <= maxI],
                 Exp[row(as.matrix(maxISeries)) <= maxI])
colnames(prePeak) <- c("Time", "Dat", "Exp")
postPeak <- initial[row(as.matrix(maxISeries)) > maxI,]
detach(initial)
#
# create columns with only a few points for pretty plots
#
#postPeak <- postPeak[row(as.matrix(postPeak[[1]]))%%3 == 0,]
colnames(postPeak) <- c("Time", "Dat", "Exp")
final.sub <- rbind(prePeak, postPeak)

# set y range
y.min <- 0.00001
y.max <- max(final.sub[,2:ncol(final.sub)], na.rm=TRUE)
upper.max <- max(upper)
y.max <- max(y.max, upper.max)

# set x range
x.max <- max(final.sub[1])

#
# test for and create graphics subdirectory
#
if (!file.exists("graphics")) dir.create("graphics")
fileName <- paste(prefix, refName, "vs", expName, sep="-")
attach(final.sub)

#
# get the similarity score
#
t <- cbind(final.sub, upper, lower)
colnames(t) <- c(colnames(final.sub),"upper","lower")

within <- t[t["lower"] < t["Exp"] & t["Exp"] < t["upper"],] # note the comma!
SM <- nrow(within)/nrow(t)
print(paste("Similarity Measure ", nrow(within), "/", nrow(t), " = ", SM))

#jpeg(paste("graphics/",fileName,".jpg",sep=""), quality=100, width=864, height=432)

#xfig(paste("graphics/",fileName,".fig",sep=""), width=12, height=6, onefile=TRUE)

#X11(fonts=c(X11Font(helvetica.narrow), X11Font(symbol)))

#png(paste("graphics/",fileName,".png",sep=""),
#    width=864, height=432,
#    fonts=c(X11Font(helvetica), X11Font(symbol)))

#png(paste("graphics/",fileName,".png",sep=""),
##    width=864, height=432) # publications
##    width=500, height=250) # parameter search

svg(paste("graphics/",fileName,".svg",sep=""),
    width=6, height=3)

par(family="sans")

# set margins
topmar <- ifelse(USETITLE == TRUE, 2.5, 0.1)
lefmar <- ifelse(YTICLABELS == TRUE, 4.7, 0.2)
lefmar <- ifelse(YLABVALUE == "", lefmar-2, lefmar)
botmar <- ifelse(XTICLABELS == TRUE, 4.5, 0.2)
botmar <- ifelse(XLABVALUE == "", botmar-2, botmar)
par(mar=c(botmar,lefmar,topmar,1)) # b,l,t,r

# margin within the plot region
par(xaxs="i")
par(mex=0.5)
par(mgp=c(3,1,0))

legstr <- c()
legcol <- c()
leglty <- numeric()
leglwd <- numeric()
legsym <- numeric()
count <- 1

# set symbol and font sizes
par(cex=0.68)
par(cex.axis=1.2)
par(cex.lab=1.2)

# make band thick and gray
par(lwd=3)
par(col=gray(0.9))
par(pch=16)

suppressWarnings(
    plot(initial[["Time"]], upper, type="l",
         log=LOGVALUE, axes=FALSE,
         ##frame.plot=TRUE,
         tck=0, ylab=YLABVALUE, xlim=c(0,x.max),
         ylim=range(c(y.min,y.max)),
         xlab=XLABVALUE)
)

legstr[count] <- colnames(upper)
legcol[count] <- par("col")
leglty[count] <- par("lty")
leglwd[count] <- par("lwd")
legsym[count] <- NULL
count <- count + 1

lines(initial[["Time"]], lower)
legstr[count] <- colnames(lower)
legcol[count] <- par("col")
leglty[count] <- par("lty")
leglwd[count] <- par("lwd")
legsym[count] <- NULL
count <- count + 1

par(col="black")
par(lwd=1)
par(lty="solid")

## # plot the band and data

if (PLOTNOMINAL == TRUE) {
  lines(final.sub[["Time"]], final.sub[["Dat"]], col="black")
  legstr[count] <- refName
  legcol[count] <- par("col")
  leglty[count] <- par("lty")
  leglwd[count] <- par("lwd")
  legsym[count] <- -1
  count <- count + 1
}

# plot the experimental data
points(final.sub[["Time"]], final.sub[["Exp"]])

legstr[count] <- "Analog"
legcol[count] <- par("col")
leglty[count] <- "blank"
leglwd[count] <- par("lwd")
legsym[count] <- par("pch")
count <- count + 1

# draw the box around the plot (perhaps replaces fram.plot=TRUE ?)
box()

if (USETITLE == TRUE)
  title(paste(prefix, expName, "SM =", format(SM,digits=4)))

if (USELEGEND == TRUE)
  legend("topright",legend=legstr,col=legcol, lty=leglty, lwd=leglwd, pch=legsym)

# minor X tics
axis(side=1, at=seq(0,x.max), tck=0.01, labels=FALSE, tick=TRUE)
# major X tics
axis(side=1, at=seq(0,x.max,10), tck=0.02, labels=XTICLABELS, tick=TRUE)

if (LOGVALUE == "y") {
  minortics <- append(seq(10^-5,10^-4,10^-5),
                      seq(10^-4,10^-3,10^-4))
  minortics <- append(minortics, seq(10^-3,10^-2,10^-3))
  minortics <- append(minortics, seq(10^-2,10^-1,10^-2))
  minortics <- append(minortics, seq(10^-1,1,10^-1))
  axis(side=2, at=minortics, tck=0.01, labels=FALSE, tick=TRUE)

  majortics <- c(10^-5, 10^-4, 10^-3, 10^-2, 10^-1)

  lab1 <- expression({10^{phantom()-5}})
  lab2 <- expression({10^{phantom()-4}})
  lab3 <- expression({10^{phantom()-3}})
  lab4 <- expression({10^{phantom()-2}})
  lab5 <- expression({10^{phantom()-1}})

  if (YTICLABELS == TRUE) {
    axis(side=2, at=majortics, tck=0.02,
         labels=c(lab1, lab2, lab3, lab4, lab5))
  } else {
    axis(side=2, at=majortics, tck=0.02, labels=FALSE)
  }

} else {
  axis(side=2, at=seq(0,y.max+0.10,0.1), tck=0.02, labels=YTICLABELS)
}

par(lty="solid")

