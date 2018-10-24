#! /usr/bin/Rscript
##
# Read multiple *.csv files and plot each column vs the 1st.
#
# Time-stamp: <2018-09-07 17:03:22 gepr>
#

SEMILOG <- F
plot.data <- T
plot.svg <- T

argv <- commandArgs(TRUE)

if (length(argv) < 2) {
    print("Usage: cmp-by-col.r <analysis .csv file> <analysis .csv file>")
    print("  e.g. cmp-by-col.r x00[1-6]_body.csv y00[1-6]_hsolute-dCV.csv")
    print("Note that columns must match across all files.")
    quit()
}

source("~/R/misc.r")

# determine # of plots
nplots <- length(argv)
plot.cols <- round(sqrt(nplots))
# add a new row if we rounded up
plot.rows <- ifelse(plot.cols >= sqrt(nplots), plot.cols, plot.cols+1)

#
# test for and create graphics subdirectory
#
if (!file.exists("graphics")) dir.create("graphics")

data <- vector("list")
titles <- vector("list")
## get component name from basename of 1st argument
argv.basename <- basename(argv[1])
compname <- substr(argv.basename,regexpr('_',argv.basename)+1,nchar(argv.basename))
fileName.base <- substr(compname, 0, regexpr('(_|.csv)', compname)-1)
expnames <- ""
print(fileName.base)
filenum <- 1
for (f in argv) {
  nxtName <- substr(basename(f),0,regexpr('_',basename(f))-1)
  titles[[filenum]] <- nxtName
  ##fileName.base <- paste(fileName.base,nxtName,sep="-")
  expnames <- paste(expnames,nxtName,sep="-")
  data[[filenum]] <- read.csv(f)

  filenum <- filenum+1
}

## assume all Time vectors are the same

columns <- colnames(data[[1]])
column.1 <- columns[1]
max.1 <- max(data[[1]][column.1])
##max.1 <- 3000 # hard-code a max X value (cycle)

pb <- txtProgressBar(min=0,max=length(columns),style=3)
setTxtProgressBar(pb,1)

for (column in columns[2:length(columns)]) {
  skip <- FALSE

  ## get max/min of this column over all data sets if it exists
  min.2 <- Inf
  max.2 <- -1 # init max.2
  for (df in data) {
    if (!is.element(column,colnames(df))) skip <- T
    else {
      min.2 <- min(min.2, min(df[column], na.rm=T))
      max.2 <- max(max.2, max(df[column], na.rm=T))
    }
  }
  if (SEMILOG && min.2 <= 0) min.2 <- 1.0
  ##if (skip) next  # skip columns that don't exist in all files

  ##print(paste("Working on",column,"..."))
  fileName <- paste("graphics/", fileName.base, "-", column, expnames,
                    ifelse(plot.data, "-wd", ""), sep="")

  if (nchar(fileName) > 255) {
    library(digest)
    fileName <- paste("graphics/", fileName.base, "-", column, digest(expnames),
                      ifelse(plot.data, "-wd", ""), sep="")
  }

  if (plot.svg) {
    svg(paste(fileName,".svg",sep=""), width=10, height=10)
  } else {
    png(paste(fileName,".png",sep=""), width=1600, height=1600)
  }

  ## set margins and title, axis, and label font sizes
  par(mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
  par(mfrow=c(plot.rows,plot.cols))


  ## plot this column from all data sets
  ndx <- 1
  for (df in data) {
    attach(df)
    if (exists(column)) {
      dat <- cbind(df$Time[df$Time<max.1],get(column)[1:length(df$Time[df$Time<max.1])])
    } else {
      dat <- cbind(df$Time[df$Time<max.1],rep(0,length(df$Time[df$Time<max.1])))
    }
    detach(df)
    colnames(dat) <- c(column.1, column)

    ## create the moving average data set
    dat.tmp <- dat
    dat.tmp[is.na(dat.tmp)] <- 0 # replace NAs with zeros
    dat.ma <- ma.cent(dat.tmp[,2:ncol(dat.tmp)], n=181)
    dat.ma <- cbind(dat.tmp[,1], dat.ma)
    dat.ma <- as.data.frame(dat.ma)
    colnames(dat.ma) <- colnames(dat)

    suppressWarnings(
        plot(dat.ma, main=titles[[ndx]], xlim=c(0,max.1), ylim=c(min.2,max.2), log=ifelse(SEMILOG,'y',''), type="l")
    )
    if (plot.data) points(dat, pch="·")
    ## minor.tick(nx=5,ny=5)
    grid()
    ndx <- ndx+1
  }

  setTxtProgressBar(pb,getTxtProgressBar(pb)+1)
}

close(pb)

