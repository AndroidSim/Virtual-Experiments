#! /usr/bin/Rscript
###
## Slices out the hcount-????.csv datat within the given band, and prints the
## average number of Hepatocytes in that band.
##
## Time-stamp: <2019-07-18 16:24:29 gepr>
###
argv <- commandArgs(TRUE)

if (length(argv) < 4 || !(argv[1] == "dCV" || argv[1] == "dPV")) {
  print("Usage: hcounts-inband.r <dCV or dPV> min max exp1 exp2 ...")
  q()
}

direction <- argv[1]
dmin <- as.numeric(argv[2])
dmin_ord <- dmin+1
dmax <- as.numeric(argv[3])
dmax_ord <- dmax+1
exps <- argv[-(1:3)]

for (exp in exps) {
  files <- list.files(path = exp, pattern=paste("hcount-",direction,"-[0-9]+.csv",sep=""), recursive=T)
  total <- 0
  for (f in files) {
    fn <- paste(exp,"/",f,sep="")
    d <- read.csv(fn, colClasses="numeric")
    if (dmin_ord > ncol(d)) next
    thismax <- ifelse(dmax_ord > ncol(d), ncol(d), dmax_ord-1)
    row <- d[1,dmin_ord:thismax]
    rowtot <- sum(row)
    ##    print(cbind(row[1,], rowtot))
    total <- total+rowtot
  }
  avg <- total/length(files)
  print(paste("\"interval\"", ",", " \"µ(#vHPCs)\"", sep=""), quote=F);
  print(paste("\"", direction,"∈[",dmin, ",", dmax, ")\", ", avg, sep=""), quote=F)
}

