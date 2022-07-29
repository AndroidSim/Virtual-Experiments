#! /usr/bin/Rscript
###
## Time-stamp: <2021-01-29 15:07:25 Ryan Kennedy>
##
## Gathers d[CP]V data on RxnProbMap elements in hepinit-*.json
## files.  Write them for each trial and a final table with experiment
## totals. This is built from hpcrpstats.r.
##
## Note that the hepinit-????.json output contains the (potentially
## false) d[CP]V IVars (i.e. instance variables) inside the Hepatocyte
## objects.  If hepInitRead=true, then the d[CP]V values will not
## reflect the true values given the HepStruct used. However, these
## RxnProb stats will *still* be correct, commensurate with the d[CP]V
## values found in these hepinit-*.json files, because they are
## calculated based on what the Hepatocyte *thinks* its d[CP]V are.
##
###
argv <- commandArgs(T)

usage <- function() {
  print("Usage: pmetabdcvdpv.r <exp directories>")
  print("  directories should contain files like hepinit-????.json")
  quit()
}
if (length(argv) < 1) usage()

library(stats)
library(matrixStats)

###
## test for and create graphics subdirectory
###
if (!file.exists("graphics")) dir.create("graphics")

exps <- argv
filebase <- "hepinit"

whole <- vector()
for (x in exps) {
	cat(paste("Processing",x,"\n"))
	if (!file.exists(x)) {
		print(paste(x,"doesn't exist."))
		quit("no")
	}
	pattern <- paste(filebase,"-[0-9]+.json",sep="")
	files <- list.files(path=x, pattern=pattern, recursive=F, full.names=T)
	pb <- txtProgressBar(min=0,max=length(files),style=3)
	setTxtProgressBar(pb,0)

	dat.all <- numeric()
	for (f in files) {
		tdat <- numeric()
		trial.name <- substr(f,regexpr("[0-9]+.json",f),regexpr(".json",f)-1)
		hi <- jsonlite::fromJSON(f)
		hi.hr <- hi$hepatocyte_records
		## assume the enzyme groups are all the same for all vHPCs!
		if (!exists("egs")) egs <- names(hi.hr[[1]]$RxnProbMap)

		vHPC.ids <- names(hi.hr[1:length(hi.hr)])
		for (id in vHPC.ids) {
			dcv <- hi.hr[id][[1]]$dCV
			dpv <- hi.hr[id][[1]]$dPV
			rps <- hi.hr[id][[1]]$RxnProbMap
			tdat.h <- cbind(dcv, dpv, rps)
			tdat <- rbind(tdat, tdat.h)
		}
		
		dat.all <- rbind(dat.all,tdat)
		write.csv(tdat, paste(x, "_pMetab-{dCV,dPV}-", trial.name, ".csv", sep=""), quote=F, row.names=F)
		setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress bar
	} # end loop over files
	close(pb)
	# convert all data from list to data frame
	cn <- colnames(dat.all)
	dat.all <- data.frame(matrix(unlist(dat.all), nrow=nrow(dat.all), ncol(dat.all)))
	colnames(dat.all) <- cn
	
	cat("\n")
	print(paste("mean of all rps = ",mean(dat.all[,"rps"])))
	
	directions <- c("dpv", "dcv")
	for (direction in directions) {
		## sort the whole data set based on direction
		dat <- dat.all[order(dat.all[,direction]),]
		dat.d <- data.frame()
		## loop through unique distances for each direction
		for (d in unique(dat[,direction])) {
			tmpdf <- dat[dat[direction] == d,]
			dat.d <- rbind(dat.d,c(d,nrow(tmpdf),mean(tmpdf[,"rps"])))
		}
		colnames(dat.d) <- c(direction,"Hamts","avgrps")
		
		#arps <- (dat.d[,"Hamts"]/sum(dat.d[,"Hamts"]))%*%(dat.d[,"avgrps"])
		#print(paste("average rps = ",arps))
		
		## plot points for all pMetab vs direction and a line for avg
		print(paste("Plotting avg(pMetab) and amt(vHPC) vs ",direction))
		outFile <- paste("graphics/",x,"_pMetab-",direction,"-plot",sep="")
		png(paste(outFile, ".png", sep=""), 1600, 1600)
		par(lwd=3, mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
		#points(dat[,"dcv"], dat[,"prs"], pch="·")
		plot(dat.all[,direction], dat.all[,"rps"], xlab=direction, ylab="reaction probability", pch="·", cex=10)
		lines(dat.d[,direction],dat.d[,"avgrps"],col="red",lwd=10)
		# make barplot of vHPC amounts
		outFile <- paste("graphics/",x,"_Hamts-rps-",direction,"-plot",sep="")
		png(paste(outFile, ".png", sep=""), 1600, 1600)
		par(mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
		#dtv <- as.vector(dt)
		barplot(dat.d[,"Hamts"], names.arg=dat.d[,direction], main=x, xlab=direction, ylab="Hamts")
		box()
	} # end loop over directions
} # end loop over experiments

quit()
