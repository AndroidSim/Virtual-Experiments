#! /usr/bin/Rscript

###
## Cumulative statistics for the SS geometries and Lobule paths.
##
## Time-stamp: <2019-07-17 18:09:25 gepr>
###
argv <- commandArgs(T)

usage <- function() {
  print("Usage: ssgeom.r <exp directories>")
  print("  directories should contain files named ssgeom-[0-9]+.csv")
  quit()
}

library(stringr) # for prefix manipulation
library(matrixStats)

if (length(argv) <= 0) usage()

prefixes <- c("ssgeom")

for (prefix in prefixes) {
  cat(paste("Working on",prefix,"\n"))

  for (expDir in argv) {
    expName <- basename(expDir)
    expDir <- paste(dirname(expDir), expName, sep="/") # removes trailing slash from the argument

    cat(paste("\t", "experiment ", expDir,"\n",sep=""))
    if (!file.exists(expDir)) {
      print(paste(expDir," doesn't exist."))
      next
    }
    trials <- list.files(path = expDir, pattern=paste(prefix,"-[0-9]+.csv",sep=""), recursive = TRUE)
    if (length(trials) <= 0) {
      cat(paste("WARNING: No",prefix,"files in",expDir,"\n"))
      next
    }

    ## per trial stats
    trial.stats <- data.frame()
    trial.names <- c()
    trial.layers <- vector("list")
    for (trial in trials) {
      trial.name <- substr(basename(trial), 0, regexpr('.csv', basename(trial))-1)
      trial.number.s <- substr(trial.name, attr(regexpr('ssgeom-',trial.name), "match.length")+1, nchar(trial.name))
      trial.number.i <- as.numeric(trial.number.s)
      trial.names <- c(trial.names, paste("Trial", trial.number.i))
      ## SS stats
      ssdat <- read.csv(paste(expName,"/",trial,sep=""))

      ## min.circ, max.circ, μ.circ, σ.circ, min.length, max.length, μ.length, σ.length
      ssdat.m <- as.matrix(ssdat)
      trial.stats <-
        rbind(trial.stats,
              c(
                  colRanges(ssdat.m)[3,1:2], colMeans(ssdat.m)[3], colSds(ssdat.m)[3], # circ
                  colRanges(ssdat.m)[4,1:2], colMeans(ssdat.m)[4], colSds(ssdat.m)[4] # length
              ))

      ## layer stats
      min.layer <- colRanges(ssdat.m)[2,1]
      max.layer <- colRanges(ssdat.m)[2,2]
      layer.stats <- data.frame()
      for (layer in min.layer:max.layer) {
        ssdat.m.layer <- ssdat[which(ssdat$Layer == layer),]
        ssdat.m.layer.m <- as.matrix(ssdat.m.layer)
        layer.stats <-
          rbind(layer.stats,
                c(
                    colRanges(ssdat.m.layer.m)[3,1:2], colMeans(ssdat.m.layer.m)[3], colSds(ssdat.m.layer.m)[3], # circ
                    colRanges(ssdat.m.layer.m)[4,1:2], colMeans(ssdat.m.layer.m)[4], colSds(ssdat.m.layer.m)[4] # length
                ))
      }
      colnames(layer.stats) <- c("min(Circ)", "max(Circ)", "μ(Circ)", "σ(Circ)", "min(Length)", "max(Length)", "μ(Length)", "σ(Length)")
      rownames(layer.stats) <- seq(min.layer, max.layer)
      trial.layers[[trial.number.s]] <- layer.stats

    } # end for (trial in trials) {
    colnames(trial.stats) <- c("min(Circ)", "max(Circ)", "μ(Circ)", "σ(Circ)", "min(Length)", "max(Length)", "μ(Length)", "σ(Length)")
    trial.stats.m <- as.matrix(trial.stats)
    ## exp SS stats
    trial.stats <- rbind(trial.stats, colSums(trial.stats.m))
    trial.stats <- rbind(trial.stats, colMeans(trial.stats.m))
    trial.stats <- rbind(trial.stats, colSds(trial.stats.m))

    rownames(trial.stats) <- c(trial.names, "Σ", "μ", "σ")
    write.csv(trial.stats, paste(expName, "_ssstats.csv",sep=""))

    ## exp layer stats
    max.layers <- -1
    for (trial.layer in trial.layers) max.layers <- max(max.layers, nrow(trial.layer))
    layer.stats <- data.frame()
    layer.names <- c()
    for (layer in 0:(max.layers-1)) {
      layers.per.trial <- data.frame()
      layer.names <- c(layer.names, paste("Layer", layer))
      ## get this layer from each trial and insert into layer.stats
      for (trial.layer in trial.layers) {
##        print(paste(layer.names,"trial"))
        this.trial.this.layer <- trial.layer[rownames(trial.layer) == layer,]
##        print(this.trial.this.layer)
        layers.per.trial <- rbind(layers.per.trial, this.trial.this.layer)
      }
      layers.per.trial.m <- as.matrix(layers.per.trial)
      layer.stats <- rbind(layer.stats, colMeans(layers.per.trial.m))
##      print(layer.stats)
    }
    ## now add the stats of the stats of the stats
    layer.stats.m <- as.matrix(layer.stats)
    layer.stats <- rbind(layer.stats, colSums(layer.stats.m))
    layer.stats <- rbind(layer.stats, colMeans(layer.stats.m))
    layer.stats <- rbind(layer.stats, colSds(layer.stats.m))

    rownames(layer.stats) <- c(layer.names, "Σ", "μ", "σ")
    colnames(layer.stats) <- colnames(trial.stats) # same column names as trial.stats // maybe should add an extra () wrapper?
    write.csv(layer.stats, paste(expName, "_layerstats.csv", sep=""))

  } ## end for (expDir in argv)

  ## total stats

} ## for (prefix in prefixes)
