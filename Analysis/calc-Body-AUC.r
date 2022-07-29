#! /usr/bin/Rscript

###
## Read multiple Body *.csv files and calculate a pre-specified
## area-under-the-curve (AUC) for each column.
##
## Time-stamp: <2020-12-07 aks>
###

argv <- commandArgs(T)

usehalfpeak <- T

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: calc-Body-AUC.r <Body CSV file1> <Body CSV file2> ...")
  print("e.g. ./bin/calc-Body-AUC.r exp_body-avg.csv")
  quit()
}
if (length(argv) < 1) usage()

shift <- function(x, n){
  return(c(x[-(seq(n))],rep(x[length(x)], n)))
  # c(tail(x, -n), rep(NA, n))
}

diffadd <- function(x) {
	return(c(diff(x),0))
}

expfunc <- function(x, alpha, beta) {
	y <- alpha*exp(beta*x)
	return(y)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

datafiles <- argv
nfiles <- length(datafiles)

for (f in datafiles) {
	print(paste("Working on", f))

	if (!file.exists(f)) {
		print(paste(f,"doesn't exist."))
		next
	}

	## parse file name
	dnamef <- dirname(f)
	bnamef <- basename(f)
	xcnamef <- unlist(strsplit(bnamef,"_"))
	xnamef <- xcnamef[1]
	cnamef <- xcnamef[2]
	cnamef <- substr(cnamef, 0, regexpr('.csv', cnamef)-1)
	fileName.base <- paste(xnamef,cnamef,sep='_')

	## read data
	dat <- read.csv(f, check.names=FALSE, colClasses="numeric")
	dat.time <- dat[,1]
	
	## initialize the AUC data structures
	columns <- colnames(dat[,2:ncol(dat)])
    qAUC <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(qAUC) <- columns
    hAUC <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(hAUC) <- columns
    iAUC <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(iAUC) <- columns
    iAUC.inf <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(iAUC.inf) <- columns
    for (cn in columns) {
		# measure %in% c("nectrig","necrotic","stressed")
		if (cn %!in% c("Compound")) next
		print(cn)
		if (all(dat[,cn] == 0)) {
			print(paste("column ",cn,"contains all zeros...skipping"))
			next
		}
		
		## prepare the data used to calculate the AUC
		## find the peak or max value of the Body profile 
		peak <- max(dat[,cn])
		ipeak <- which.max(dat[,cn])
		## new data set as subset from old using index-at-peak
		ndat <- dat[ipeak:nrow(dat),]
		ndat.time <- ndat[,1]
		aucdat <- ndat
		aucdat.time <- ndat.time
		if (usehalfpeak) {
			halfpeak <- 0.5*peak
			if (ndat[nrow(ndat),cn] > halfpeak) {
				print(paste("last data value",ndat[nrow(ndat),cn]," is > halfpeak",halfpeak))
				print("skipping this AUC calculation")
				next
			}
			istart <- which.min(abs(ndat[,cn] - halfpeak))
			iend <- length(ndat[,cn])
			## further subset of data using index-at-half-peak from new data set
			aucdat <- ndat[istart:iend,]
			aucdat.time <- aucdat[,1]
		} 
		
		## calculate AUC by quadrature
		i <- 2:nrow(aucdat)
		qAUC[,cn] <- ((aucdat[i,cn] + aucdat[i-1,cn])/2) %*% (aucdat.time[i] - aucdat.time[i-1])
		#i <- (istart+1):iend
		#AUC[,cn] <- ((dat[i,cn] + dat[i-1,cn])/2) %*% (dat.time[i] - dat.time[i-1])
		#hAUC[,cn] <- t((halfdat[,cn] + halfdat[-1,cn])/2) %*% diff(halfdat.time)
		#hAUC[,cn] <- t((halfdat[,cn] + shift(halfdat[,cn],1))/2) %*% diffadd(halfdat.time)
		
		## prepare the data to be fit
		#dat.x <- halfdat[,1]
		#dat.y <- halfdat[,cn]
		dat.x <- aucdat[,1]
		dat.y <- aucdat[,cn]
		dat.df <- data.frame(x = dat.x, y = dat.y)
		plot(dat.x,dat.y)
		
		# Prepare a good inital state by fitting a linear model
		#theta.0 <- max(dat.y) * 1.1
		#model.0 <- lm(log(- y + theta.0) ~ x, data=dat.df)
		#a <- 1
		#b <- 1
		#model.0 <- lm(y ~ I(a*x+b), data=dat.df)
		model.0 <- lm(y ~ x, data=dat.df)
		#alpha.0 <- -exp(coef(model.0)[1])
		alpha.0 <- coef(model.0)[1]
		beta.0 <- coef(model.0)[2]
		#print(summary(model.0))
		lines(dat.df$x, predict(model.0, list(x = dat.df$x)), col = 'red', lwd = 3)
		
		#start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
		#start <- list(alpha = alpha.0, beta = beta.0)
		start <- list(alpha = alpha.0, beta = 0.000001)

		# Fit the exponential model
		#model <- nls(y ~ alpha * exp(beta * x) + theta , data = dat.df, start = start)
		model <- tryCatch(nls(y ~ alpha * exp(beta * x), data = dat.df, start = start),
		error = function(e) { print(paste("error for column ",cn," in nls exponential fit.")); print(e$message)})
		if (is.list(model)) {
			alpha <- coef(model)[1]
			beta <- coef(model)[2]
			#print(summary(model))
			#print(paste("alpha = ",alpha))
			#print(paste("beta = ",beta))
			lines(dat.df$x, predict(model, list(x = dat.df$x)), col = 'blue', lwd = 3)
			
			# integrate fitted expontential to calculate AUC
			#x <- dat.x
			result <- integrate(expfunc,alpha=alpha,beta=beta,lower=dat.x[1],upper=dat.x[length(dat.x)])
			inf.result <- tryCatch(integrate(expfunc,alpha=alpha,beta=beta,lower=dat.x[1],upper=Inf), 
			error = function(e) { print(paste("error for column ",cn," in integral with Inf upper limit.")); print(e$message)})
			#inf.result <- integrate(expfunc,alpha=alpha,beta=beta,lower=dat.x[1],upper=Inf)
		
			iAUC[,cn] <- result$value
			if (is.list(inf.result)) {
				iAUC.inf[,cn] <- inf.result$value
			}
		}
    }
	
	#m <- nls(y ~ I(a*exp(-b*x)+c), data=df, start=list(a=max(y), b=1, c=10), trace=T)
	#y_est<-predict(m,df$x)
	#exponential.model <- nls(val ~ a*exp(b*time), start=c(b=-0.1,h=30))
	
	#cat("AUC results by discrete quadrature\n")
    #print(qAUC)
    #cat("AUC results by adaptive quadrature using R integrate() function with definite limits\n")
    #print(iAUC)
    #cat("AUC results by adaptive quadrature using R integrate() function with Inf upper limit\n")
    #print(iAUC.inf)
    
    AUC.out <- rbind(qAUC,iAUC,iAUC.inf)
	rownames(AUC.out) <- c("qAUC","iAUC","iAUC.inf")
	if (usehalfpeak) {
		write.csv(AUC.out,paste(fileName.base,"-Body-AUC-hp.csv",sep=""),row.names=T)
	} else {
		write.csv(AUC.out,paste(fileName.base,"-Body-AUC-p.csv",sep=""),row.names=T)
	}
	
}

quit()
