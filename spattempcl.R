##########################################
# File to fit spatial-temporal models
#
# Cleaned 12/13/13
# Revised 8/15/10 to average collocated monitors instead
# of jitter lat/long
# Revised 8/8/11 to work with pm25 monitors
#
# dependence on raw data
#
# Creates STmodellog.RData (too large for computer)
#
##########################################

#Load PM2.5 data
#monloc<-.readRDS("/home/bst/student/jkrall/airpollution/STmodelMUTAU/monitor_locations.rds")
#> dim(monloc)
#[1] 19417     3
#> head(monloc)
#  Latitude Longitude    monitor
#1 32.43746 -86.47289 01001.0001
#2 32.42833 -86.44361 01001.0002
#3 32.33266 -86.79152 01001.0003
#4  0.00000   0.00000 01003.0001
#5 30.55237 -87.70691 01003.0002
#6 30.55528 -87.71361 01003.0003


#leves<-.readRDS("/home/bst/student/jkrall/airpollution/STmodelMUTAU/monitor_locations.rds")
# length(leves)
#[1] 1685
#> head(leves[[1]])
#        Date PM25
#1 2000-01-01   NA
#2 2000-01-02   NA
#3 2000-01-03   NA
#4 2000-01-04   NA
#5 2000-01-05   NA
#6 2000-01-06   NA
#> head(names(leves))
#[1] "01003.0010" "01027.0001" "01033.1002" "01049.1003" "01053.0002"
#[6] "01055.0010"


#Load required packages
lloc <- "/home/bst/student/jkrall/Rlibrary"
library(splancs)
library(geoR)
library(splines)
library(mvtnorm)
library(geosphere, lib.loc = lloc)
library(matrixcalc, lib.loc = lloc)
library(timeDate, lib.loc = lloc)


getd <- function(dat) {
	whcc <- which(complete.cases(dat) & dat[, 3] > 0)
			
	dow <- (dayOfWeek(timeDate(dat[, 1])))
	
	dat <- data.frame(dat, dow)
	colnames(dat) <- c("date", "mon", "cons", "dow")
	dat <- dat[whcc, ]
	dat[, 3] <- log(dat[, 3])
	
	dates2 <- seq(as.Date("2000-01-01", origin = "1970-01-01"), 
		to = as.Date("2005-12-31", origin = "1970-01-01"), by = 1)
	dates2 <- data.frame(dates2)
		
	unmon <- unique(dat$mon)
	lm1 <- vector(, length = nrow(dat))
	for(j in 1 : length(unmon)) {
		
		wh1 <- which(dat$mon == unmon[j])
		dat1 <- dat[wh1, ]
		if(nrow(dat1) > 50) {
			dat1 <- merge(dat1, dates2, by.x = "date", by.y = "dates2", 
				all.y = T)
		
		
			temp <- try(lm(cons ~ factor(dow) + 
				ns(date, df = 7 * 6), data = dat1)$resid)
			if(class(temp) != "try-error") {

					lm1[wh1] <- temp
			}else{
				print(unmon[j])
				lm1[wh1] <- rep(NA, length(wh1))
				} 
		}
		
	}
	list(lm1, whcc)

}

detrenddat <- function(monscons, class = "cons") {
	
	
	if(class != "pm") {
		monsconsd <- monscons[, 1 : 12]
		for(i in 5 : 12) {
			dat <- monscons[, c(1, 2, i)]
			temp <- getd(dat)
			monsconsd[temp[[2]], i] <- temp[[1]]
		}
	}else{
		monsconsd <- monscons
		dat <- monscons[, c(1, 2, 5)]
		temp <- getd(dat)	
		monsconsd[temp[[2]], 5] <- temp[[1]]
		
		
		}

	monsconsd
}


getconsdata <- function(monscons, namecomponent) {

	data <- cbind(monscons[, c(1 : 4)], monscons[, namecomponent])
	data[, 6] <- substr(data[, 2], 1, 9)	
	data
}

#Function to maximize likelihood
makeNLL <- function(x, D, tau = T, type = "matern"){
	    function(p){
	    	#set parameters
			sigma <- p[1]
			phi <- p[2]
                        if(type == "matern") {
                             kappa <- p[3]
                             mu <- p[4]
                             n <- 5

                        }else{
                            mu <- p[3]
                            n <- 4
                        }
      	 	if(tau == T) {
      			tau <- p[n]
      		}else {
      			tau <- 0
      		}
		
		
			dmv <- vector(,length = length(x))
			S <- vector(mode="list",length=length(x))
			for(i in 1 : length(x)) {
				#Set matern covariance of distances
                                if(type == "matern") {
                                    spatcov <- matern(D[[i]], phi, kappa)
                                }else{

                                    spatcov <- exp(-phi * D[[i]])
                                }
				S[[i]] <- sigma^2 * spatcov + diag(x = 1,
					nrow = length(D[[i]][, 1])) * tau^2  
				#Find MVN lhood for mean mu and var S
				tpd <- try(is.positive.definite(S[[i]]))
                if(class(tpd) != "try-error"){
                	if(tpd) {
				dmv[i] <- dmvnorm(x[[i]],
					rep(mu, length(x[[i]])),
                                                  S[[i]], log = TRUE)
                       }else{dmv[i] <- -10^10}
                              }else{dmv[i] <- -10^10}
			}
		print(p)
		print(-sum(dmv))
                        out <- -sum(dmv)
#                        if(is.nan(-sum(dmv))){browser()}
                        if(is.nan(-sum(dmv))){
                          out <- 10^10
                        }
                           out
                      }
	}


		
		
		
getalldata <- function(data, namecomponent, region = NULL, detrend = T) {
	
	if(namecomponent != "PM2.5" & namecomponent != "PM") {
		datasub <- getconsdata(data, namecomponent)
	}else {
		datasub <- data
		datasub[, 6] <- as.character(data[, 2])
	}

	#Find complete cases where level>0
	datasub <- datasub[complete.cases(datasub), ]

	if(detrend == F) {
		   datasubplus <- datasub[which(datasub[, 5] > 0),]
	}	
	#Subset data based on region of interest
	if(!is.null(region)) {
		if(region == "nw") {
			datasubplus <- 
				datasubplus[which(datasubplus[, 3] > 41.9 & 
				datasubplus[, 4] <= -109),] 
		}else if(region == "nmw") {
			datasubplus <- 
				datasubplus[which(datasubplus[, 3] > 40.9 & 			
				datasubplus[, 4] > -109 & 
				datasubplus[, 4] <= -86.11),]
		}else if(region == "ne") {
			datasubplus <- 
				rbind(datasubplus[which(datasubplus[, 3] > 36.45 & 
				datasubplus[, 4] > -86.11),],
		 		datasubplus[which(datasubplus[, 3] > 36.45 & 			
		 		datasubplus[, 3] <= 40.9 & 
		 		datasubplus[, 4] > -89.3 & 
		 		datasubplus[, 4] <= -86.11),])
		 }else if(region == "sw") {
		 	datasubplus <- 
		 		datasubplus[which(datasubplus[, 3] <= 41.9 & 			
		 		datasubplus[, 4] <= -109),]
		 }else if(region == "smw") {
		 	datasubplus <- 
		 		datasubplus[which(datasubplus[, 3] <= 40.9 & 			
		 		datasubplus[, 4] > -109 & 
		 		datasubplus[, 4] <= -89.3),]
		 }else	if(region == "se") {
		 	datasubplus <- 
		 		datasubplus[which(datasubplus[, 3] <= 36.45 & 			
		 		datasubplus[, 4] > -89.3),]
		 	}
 	}
 	
 	datasubplus
	
}	



	
stmodlog <- function(data, namecomponent, region = NULL, n, s, u, 
	l = c(0.0001, 0.0001, 0.0001, -100, 0.0001),
	tau = T, detrend = F, type1 = "matern") {
		
		
	if(namecomponent != "PM2.5") {
		datasub <- getconsdata(data, namecomponent)
	}else {
		datasub <- data
		datasub[, 6] <- as.character(data[, 2])
	}

	#Find complete cases where level>0
	datasub <- datasub[complete.cases(datasub), ]
	datasubplus <- datasub[which(datasub[, 5] > 0),]

	#Subset data based on region of interest
	if(!is.null(region)) {
		if(region == "nw") {
			datasubplus <- 
				datasubplus[which(datasubplus[, 3] > 41.9 & 
				datasubplus[, 4] <= -109),] 
		}else if(region == "nmw") {
			datasubplus <- 
				datasubplus[which(datasubplus[, 3] > 40.9 & 			
				datasubplus[, 4] > -109 & 
				datasubplus[, 4] <= -86.11),]
		}else if(region == "ne") {
			datasubplus <- 
				rbind(datasubplus[which(datasubplus[, 3] > 36.45 & 
				datasubplus[, 4] > -86.11),],
		 		datasubplus[which(datasubplus[, 3] > 36.45 & 			
		 		datasubplus[, 3] <= 40.9 & 
		 		datasubplus[, 4] > -89.3 & 
		 		datasubplus[, 4] <= -86.11),])
		 }else if(region == "sw") {
		 	datasubplus <- 
		 		datasubplus[which(datasubplus[, 3] <= 41.9 & 			
		 		datasubplus[, 4] <= -109),]
		 }else if(region == "smw") {
		 	datasubplus <- 
		 		datasubplus[which(datasubplus[, 3] <= 40.9 & 			
		 		datasubplus[, 4] > -109 & 
		 		datasubplus[, 4] <= -89.3),]
		 }else	if(region == "se") {
		 	datasubplus <- 
		 		datasubplus[which(datasubplus[, 3] <= 36.45 & 			
		 		datasubplus[, 4] > -89.3),]
		 	}
 	}
	
	#Find unique dates
	dates <- sort(unique(datasubplus[, 1]))
	
	#Organize component data by date
	wtvec <- vector(mode = "list", length = length(dates))
	len <- vector(, length = length(dates))
	for(i in 1 : length(dates)) {
		#Which data matches with date i
		wtvec[[i]] <- datasubplus[which(datasubplus[, 1] == 
			dates[i]),c(2, 5, 6)]
		#Number of observations associated with date i 
		len[i] <- length(wtvec[[i]][, 1])
		
		
		#	*AVG*
		#The following coincides with code in complexfun419.R
		#Find mean of colocated monitors                                      
   		newlevels <- tapply(wtvec[[i]][, 2],
   			wtvec[[i]][, 3], mean)
   		
   		#remove zeroes before logging	
   		for(j in 1 : length(newlevels)) {
			if(newlevels[j] == 0){
				newlevels[j] <- NA
				}
			}
   		
   		#wtvec is list, each element is matrix of 
   			#monid (no collocated) and logged levels *AVG*
		if(detrend == F) {
	   		wtvec[[i]] <- cbind(as.numeric(names(newlevels)),
	   			log(newlevels))
	   	}else{
	   		wtvec[[i]] <- cbind(as.numeric(names(newlevels)),
	   			(newlevels))
	   		
	   		}

   		wtvec[[i]] <- wtvec[[i]][which(complete.cases(wtvec[[i]]) == TRUE), ]
		}
	
	#Limit data to days with at least n observations
	wtvecn <- wtvec[which(len >= n)]
	
	
	#Reduce wtvecn to only logged component data
	wt <- vector(mode = "list", length = length(wtvecn))
	for(i in 1 : length(wtvecn)) {
		wt[[i]] <- wtvecn[[i]][, 2]
		}



	#Mon ids, latitude, longitude
	moninfo <- cbind(as.numeric(datasubplus[, 6]),
		datasubplus[, 4], datasubplus[, 3])	

	#Unique monitor days
	uniquemon <- unique(moninfo)

	#cons	
	#Latitude and longitude of unique monitors
	moncoord <- matrix(c(as.numeric(uniquemon[, 2 : 3])),
		byrow = FALSE, ncol = 2)
		

	#Find distances between all monitors
	# distmat <- as.matrix(dist(moncoord, method = "euclidean"))
	#better distance, in kM
	 nr <- nrow(moncoord)
	 distmat <- matrix(nrow = nr, ncol = nr)
	 for(i in 1 : nrow(moncoord)) {
		 distmat[i, ] <- distMeeus(moncoord, moncoord[i, ]) / 1000
		 distmat[i, i] <- 0
	 }
	#distmat <- as.matrix(dist(moncoord))
	#Name columns/rows with monitor ids
	colnames(distmat) <- uniquemon[, 1]
	rownames(distmat) <- uniquemon[, 1]


	
	#For each day, find distmat for monitors that report
	matchdist <- vector(mode="list", length = length(wtvecn))
	for(i in 1 : length(wtvecn)) {
		wmon<-wtvecn[[i]][, 1]
		class(wmon)<-"character"
		dm <- distmat[wmon, wmon]
		matchdist[[i]] <- dm
		}
		
	#wt is observations, matchdist is matrix of distances 
		#for day i
	#Maximize likelihood
	nll <- makeNLL(wt, matchdist, tau = tau, type = type1)
        ps <- c(1, 10, 1, 1, 1)
        if(type1 != "matern") {
            l <- l[-3]
            u <- u[-3]
            ps <- ps[-3]
            ps[2] <- 0.01
            s <- s[-3]
            l[2] <- 0.001
        }

        if(type1 == "matern" & namecomponent == "PM2.5" & is.null(region)) {
            ps[2] <- 1
            print(ps)

        }

        if(type1 == "matern" & namecomponent == "silicon" & region == "nmw") {
            #u[3] <- 100
            ps[2] <- 1
        }
        
	#Optimize likelihood
	opt <- optim(s, nll, method = "L-BFGS-B",
		lower = l,
		upper = u, hessian = TRUE, control = list(maxit = 500, parscale = ps))

        if(opt$conv != 0) {
            ps[2] <- 1
            if(type != "matern") {
                ps[2] <- 0.001
            }
            opt <- optim(s, nll, method = "L-BFGS-B",
                      lower = l,
                      upper = u, hessian = TRUE, control = list(maxit = 500, parscale = ps))
            if(opt$conv != 0){
                print("error: did not converge")
                print(c(opt$conv, opt$message))
            }
        }

        
	listall <- list(datasubplus, opt)
	names(listall) <- c("datasubplus", "opt")
	listall
	
	
	# list(wt, matchdist, distmat, wtvecn)
}


# save.image(file="STmodellog.RData")
