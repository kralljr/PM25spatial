#####################################
# Functions to make ST/simple average AP predictions
#
# Edited 12/13/13
#
#
######################################


#Error when run 8/18
#Warning message:
#'initTLNise' is deprecated.
#Use ''seed' argument directly in 'tlnise'' instead.
#See help("Deprecated")
#



#Function to find distances between matrix of points 	#(pts) and one point (xy)
eucdist <- function(xy, pts) {
	pts <- as.matrix(pts, ncol = 2)
	xy <- as.matrix(xy)
	#Sweep difference over coordinates
	subs <- (sweep(pts, 2, xy, "-"))
	
	#Sum squared rows (differences), then take square root 
	subsum <- sqrt(rowSums(subs^2))
	subsum
	}
	






	
#Function to create simple average estimates
avgfun <- function(data, city) {
	
	#Get unique dates represented by the monitors
	data <- data[order(data[, 1]), ]
	date <- unique(data[, 1])
	
	#Apply mean over all monitors for each date, 		
	avgs <- tapply(data[, 5], data[, 1],
		mean, na.rm = TRUE)
	
	cbind(date, avgs)
}





#############
##Predicted##
#############

#Function to find uniform points within city
unifs <- function(N, city) {
	numcit <- which(names(bboxcor) == city)
	
	for (i in 1 : 2) {
		#Find points in city bounding box
		unifpts <- matrix(c(runif(N, 
			min = bboxcor[[numcit]][1],
			max = bboxcor[[numcit]][2]),
			runif(N, min = bboxcor[[numcit]][3],
			max = bboxcor[[numcit]][4])),
			byrow = FALSE, ncol = 2)
				
		#Which of these points are actually in the city
		unifpts <- unifpts[which(inout(unifpts,
			citypolymat[[numcit]]) == TRUE), ]
			
		#How many points are actually in the city
		n <- length(unifpts[, 1])
		
		#Redo with a larger number of points
		N <- N / n * N
		}
		unifpts
	}





########
#Function to find data in region
distfun <- function(rdata, cons) {

	#lat long of unique monitors
	monssub <- substr(rdata[, 2], 1, 9)
	
	if(cons == "PM") {
		monssub <- rdata[, 2]
	}
	rdata <- unique(data.frame(rdata[, c(4, 3)], monssub))
	rdata <- rdata[order(rdata[, 3]), ]

	monitors <- rdata[, c(1, 2)]
	monnames <- rdata[, 3]
	
	# distmat <- as.matrix(dist(monitors, method = "euclidean"))
	nr <- nrow(monitors)
	distmat <- matrix(nrow = nr, ncol = nr)
	for(i in 1 : nrow(monitors)) {
		distmat[i, ] <- distMeeus(monitors, monitors[i, ]) / 1000
		distmat[i, i] <- 0
	}
	
	
	colnames(distmat) <- monnames
	rownames(distmat) <- monnames
	list(distmat, monitors)
	}







	
	
	
	
#######	
#Find dates with >n observations in region
ndates <- function(dates, rdata, num = 3) {
	undate <- sort(dates)
	
	for(i in 1 : length(undate)) {
		x <- length(rdata[which(rdata[, 1]
			== undate[i]), 1])
		if (x < num) {
			undate[i] <- NA
		}
	}
	undate <- undate[!is.na(undate)]
	undate
}





	
#Function to calculate the matern function 
matmaternf <- function(mat, pars) {
	matrixmatern <- pars[1]^2 * matern(mat,
		pars[2], pars[3])
	matrixmatern
	}
	




########
#Function to get pollutant levels organized
#List by unique date of monitor ids and pollutant levels for
	#monitors with values on that day
wtfun <- function(undate, rdata){
	wt <- list()
	wtcomp <- list()
	
	#for each unique date
	for (i in 1 : length(undate)) {
		
			#Look at IDs and component levels (not logged)
			#monitors, levels, monitors short
			wt[[i]] <- rdata[which(rdata[, 1] == 
				undate[i]), c(2, 5, 6)]
			
			#Find mean of colocated monitors
			newlevels <- tapply(wt[[i]][, 2], wt[[i]][, 3], mean)
			
			#Eliminate zeroes before log
			for (j in 1 : length(newlevels)) {
				if(!is.na(newlevels[j])) {
					if(newlevels[j] == 0) {
						newlevels[j] <- 0.00001
					}
				}else{
					browser()
					}
			}
			
			wt[[i]] <- data.frame(as.character(unique(wt[[i]][, 3])),
				log(newlevels))
			wt[[i]] <- wt[[i]][complete.cases(wt[[i]]), ]
			
			#only constituent
			wtcomp[[i]] <- wt[[i]][, 2]
			
			#Name by ID
			names(wtcomp[[i]]) <- wt[[i]][, 1]
		} 
		
	#Name by unique dates
	names(wt) <- undate
	names(wtcomp) <- undate
	wtcomp
	}
	
	
	
	
	
	
	
########	
h11fun<-function(distmat, pars, unmonid, undate, wtcomp, type = "matern") {
	#Calculate matern function on distances 
	#between all monitors in region
	
	if(type == "matern") {
		corr <- matmaternf(distmat, pars) 
		n <- 5
	}else{
		corr <- pars[1]^2 * exp(-pars[2] * distmat)
		n <- 4
		}
	
	h11 <- corr + diag(x = 1, 
		nrow = length(distmat[, 1])) * pars[n]^2 

	#Name columns/rows by all monitor IDs
	colnames(h11) <- unmonid
	rownames(h11) <- unmonid
	
	
	#Order h11 by correct monitor IDs
	h11sub <- list()
	#for each date
	for(i in 1 : length(undate)) {
		h11sub[[i]] <- h11[names(wtcomp[[i]]), names(wtcomp[[i]])]
		}	
	h11sub
	}
	
	
	
	
########	
h12fun <- function(monitors, unifpts, pars, 
	unmonid, undate, wtcomp, type = "matern") {
		
		
	#Calculate distances/matern between random points and monitors
	materndist <- matrix(ncol = length(unifpts[, 1]),
		nrow = length(monitors[, 1]))
	for (i in 1 : length(monitors[, 1])) {
		# materndist[i, ] <- eucdist(monitors[i, ], unifpts) 
		materndist[i, ] <- distMeeus(monitors[i, ], unifpts)/1000
		materndist[i, i] <- 0
	}
	
	
	if(type == "matern") {
		h12mat <- matmaternf(materndist, pars) 
	}else{
		h12mat <- pars[1]^2 * exp(-pars[2] * materndist)
		}
	n <- length(unifpts[, 1])
	h12 <- 1 / n * rowSums(h12mat)

	#Name columns/rows by all monitor IDs
	#Order h12 by correct names
	names(h12) <- unmonid
	h12sub <- vector(mode = "list", )
	for (i in 1 : length(undate)) {
		h12sub[[i]] <- h12[names(wtcomp[[i]])]
	}
	h12sub
}	








#########
#Main Function
#########

spatialpred <- function(data, pars, dates, city, N,
	cons, num = 3, type = "matern") {
	#Find number corresponding to the city
	unifpts <- unifs(N, city)
	
	#Get data with logs for region from STmodel
	# not for regional
	# data <- out$datasubplus 
	
	#Creation of dist matrix for monitors (no jitter)
	#reorder and discard duplicates
	temp <- distfun(data, cons)
	distmat <- temp[[1]]
	monitors <- temp[[2]]

	#Monitors in the regions
	undate <- ndates(dates, data, num)
	unmonid <- rownames(distmat)
		
	wtcomp <- wtfun(undate, data)
	h11sub <- h11fun(distmat, pars, unmonid, undate, wtcomp, type)
	h12sub <- h12fun(monitors, unifpts, pars, 
		unmonid, undate, wtcomp, type)
		
	if(type == "matern") {
		n <- 4
	}else{
		n <- 3
		}
	
	#Calculate predicted values
	xgivw <- vector(, length = length(undate))
	for (i in 1 : length(undate)) {
		xgivw[i] <- pars[n] + h12sub[[i]] %*% solve(h11sub[[i]]) %*%	(wtcomp[[i]] - pars[n])
	}
	cbind(undate, exp(xgivw))
}








