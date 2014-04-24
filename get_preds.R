


# # > rm(list = ls())
# > load("/Users/jennakrall/Dropbox/PM25cons_mort/PMcons_mort_analysis/datfinal_9aug11.RData")
# > load("/Users/jennakrall/Documents/Old/sourceapp/Prednaivedata_17aug11.RData")
# > 
# > load("/Users/jennakrall/Dropbox/PM25cons_mort/pm25cons_mort_st/spatial_model/items_pm25cons.RData")
# > monscons <- readRDS("/Users/jennakrall/Desktop/cache/monitor-subset.rds")
# > save(list = ls(), file = "all_info_forpred.RData")



getpreds <- function(fn1, PARsufcons, PARsufPM, 
	rdatapath, outpath, type = "reg") {
		
		
	#for running locally
	
	require(splancs)
	require(gpclib)
	require(splines)
	require(mvtnorm)
	require(tlnise)
	require(tsModel)
	require(maps)
	
	
	consall <- c("OC_K14" , "Elemental_Carbon", 
		"silicon", "Sodium_Ion", "SULFATE", "NITRATE", 
		"AMMONIUM", "PM25_Mass", "PM")
		
	cities <- yesmons[-which(yesmons %in% c("taco", "ftwa"))]  
	
	#get monitors for cities PMairs
	unmons <- unique(datfinal[, c(2, 4, 3)])
	colnames(unmons) <- c("mon", "x", "y")
	monsPM <- list()
	for(i in 1 : length(cities)) {
		coord <- citypolymat[[cities[i]]]
		colnames(coord) <- c("x", "y")
		monsPM[[i]] <- unmons[which(inout(unmons[, c(2, 3)], coord)), 1]
	}
	names(monsPM) <- cities	
	
	seeds1 <- c(6991,2094,9008,411,1434,5146,4293,4679, 542)
	for(i in 1 : 9) {

		seed <- i
		seedUSE <- seeds1[seed]
		cons <- consall[i]
		print(cons)

		#name output
		fn <- paste0(fn1, cons)

		#Create list per component of naive/predicted
		spatialavg <- list()
		avg <- list()
		
		
		set.seed(seedUSE)
		#Loop through all cities
		for(j in 1 : length(cities)){
			
			city <- cities[j]
			regions <- areawh[city]
			mons <- monsin[[city]]
		
			if(cons == "PM") {
				
				#get parameters
				dfile <- paste0(PARsufPM, regions, ".RData")
				load(file.path(rdatapath, dfile))
				pars <- out$opt
				pars <- pars$par

				#get data  for monitor
				data1 <- datfinal
				mons <- monsPM[[city]]
				
	
			}else{
	
				#get parameters
				dfile <- paste0(PARsufcons, regions, cons, ".RData")
				load(file.path(rdatapath, dfile))
				pars <- out$opt
				pars <- pars$par
				
				#get data for monitor
				data1 <- monscons
				

			
			}
			
			#get in correct format
			data <- data1[(data1[, 2] %in% mons), ]	
			data <- getalldata(data, consall[i])
			
			#get community averages
			avg1 <- avgfun(data, city)
					
			if(type == "reg") {
				type1 <- regions
			}else{
				type1 <- NULL
				}
			data2 <- getalldata(data1, consall[i], type1)			
					
			#give dates from avg, get predicted from spatial		
			spatialavg[[j]] <- spatialpred(data2, pars, avg1[, 1],
					city, 5000, cons)
		
			#Find observation number for matching dates
			x <- match(spatialavg[[j]][, 1], avg1)
			x <- x[!is.na(x)]
			avg[[j]] <- avg1[x, ]
			# }
		}#end loop over cities
		
			
		names(spatialavg)<- cities
		names(avg)<- cities
		
		save(spatialavg, avg, file = file.path(outpath, paste0(fn, ".RData")))
		
	}
		
}