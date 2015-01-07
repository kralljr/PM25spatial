#####
# check data for dtmonscons
load("/Users/jennakrall/Dropbox/PM25cons_mort/pm25cons_mort_st/spatial_model/dtmonscons.RData")


getregion <- function(latlong) {
	lat <- latlong[2]
	long <- latlong[1]
	if(lat > 41.9  & long <= -109) {
		region <- "nw"
	}else if(lat > 40.9 & long > -109 & long <= -86.11) {
		region <- "nmw"
	}else if((lat > 36.45 & long > -86.11) | (lat > 36.45 & 
		lat <= 40.9 & long > -89.3 & long <= -86.11)) {
			region <- "ne"
	}else if(lat <= 41.9 & long <= -109) {
		region <- "sw"
	}else if(lat <= 40.9 & long > -109 & long <= -89.3) {
		region <- "smw"
	}else if(lat <= 36.45 & long > -89.3) {
		region <- "se"
	}else{
		stop("error")
	}
	
	region
}


tabs <- matrix(nrow = 8, ncol = 6)
for(i in 5 : 12) {
	dat <- dtmonscons[, c(1:4, i)]
	dat <- dat[complete.cases(dat), ]
	unmon <- unique(dat[, c(2, 4, 3)])
	regs <- vector()
	for(j in 1 : nrow(unmon)) {
		regs[j] <- getregion(unmon[j, c(2, 3)])
	}
	tabs[i - 4, ] <- (table(regs))
}
rownames(tabs) <- colnames(dtmonscons)[5:12]
colnames(tabs) <- sort(regions)
xtable(tabs)


# # 		if(region == "nw") {
			# datasubplus <- 
				# datasubplus[which(datasubplus[, 3] > 41.9 & 
				# datasubplus[, 4] <= -109),] 
		# }else if(region == "nmw") {
			# datasubplus <- 
				# datasubplus[which(datasubplus[, 3] > 40.9 & 			
				# datasubplus[, 4] > -109 & 
				# datasubplus[, 4] <= -86.11),]
		# }else if(region == "ne") {
			# datasubplus <- 
				# rbind(datasubplus[which(datasubplus[, 3] > 36.45 & 
				# datasubplus[, 4] > -86.11),],
		 		# datasubplus[which(datasubplus[, 3] > 36.45 & 			
		 		# datasubplus[, 3] <= 40.9 & 
		 		# datasubplus[, 4] > -89.3 & 
		 		# datasubplus[, 4] <= -86.11),])
		 # }else if(region == "sw") {
		 	# datasubplus <- 
		 		# datasubplus[which(datasubplus[, 3] <= 41.9 & 			
		 		# datasubplus[, 4] <= -109),]
		 # }else if(region == "smw") {
		 	# datasubplus <- 
		 		# datasubplus[which(datasubplus[, 3] <= 40.9 & 			
		 		# datasubplus[, 4] > -109 & 
		 		# datasubplus[, 4] <= -89.3),]
		 # }else	if(region == "se") {
		 	# datasubplus <- 
		 		# datasubplus[which(datasubplus[, 3] <= 36.45 & 			
		 		# datasubplus[, 4] > -89.3),]
		 	# }