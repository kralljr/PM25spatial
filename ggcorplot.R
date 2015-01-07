library(ggplot2)

#define a helper function (borrowed from the "ez" package)
ezLev=function(x,new_order){
	for(i in rev(new_order)){
		x=relevel(x,ref=i)
	}
	return(x)
}

ggcorplot = function(data,var_text_size,cor_text_limits){
	z = data
	# names(z)=c('x','y','x_lab','y_lab')
	names(z)=c('cor','x_lab','y_lab', "type")
	
	z$x_lab = ezLev((z$x_lab), consallmatch[, 2])
	z$y_lab = ezLev((z$y_lab), consallmatch[, 2])
	z=z[z$x_lab!=z$y_lab,]

	f = facet_grid(y_lab~x_lab,scales='free_y')
    o = theme(
    plot.background = element_blank()
   ,panel.grid.major = element_blank()
   ,panel.grid.minor = element_blank()
   ,panel.border = element_blank()
   ,panel.background = element_blank()
  ) 
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


	# size_scale = scale_size(limits = c(0,1))
	return(
		ggplot(z, x = cor, group = type, colour = type)+
		xlab("Correlation") +
		ylab("Density") +
		# points_layer+
		geom_density(aes(x=cor, colour=type),fill=NA) +
		geom_vline(aes(yintercept = 0), colour = "grey80", 
		linetype = "dashed") +
		# xlim(-1, 1) +
		    # scale_colour_hue(drop = FALSE, name="",
                     # breaks=c("Trad.", "Spatial"),
                     # labels=c("Trad.", "Spatial"),
                     # values = cbPalette[-1]) + 
		scale_x_continuous(limits=c(-0.5, 1), breaks = c( 0, 0.5)) +
		f+
		o+
		# size_scale + 
		scale_colour_manual(values=cbPalette[-1], name = "") +
		geom_hline(aes(yintercept = 0), colour = "grey80", size = 0.7) 
	)

}