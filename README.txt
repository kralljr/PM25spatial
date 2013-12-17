spattempcl.R
Functions for fitting spatial model for PM2.5 constituents
1. load libraries
2. getconsdata function: select constituent of interest from data
3. makeNLL function: function to create likelihood function
4. stmodlog function: organize data to then optimize likelihood

get_spatial_cov.R
1. eucdist: function to get Euclidean distance between one point and matrix of pts
2. avgfun: function to get average of concentrations in a community
3. unifs: function to get uniformly distributed points in a community
4. distfun: function to get distances between monitors in community
5. ndates: function to get unique dates with more than 3 monitors reporting
6. matmaternf: function to get matern covariance
7. wtfun: function to get logged constituent concentrations for each monitor
8. h11fun: function to get variance for monitors
9. h12fun: function to get covariance between unifpts and monitors
10. spatialpred: function to get predicted concentrations (per day) for a community 
