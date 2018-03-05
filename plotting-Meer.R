

### R script for plotting functions

# pie table
pie(table( c(rep('>1 Mb',times = 437), rep('≤1 Mb',times = 1358-437) ) ),
    main = 'Non-overlapping Enhancers \nDistance to Nearest Oncogene', 
    cex.main =1.5, # change title size
    radius = 2, # change size of pie
    col = c('black','white')) # change colour of pie
