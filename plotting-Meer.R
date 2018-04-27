

### R script for plotting functions

# pie table
getPrettyPie = function(pieDF, title) {
    ### takes in DF with V1 = vector of labels, V2 = counts
    pie(table( c(rep(pieDF[,1],times = pieDF[,2]) ) ),
    main = title, 
    cex.main =1.5, # change title size
    radius = 2, # change size of pie
    col = unique(colors())[1:nrow(pieDF)] # change colour of pie
       ) 

}
