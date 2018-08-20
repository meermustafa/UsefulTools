# Meer Mustafa
# 4.27.18


### R script for plotting functions

# colors ----

# R base colors
colors()
# http://research.stowers.org/mcm/efg/R/Color/Chart/


library(RColorBrewer)
RColorBrewer::display.brewer.all()

library(wesanderson)
wesanderson::wes_palettes

#http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
#install.packages('mapsdata')
library(mapsdata)
library(maps)
maps::
china = map_data()



# plot panelling ----
gridExtra::grid.arrange
ggdraw
# ggdraw(): Initialize an empty drawing canvas
# draw_plot(): Places a plot somewhere onto the drawing canvas.
# draw_plot_label(): Adds a plot label to the upper left corner of a graph. It can handle vecgors of labels with associated coordinates.

ggdraw() + 
  draw_plot(bxp, x = 0, y = 0.5, width = 0.5 , height = 0.5) +
  draw_plot(dp, x = 0.5, y = 0.5, width = 0.5, height = 0.5) + 
  draw_plot(lp, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A","B","C"), x= c(0, 0.5, 0), y = c(1,1,0.5), size = 15)




# pie chart ----
getPrettyPie = function(pieDF, labels, title) {
    ### takes in DF with V1 = vector of labels, V2 = counts
    pie(table( c(rep(pieDF[,1], times = pieDF[,2]) ) ),
    labels = pieDF[,1],
    main = title, # title
    cex.main = 2.8, # change title size
    #cex.axis = 1.2,
    radius = 1.2, # change size of pie
    col = sample(x = unique(colors()), size = nrow(pieDF), replace = F) # change colour of pie
       ) 

}



