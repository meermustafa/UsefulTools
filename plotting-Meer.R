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




# ggplot ----

library(ggplot2)
library(ggrepel)

ggplot(df, aes(x = x, 
                   y = y
                   #, colour = SampleName
                   # ,label = UniqueID # use this or ggrepel
                   # ,group = UniqueID # group controls between what factorized group to draw lines between
                   )) +  
  
  geom_line( aes(color = Label), size = 1 ) + 
  
  geom_point( 
    alpha = 1/50,
    aes(color = SampleName, size = ValidationScore
      #,alpha = 1/20
    )
  ) + 
  
  # when labeling data points on a line plot - label them at the end of the line with added spacing 
  geom_dl(aes(label = UniqueID), method = list(dl.combine(
    "first.points", # choose at which points to label
    "last.points"), cex = 1, # label which points, what size
    dl.trans(x = x + 0.2), "last.points" ) # space between label and data point
  ) +

  # label data points by repelling to avoid overlap
  geom_label_repel(aes(label =  ifelse(Rank < 6 | Rank >= (allrankstotake-4) , 
                                       yes = as.character(Oncogenes),
                                       no = '') # label based on data
                      ,fontface = 'italic' # spacial characters for the label letters
                      ,color = 'indianred3' # color the label letters
                      )
                     ,check_overlap = T
                     ,hjust = 0
                     ,vjust=-.2
                     ,size = 4
                     ,font = 3
                     ,direction = "y"
                     ,nudge_x = -0.05
                     ,segment.size = 0.5
                     ,segment.color = "grey50"
                    , parse = TRUE # for italics comprehension
                    , element_text(face = 'italics')
                    ) +

  # pick x label ticks and limits
  scale_x_continuous(
    breaks=c(128000000,128500000,129000000,129500000,130000000,130500000) # tick makrs
    ,labels=c("128000000","128500000","129000000",'129500000','130000000','130500000') # labels
    ,limits = c(127500000, 131000000) # limits
  ) +
  
  
  # color manual via hex code
  scale_color_manual(values = c('#FF0000', # red
                                "#008000", # green
                                '#800080', # violet
                                "#0000FF", # blue
                                '#000000' # black
  )
  ) +
  
  theme_classic() +
  
  
  labs (title = paste0()
        ,subtitle = paste0()
        ,x = ""
        ,y = ""
        ,color = 'Legend Color Title'
        ,size = 'Legend Size Title'
  ) +
  
  # change legend data size and color transparency, despite the conditions in the plot
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))
         ,size = guide_legend(override.aes = list(alpha = 1))
  )+
  
  theme(title = element_text(size=22), 
        legend.title = element_text(size = 18), 
        legend.text = element_text(size=20), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_blank()
        ,axis.title.x = element_text(size=20, face="bold")
        ,axis.title.y = element_blank()
  )




# fancy R base plots ----
plot(x = xvalues, 
     y = yvalues,
     main = paste0(), 
     ylab = paste0(), 
     xlab =  paste0(),
     pch = 19,
     type = 'p',
     xlim = c(), 
     ylim = c(),
     cex= 0.7, cex.lab = 1.6, cex.axis = 1.5, cex.main = 1.8, col = alpha(cols[fileIndex],0.4)
)

# add line between points
lines(x = df$sgRNA_start[df$sgRNA_start > leftBoundary & df$sgRNA_start < rightBoundary], 
      y = yvalues, 
      type = 'l', 
      lwd = 0.3)

# text for a line
text(x = xmarker, y = 1.1, labels = '', cex = 1.5)
abline(v = xmarker, lty = 1, lwd = 0.5)

# filling in area of R plots ----
xx = c(xvalues, 
       rev(xvalues)
       )
yy = c(rep(0, length(yvalues)), rev(yvalues))
polygon(xx, yy, col = alpha(cols[fileIndex],0.3), border = NA)

abline(h=1, lty = 2, lwd = 0.5)




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



