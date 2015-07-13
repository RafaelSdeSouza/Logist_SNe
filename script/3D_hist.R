require(plot3D)
require(gdata)
require(Cairo)
data.1= read.table(file="..//data//multinomial_cat.dat",header=TRUE,na.strings = "",sep="\t")
data.1$Galtype2<-trim(data.1$Galtype2)

dat3D<-data.1[,c(1,4)]
VV<-as.matrix(table(dat3D))
labels <- c('Ia', 'Ib', 'Ibc','Ic', 'II', 'IIn')
labels2<-c("E", "E/S0", "Im","S","S0")


# transparent colors

CairoPDF(file="Hist3D.pdf",height = 8,width = 8)
hist3D(z = VV, scale = FALSE, expand = 0.00175, bty = "g",
       col = jet2.col(200, alpha = 0.5), border = "black",xlab="SNe",ylab="galaxy",
       zlab=list("Number of SNe",rot = -90),axes=FALSE,  space = 0.05,ticktype = "detailed",
       d = 1, cex.axis = 1e-10,tck = -.01,colkey=list(length=0.35,dist=-0.075,
                                                      tick=FALSE))
text3D(x = seq(0,1,0.2), y = rep(-0.25, 6), z = rep(0, 6),
       labels = labels,
       add = TRUE, adj = 0)
text3D(y = seq(0,1,0.25), x = rep(1.15, 5), z = rep(0, 5),
       labels = labels2,
       add = TRUE, adj = 0)

dev.off()



# save plotting parameters
pm <- par(mfrow = c(2, 2))
pmar <- par(mar = c(5.1, 4.1, 4.1, 2.1))

## =======================================================================
##  colorkey as argument of a plot3D function
## =======================================================================
# default, colkey = NULL: adds colkey because multiple colors
image2D(z = volcano)  

# default, colkey = NULL: no colkey because only one color
image2D(z = volcano, col = "grey", shade = 0.2, contour = TRUE)  

# colkey = FALSE: no color key, no extra space foreseen
image2D(z = volcano, colkey = FALSE)

# colkey = list(plot = FALSE): no color key, extra space foreseen
image2D(z = volcano, colkey = list(plot = FALSE, side = 3))
colkey (side = 3, add = TRUE, clim = range(volcano))


## =======================================================================
##  colorkey in new plot
## =======================================================================

colkey(side = 1, clim = c(0, 1), add = FALSE, clab = "z", 
       col.clab = "red", adj.clab = 0)
colkey(side = 2, clim = c(0, 1), clab = "z", length = 0.5, width = 0.5)
colkey(side = 3, clim = c(0, 1), lwd = 3, clab = c("a","b","c","d"), 
       line.clab = 5)
colkey(side = 4, clim = c(1e-6, 1), clog = TRUE, 
       clab = "a very long title in bold and close to the key", 
       line.clab = 1, side.clab = 2, font.clab = 2)
