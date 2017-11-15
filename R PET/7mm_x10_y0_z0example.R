setwd("/the/location/of/the/input/files/")

library(MASS)         # for kde2d
library(PET)
library(reshape2)     # melt
library(scales)       # for rescale
library(ggplot2)
library(RColorBrewer) # colors for heatmaps and plots

# Set colors
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(1024)
# Not a true heatmap as interpolated (kernel density estimation)

# GEOMETRY OF PET
# 1 layer with radius 437.3 mm
# 384 strips
# size of single strip: 7 mm x 19 mm x 500 mm 
# size of single strip: 4 mm x 19 mm x 500 mm 

N <- 384L   # No of detectors
R <- 43.73 # PET radius, cm
Delta.s <- R*pi/N # why doubled???

half.index <- as.integer(round(R/Delta.s)) # half of No. of samples by displacement
s.index <- as.integer(seq(-half.index,half.index,1))
# angle approximate
appr.displ <- s.index*Delta.s

######################### Import data ###########################################
# important! the axis X and Y are mixed up in the example (last two columns), 
# hence the correct location of the source will be at X=0 cm, Y=10 cm
lor.data <- read.table("PSF_384strips_370kBq_600s_L050_x10_y0_z0", header = FALSE)
lor.norm <- data.frame(round(lor.data[,10],3),round(lor.data[,11],3))
# add 2 points to create full scale of sinogram (-43.73,43.73)x(0,pi)
lor.norm <- rbind(lor.norm,
                  c(-43.73,0),
                  c(43.73,pi))
names(lor.norm) <- c("displacement","angle")
lor.norm <- lor.norm[order(lor.norm$angle,lor.norm$displacement),]
# rearrange
rownames(lor.norm) <- 1:nrow(lor.norm)

# check dimensions
# length.s <- length(unique(lor.norm$displacement))
# length.fi <- length(unique(lor.norm$angle))

########### TEMPORARY GRID ##############
#s.grid <- (43.73/250)*c(0,seq(1:250))
s.grid <- appr.displ[-(1:244)]
fi.grid <- (pi/384)*c(0,seq(1:384))

full.s.grid <- c(sort(-s.grid[-1]),s.grid)

###############################################################
# KDE algorithm {MASS}
# Default call 
# Adjust binning (interpolate - can be computationally intensive for large datasets)
# 245 points - displacement, 284 - angle
k <- kde2d(lor.norm$displacement, lor.norm$angle, 
           h = c(diff(appr.displ)[1]*2,
                 diff(fi.grid)[1]*2), 
           n = c(245,384))
# the best bandwidth h=c(diff(appr.displ)[1]*2,diff(fi.grid)[1]*2) - double size of sinogram pixel 

# sinogram image
image(k$x,k$y,k$z, xlim = c(0,15),col=r,xlab = "s, cm", ylab = expression(phi),
      main = "Original sinogram")

# ########################################################################
# perform iradon
# sinogram from Kernel density estimation (kde2d {MASS})
t <- Sys.time()
# 1000 x 1000 px - corresponds to ~4x zoom (interpolation used)
# iradon - inverse radon (requires sinogram to be transposed! Hence t(k$z) is an input)
irMy <- iradon(t(k$z),
                XSamples = 1000,
                YSamples = 1000,
                FilterTyp = "Ramp")
Sys.time() - t # how much time did the reconstruction take
xyMy <- seq(-max(k$x)/sqrt(2), max(k$x)/sqrt(2), length.out = 1000)

# plots of restored images
# ggplot2 rescale to dataframe again
df <- melt(irMy$irData)
df$Var1 <- rep(xyMy,1000)
df$Var2 <- rep(xyMy,each = 1000)
df$value <- rescale(df$value) # Intensity from 0 to 1
names(df)[3] <- "Intensity"
ggplot(df,aes(Var1,Var2,fill=Intensity)) + geom_tile() + 
  scale_fill_gradientn(colours=r) +
  xlim(-10,10) + ylim(-4,16) + ggtitle("FBP with KDE used") +
  labs(x="X, cm",y="Y, cm") 