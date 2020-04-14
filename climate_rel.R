#library(stargazer)

# opening tif files for annual mean temperature and annual precipitation
  # annual mean temperature from WorldClimd database, code BIO01
temp <- raster("./DATA/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
res(temp) # res as [1] 0.1666667 0.1666667
# changing res to match that of cardamom and butler (1,1)
temp.1 <- aggregate(temp, fact = 1/res(temp))
res(temp.1) # now is 1,1
plot(temp.1)
png("./figures/")
temp.plot <- levelplot(temp.1, par.settings=BuRdTheme, margin=F,
                       par.strip.text=list(font=2, cex=1.2))
temp.plot

  # annual precipitation from WorldClimd database, code BIO12
ppt <- raster("./DATA/wc2.1_10m_bio/wc2.1_10m_bio_12.tif")
res(ppt) # res as [1] 0.1666667 0.1666667
# changing res to match that of cardamom and butler (1,1)
ppt.1 <- aggregate(ppt, fact = 1/res(ppt))
res(ppt.1) # now is 1,1
plot(ppt.1) 
bluescale <- colorRampPalette(c("#FFF8DC", "#7FCDBB",
                                "#40B6C4", "#2C7FB8" , "#253494"))
myTheme <- rasterTheme(region= dichromat(bluescale(10), "deutan"))
ppt.plot <- levelplot(ppt.1, par.settings=myTheme, margin=F)
ppt.plot

# saving panelled plots of climate data
png("./figures/panel.climate.png", width = 40, height = 30, units = "cm",
    res = 500)
(panel.climate <- ggarrange(ppt.plot, temp.plot, nrow = 2, ncol = 1))
text.temp <- grid::grid.text('Mean annual temperature (°C)', rot=90,
                             y=unit(0.25, "npc"), 
                             x=unit(0.90, "npc"))
text.ppt <- grid::grid.text('Annual precipitation (mm)', rot=90,
                            y=unit(0.8, "npc"), 
                            x=unit(0.90, "npc"))
dev.off()

# turn the rasters to data frames, renaming column of values and removing NAs
temp.df <- as.data.frame(temp.1, xy=T) %>%
  rename(temp.C = wc2.1_10m_bio_1) %>%
  filter(temp.C!=0)
ppt.df <- as.data.frame(ppt.1, xy=T) %>% 
  rename(ppt.mm = wc2.1_10m_bio_12) %>%
  filter(ppt.mm!=0)

# joining dataframes with cardamom and butler 
# joining ppt and temp
climate <- left_join(ppt.df,temp.df)
# joining climate data frame with cardamom and butler 
# sla mean
cl.sla <- left_join(climate,joined_sla) %>%
  filter(cardamom!=0, butler!=0) %>% # filtered out NAs when joined with clim data
  mutate(`Cardamom`=cardamom, `Butler`=butler) %>%
  gather(key=dataset, value=sla.mean,-x,-y,-ppt.mm,-temp.C,-cardamom,-butler) 
# sla stdev
cl.std <- left_join(climate,joined.std.nafree) %>%
  filter(cardamom_std !=0, butler_std!=0) %>%# filtered out NAs when joined with clim 
  mutate(`Cardamom`=cardamom_std, `Butler`=butler_std) %>%
  gather(key=dataset, value=sla.stdev,-x,-y,-ppt.mm,-temp.C,
         -cardamom_std, -butler_std)



# density graphs overlapped 
(density.mean <- ggplot(cl.sla,aes(x=ppt.mm,y=sla.mean,fill=dataset)) + 
  stat_density2d(geom="tile", aes(fill = dataset, alpha=..density..), 
                 contour=FALSE) + 
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values=c("#0072B2", "#D55E00"))+
  theme_minimal())
 ggsave("./figures/density.gl.mean.png",density.mean,width = 40,height = 30,
       units = "cm",dpi = 500)


density.mean.df <- ggplot_build(density.mean)$data[[1]]
cb.mean <- merge(density.mean.df,cl.sla, by.x=c("x","y"),by.y=c("ppt.mm","cardamom"))
 

annotation <- data.frame(
  x = 1000,
  y = 20,
  label = "54% overlap"
)
(contour.mean<-ggplot(cl.sla,aes(x=ppt.mm,y=sla.mean,color=dataset)) + 
  stat_density2d(geom="density2d", aes(color = dataset,alpha=..level..),
                 size=2,
                 contour=TRUE) + 
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  labs(color=" ",x="\n Annual precipitation (mm)",
       y= "SLA mean (m2.kg-1)\n")+
  theme_minimal()+
  annotate("rect", xmin = 0, xmax = 1150, ymin = 5, ymax = 20,
             alpha = .2) +
  geom_label(data=annotation, aes( x=x, y=y, label=label),                 , 
             color="black", 
             size=4,fontface="bold"))
ggsave("./figures/contour.gl.mean.png",contour.mean,width = 30,
       height = 15,units = "cm",dpi = 500)

density.mean.data <- as.data.table(ggplot_build(density.mean)$data[1])
butler.density <- density.mean.data[group==1]
cardamom.density <- density.mean.data[group==2]
butler.density <- butler.density %>% dplyr::select(density) %>%
  rename(butler=density)
cardamom.density <- cardamom.density %>% dplyr::select(density) %>%
  rename(cardamom=density)
density <- cbind(butler.density,cardamom.density)
density <- list(cardamom = cardamom.density$cardamom,
              butler = butler.density$butler) 
density.overlap <- my.overlap.sla(density,plot = F) 
density.overlap$OV
#cardamom-butler 
# 0.5378 = 54% overlap 

library(plotly)
library(MASS)

# Compute kde2d
kd <- with(cl.sla, MASS::kde2d(ppt.mm, cardamom, n = 50))
kd.b <- with(cl.sla,MASS::kde2d(ppt.mm,butler, n=50))

# Plot with plotly
(plot3d <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% 
    add_trace(data = kd.b, x=kd.b$x,y=kd.b$y,z=kd.b$z)%>%
    add_surface())


library(reshape2) # For melt function
library(scales)
# Calculate the common x and y range 
xrng = range(c(cl.sla$ppt.mm, cl.sla$ppt.mm))
yrng = range(c(cl.sla$cardamom, cl.sla$butler))

# Calculate the 2d density estimate over the common range
d1 = kde2d(cl.sla$ppt.mm, cl.sla$cardamom, lims=c(xrng, yrng), n=200)
d1[["lon"]] <- cl.sla$x
d2 = kde2d(cl.sla$ppt.mm, cl.sla$butler, lims=c(xrng, yrng), n=200)

# Confirm that the grid points for each density estimate are identical
identical(d1$x, d2$x) # TRUE
identical(d1$y, d2$y) # TRUE

# Calculate the difference between the 2d density estimates
diff12 = d1 
diff12$z = d1$z - d2$z

## Melt data into long format
# First, add row and column names (x and y grid values) to the z-value matrix
rownames(diff12$z) = diff12$x
colnames(diff12$z) = diff12$y

# Now melt it to long format
diff12.m = melt(diff12$z, id.var=rownames(diff12))
names(diff12.m) = c("ppt.mm","sla.mean","Density")

library(scales)
# Plot difference between two densities
ggplot(diff12.m, aes(ppt.mm, sla.mean, z=Density, fill=Density)) +
  geom_tile() +
  stat_contour(aes(colour=..level..), binwidth=0.001) +
  scale_fill_gradient2(low="blue",mid="white", high="red", midpoint=0) +
  scale_colour_gradient2(low=muted("blue"), mid="white", high=muted("red"), midpoint=0) +
  coord_cartesian(xlim=xrng, ylim=yrng) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  guides(colour=FALSE)

### package colorplaner ####
devtools::install_github("wmurphyrd/colorplaner")
library(colorplaner)
col_func <- function(x, y){
  x[x == 0] <- 0.000001
  y[y == 0] <- 0.000001
  x[x == 1] <- 0.999999
  y[y == 1] <- 0.999999
  # upper or lower triangle?
  u <- y > x
  # Change me for different hues.
  hue <- ifelse(u, 0.3, 0.8)
  # distace from (0,0) to (x,y)
  hyp <- sqrt(x^2 + y^2) 
  # Angle between x axis and line to our point
  theta <- asin(y / hyp)
  # Angle between 45 degree line and (x,y)
  phi <- ifelse(u, theta - pi/4, pi/4 - theta)
  phi <- ifelse(phi < 0, 0, phi)
  # Distance from 45 degree line and (x,y)
  s <- hyp * sin(phi) / sqrt(2)
  # Draw line from (x, y) to 45 degree line that is at right angles.
  # How far along 45 degree line, does that line join.
  v <- 1 - hyp * cos(phi) / sqrt(2)
  # Get hsv values.
  sapply(seq_along(x), function(i) hsv(hue[i], s[i], v[i]))
}

(cardamom.ppt.tile <- ggplot(data=cl.sla,aes(x,y,fill=ppt.mm,fill2=cardamom))+
  geom_tile()+
    scale_fill_colourplane(name = " ",
                           na.color = "white",
                           color_projection = "interpolate",
                           vertical_color = "#006400",
                           horizontal_color = "#68228B", 
                           zero_color = "#E8E6F2",
                           axis_title = "Annual precipitation (mm)",
                           axis_title_y = "SLA mean (m2.kg-1)",
                           breaks = c(1000,3000,5000),
                           limits_y = c(0,63))+
 theme(panel.background = element_rect(fill =  "white"),
        panel.grid = element_blank())+
  ggtitle("Cardamom dataset")) # caption function doesnt work
#c("#F0FFFF", "#FFF8DC", "#E0FFFF")
ggtern::ggsave("./figures/global_analysis/cardamom.ppt.tile.png",cardamom.ppt.tile,
               width = 30,height = 15,units = "cm",dpi = 500)
# interesting point to remember, with the colour gradient used for the map, 
# the working function for saving the plot is ggsave from ggtern package,
# NOT from ggplot2 

(butler.ppt.tile <- ggplot(data = cl.sla,aes(x,y,fill=ppt.mm,fill2=butler))+
    geom_tile()+
    scale_fill_colourplane(name = " ",
                          na.color = "white",
                          color_projection = "interpolate",
                          vertical_color = "#006400",
                          horizontal_color = "#68228B", 
                          zero_color = "#E8E6F2",
                          axis_title = "Annual precipitation (mm)",
                          axis_title_y = "SLA mean (m2.kg-1)",
                          breaks = c(1000,3000,5000),
                          limits_y = c(0,63))+
    theme(panel.background = element_rect(fill =  "white"),
          panel.grid = element_blank())+
    ggtitle("Butler dataset"))
ggtern::ggsave("./figures/global_analysis/butler.ppt.tile.png",butler.ppt.tile,
               width = 30,height = 15,units = "cm",dpi = 500)

cardamom.ppt.tile.df <- ggplot_build(cardamom.ppt.tile)$data[[1]] %>% 
dplyr::select(fill,fill2,x,y)
butler.ppt.tile.df <- ggplot_build(butler.ppt.tile)$data[[1]] %>%
  dplyr::select(fill,fill2,x,y)%>%
  rename("butler.fill"=fill,"butler.fill2"=fill2)
combine.tiles <- left_join(cardamom.ppt.tile.df,
                           butler.ppt.tile.df)

# find matching colours 
combine.tiles <- combine.tiles %>%
  mutate(match.col = fill==butler.fill)

unique.col <- unique(c(combine.tiles$fill, combine.tiles$butler.fill))
combine.tiles <- combine.tiles %>%
  mutate(fill.n=as.numeric(factor(combine.tiles$fill, levels=unique.col)),
         fill2.n = as.numeric(factor(combine.tiles$butler.fill,levels = unique.col))) 
max(combine.tiles$fill.n)
max(combine.tiles$fill2.n)
(comb.tile <- ggplot(combine.tiles, aes(x,y,fill=fill.n,fill2=fill2.n))+
    geom_tile() +
    scale_fill_colourplane(name = " ",
                           na.color = "white",
                           color_projection = col_func,
                           limits = c(0,2196),
                           limits_y = c(0,2196),
                           axis_title = "cardamom",
                           axis_title_y = "butler",
                           breaks = c(500,1500),
                           breaks_y = c(500,1500)) +
    theme(panel.background = element_rect(fill =  "white"),
          panel.grid = element_blank())+
    ggtitle("Comparison of butler against cardamom for results from relationship 
            with precipitation (mm)"))
ggsave("./figures/global_analysis/comp.ppt.mean.png",comb.tile,
       width = 30,height = 15,units = "cm",dpi = 500)





