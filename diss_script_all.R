###############################################################################
#                 BSc Hons Dissertation                                       #
#                 Assessing degree of consistency between                     # 
#                 functional leaf trait estimates from diff methods           #
###############################################################################

## things to do 
  
# MODEL = BUTLER 
# OBSERVATIONS = CARDAMOM = REFERENCE

## STIPPLING!!!!!!!! to represent uncertainty (units same as that of mean)

# scatter plot - make heatmap --| ASK CODING CLUB / MADE A HEATSCATTER (BETTER REPR) / ASK IF I COULD DO IT WITH THE MAP OF THE WORLD?

# for results (to put in table and have the figures showing)
# matrix of map, see how they differ from each other by each grid point 
# calculation of bias / RMSE / R2 - visual and tabular repr.

## Packages ----
install.packages("ncdf4")
install.packages("ggplot2")
install.packages("raster")
install.packages("SimDesign")
install.packages("Metrics")
install.packages("gplots")
install.packages("diffeR")
install.packages("PerformanceAnalytics")
install.packages("hrbrthemes")
install.packages("viridisLite")
install.packages("LSD")
install.packages("rworldmap")
install.packages("Hmisc")
install.packages("formattable")
install.packages("sf")
install.packages("magick")
install.packages("webshot")
webshot::install_phantomjs()
install.packages("overlapping")

library(ncdf4)
library(RColorBrewer)
library(ggplot2)
library(raster)
library(tidyverse)
library(tidyr)
library(SimDesign)
library(Metrics)
library(gplots)
library(reshape2)
library(purrr)
library(grid)
library(diffeR)
library(PerformanceAnalytics)
library(viridis)
library(hrbrthemes)
library(lattice)
library(data.table)
library(ggpubr)
library(LSD)
library(MASS) # do i need it?
library(rworldmap)
library(Hmisc)
library(sp)
library(ggExtra)
library(formattable)
library(sf)
library(kableExtra)
library(rasterVis)
library(overlapping)


### opening netcdf, to data frame ----
# opening the .nc files 
cardamom_sla <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc", 
                       varname="sla")
cardamom_75th <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc", 
                        varname = "75th_percentile")
cardamom_25th <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc",
                        varname = "25th_percentile")
cardamom_95th <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc",
                       varname = "95th_percentile")
cardamom_sla_std <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc", 
                           varname="Standard_Deviation")
butler_sla <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", 
                     varname="sla")
butler_sla_std <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", 
                         varname="sla_std")

# making .nc files into data frames 
cardamom_sla_df <- raster::as.data.frame(cardamom_sla, xy = TRUE)
cardamom_sla_std_df <- raster::as.data.frame(cardamom_sla_std, xy=TRUE)
butler_sla_df <- raster::as.data.frame(butler_sla, xy = TRUE) 
butler_sla_std_df <- raster::as.data.frame(butler_sla_std, xy=TRUE)

# basic data manipulation
cardamom_sla_df <- cardamom_sla_df %>%
  rename("cardamom" = sla)
cardamom_sla_std_df <- cardamom_sla_std_df %>%
  rename("cardamom_std" = Standard_Deviation)
butler_sla_df <- butler_sla_df %>%
  rename("butler" = specific.leaf.area)
butler_sla_std_df <- butler_sla_std_df %>%
  rename("butler_std" = specific.leaf.area)
joined_sla <- left_join(cardamom_sla_df, butler_sla_df) 
joined_sla_std <- left_join(cardamom_sla_std_df, butler_sla_std_df)

### visualising sla from both datasets SAME SCALE ----

  # SLA
png("./figures/plot_sla_SAMESCALE.png", width = 55, height = 35, 
    units = "cm", res = 200)
#par(mfrow=c(1,2), oma = c(0,3,8,0) + 0.1, mar = c(7,0,2,8) + 0.1)
par(mfrow=c(2,2), oma = c(0,2,4,2) + 0.1)
plot(cardamom_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")), 
     zlim=c(0,65), main="Cardamom\n", 
     ylab="Latitude")
plot(butler_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")), 
     zlim=c(0,65),  main="Butler\n",
     legend.args=list(text='\n Specific Leaf Area (m2.kg-1)', side=4, 
                      font=1, line=2.3)) # doesnt show ylab for some reason...

  # SLA STD
#png("./figures/plot_slastd_SAMESCALE.png", width = 50, height = 20,
 #   units="cm", res = 300)
#par(mfrow=c(1,2), oma = c(0,2,4,2) + 0.1)

setMinMax(butler_sla_std[[1]]) #values     : 2.167324, 17.66437  (min, max)
plot(cardamom_sla_std[[1]], asp=NA, col=rev(brewer.pal(10,"RdBu")), 
     zlim=c(0,75), xlab="\nLongitude", 
     ylab="Latitude")
plot(butler_sla_std[[1]], asp=NA, col=rev(brewer.pal(10,"RdBu")),
     zlim=c(0,75), xlab="\nLongitude",
     legend.args=list(text="\nSpecific Leaf Area Standard Deviation (m2.kg-1)",
                      side=4, font=1, line=2.3))
dev.off()

### visualising diff between CARDAMOM AND BUTLER ----
breakpoints <- c(-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45, 50,
                     55,60,65)
colors <- c("#27408B", "#36648B", "#4876FF", "#8DEEEE","#000000",
            "#FFEBCD", "#FFEC8B", "#EEDC82", "#FFD700", "#EEB422", "#CD9B1D", 
            "#8B6914", "#EE6363","#CD5555",
            "#CD2626",  "#B22222", "#8B1A1A","#000000")
png("./figures/plot_diff-CARD-BUTL.png", width = 50, height = 20, units = "cm", 
    res = 200)
par(mfrow=c(1,2), oma = c(0,2,4,2) + 0.1)
plot(cardamom_sla[[1]]-butler_sla[[1]], asp=NA, breaks=breakpoints, 
     col=colors, xlab="\nLongitude", ylab="Latitude", 
     legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                        side=4, font=1, line=2.3),
     main="Cardamom-Butler SLA Mean\n")

# could re-do it finding out a way to automatise this with package like rcolourbrewer?
plot(cardamom_sla_std[[1]]-butler_sla_std[[1]], asp=NA, breaks=breakpoints,
     col=colors, xlab="\nLongitude", 
     legend.args = list(text="\n\nSLA StDev (m2.kg-1)", 
                        side=4, font=1, line=2.3),
     main="Cardamom-Butler SLA StDev\n")
dev.off()



### STIPPLING FOR BUTLER SLA MEAN (25th-75th and 25th-95th percentiles) ----
# sla mean
gl_stip_locs <- (butler_sla[[1]] >= cardamom_25th[[1]])*
  (butler_sla[[1]]<=cardamom_75th[[1]])
gl_stip_locs <- rasterToPoints(gl_stip_locs)
gl_stip_locs <- gl_stip_locs[gl_stip_locs[, "layer"] == 1,]
gl_stip_locs_95 <- (butler_sla[[1]]>=cardamom_25th[[1]])*
  (butler_sla[[1]]<=cardamom_95th[[1]])
gl_stip_locs_95 <- rasterToPoints(gl_stip_locs_95)
gl_stip_locs_95 <- gl_stip_locs_95[gl_stip_locs_95[, "layer"] ==1,]
diff_pc <- as.data.frame(gl_stip_locs_95) %>% 
  setdiff(as.data.frame(gl_stip_locs)) %>% 
  as.matrix
png("./figures/stippling_world.png", width = 40, height = 25, 
    units = "cm", res = 500)
plot(butler_sla, asp = NA, col = rev(brewer.pal(10, "RdYlBu")),
     xlab="\nLongitude", ylab="Latitude", 
     legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                        side=4, font=1, line=2.3),
     main="Butler Sla Mean (stippling)\n")
points(stip_locs, pch = 18, cex=0.5)
points(diff_pc, pch = 23, cex = 0.7, col = "darkgreen", bg="green")
dev.off()

# sla stdev 
stip_locs_std <- (butler_sla_std[[1]] >= cardamom_25th[[1]])*
  (butler_sla_std[[1]]<=cardamom_75th[[1]])
stip_locs_std <- rasterToPoints(stip_locs_std)
stip_locs_std <- stip_locs_std[stip_locs_std[, "layer"] == 1,]
stip_locs_std_95 <- (butler_sla_std[[1]]>=cardamom_25th[[1]])*
  (butler_sla_std[[1]]<=cardamom_95th[[1]])
stip_locs_std_95 <- rasterToPoints(stip_locs_std_95)
stip_locs_std_95 <- stip_locs_std_95[stip_locs_std_95[, "layer"] ==1,]
diff_pc_std <- as.data.frame(stip_locs_std_95) %>% 
  setdiff(as.data.frame(stip_locs_std)) %>% 
  as.matrix # no difference - the points lay exclusively within 25-75pc here (for stdev)
png("./figures/stippling_std_world.png", width = 40, height = 25, 
    units = "cm", res = 500)
plot(butler_sla_std, asp = NA, col = rev(brewer.pal(10, "RdYlBu")),
     xlab="\nLongitude", ylab="Latitude", 
     legend.args = list(text="\n\nSLA StDev (m2.kg-1)", 
                        side=4, font=1, line=2.3),
     main="Butler Sla StDev (stippling)\n")
points(stip_locs_std, pch = 18, cex=0.5)
dev.off()


### 1) DATA EXPLORATION 1 - SERIES OF PLOTS ----

  # scatterplots (correlation chart) ----
png("./figures/corr_chart_sla.png", width = 50, height = 30, units = "cm", 
    res = 200)
corr_char_sla <- chart.Correlation(joined_sla_nocoord, histogram=TRUE, pch=19)
dev.off()
png("./figures/corr_chart_slastd.png", width = 50, height = 30, units = "cm", 
    res = 200)
corr_char_sla_std <- chart.Correlation(joined_sla_std_nocoord, histogram = TRUE, pch=19)
dev.off()

# both too many data points to be sure of what's going on - too clustered 

  # 2D density chart (correlation charts) ----

# raster function 
(heatmap_raster<- ggplot(joined_sla_nocoord, aes(cardamom, butler))+
    stat_density_2d(aes(fill=..density..), geom = "raster", contour = FALSE)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)))
ggsave("./figures/2d_density.png", last_plot(), width = 30, height = 20, units="cm", 
       dpi = 300)

(heatmap_polygon <- ggplot(joined_sla_noNA, aes(cardamom, butler))+
    stat_density_2d(aes(fill=..level..), geom = "polygon")+
    theme_classic()) 
# interesting that here the values for cardamom (spanning beyond 20 m2.kg-1) is not being represented - maybe the plot automatically only chooses those points which are shared closer between the two variables?
ggsave("./figures/2d_density_polygon.png", last_plot(), width = 30, height = 20,
       units = "cm", dpi = 300)

  # hexbin version of 2d map ----
  
    # SLA MEAN
(hexbin_sla_corr <- ggplot(joined_sla, aes(x=cardamom, y=butler) ) +
  geom_hex(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_classic()+
  xlab("\nCardamom")+
  ylab("Butler")+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("Mean SLA (m2.kg-1)\n"))
ggsave("./figures/hexbin_map_SLAcount200.png", hexbin_map_corr, width = 30, 
       height = 20, units = "cm", dpi = 300)

    # SLA STDEV
(hexbin_slastd_corr <- ggplot(joined_sla_std, aes(x=cardamom_std, y=butler_std) ) +
    geom_hex(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    theme_classic()+
    xlab("\nCardamom")+ 
    theme(plot.title = element_text(face="bold"))+
   ylab(" ") + ggtitle("StDev SLA (m2.kg-1)\n"))
ggsave("./figures/hexbin_map_SLASTDcount200.png", hexbin_slastd_corr, width = 30, 
       height = 20, units = "cm", dpi = 300)

(hexbin_panelled <- ggarrange(hexbin_sla_corr, hexbin_slastd_corr, ncol = 2))
ggsave("./figures/hexbin_map_panel.png", hexbin_panelled, width = 50, 
       height = 20, units = "cm", dpi = 300)

  # heatmap (with lattice package and viridis package for colour palette) ASK CODING CLUB! ----

  # modifying the dataset (joined_sla) to wide format for the matrix

joined_wide <- acast(joined_sla, x+y~cardamom, value.var = "butler")
  # this gives back a huge matrix - 6,4gb - cannot be processed to be visualised...

levelplot(joined_matrix, col.regions = terrain.colors(100)) # try cm.colors() or terrain.colors()

joined_matrix <- as.matrix(joined_sla_nocoord)
similarity.matrix <- apply(joined_matrix, 2, function(x)rowSums(x==joined_matrix))
diag(similarity.matrix)<-0

#### HEATSCATTERS ####
  # creating numeric vectors of sla mean and sla stdev to input in heatscatter ----
cardamom_sla_num <- cardamom_sla_df$sla
butler_sla_num <- butler_sla_df$specific.leaf.area
cardamom_sla_std_num <- cardamom_sla_std_df$Standard_Deviation
butler_sla_std_num <- butler_sla_std_df$specific.leaf.area

  # LSD package ----
png("./figures/heatscatter_sla_panelled.png", width = 50, height = 25, 
    units = "cm", res = 300)
par(mfrow=c(1,2), oma = c(0,2,4,2) + 0.1)
  # SLA MEAN
(heatscatter_sla <- heatscatter(cardamom_sla_num, butler_sla_num, pch = 19, 
                                cexplot = 0.5, colpal="spectral", 
                           #disco() for all color options / could set to colorblind
                                add.contour=TRUE, main = "SLA Mean\n",
                                xlab="\nCardamom", 
                                ylab="\nButler"))

  # SLA STDEV
(heatscatter_slastd <- heatscatter(cardamom_sla_std_num, butler_sla_std_num, 
                                  pch = 19, cexplot = 0.5, colpal = "spectral",
                                  add.contour = TRUE, main = "SLA StDev\n", 
                                  xlab="\nCardamom",
                                  ylab=" "))
dev.off()

  # heatscatter with ggplot2 - NOT FINISHED TRYING (SHOULD I ADD A VERTICAL LINE ACROSS IT?) ----

d <- na.omit(data.frame(x=cardamom_sla_num, y=butler_sla_num))
d <- kde2d(d$x, d$y)
d_df <- data.frame(with(d, expand.grid(x,y)), as.vector(d$z))
names(d_df) <- c("d_x","d_y","d_z")
model_fit <-loess(d_z~d_x*d_y, data = d_df)
d$pointdens <- predict(model_fit, newdata=data.frame(d_x=d$x, d_y=d$y))
(heatscatter_ggplot <- ggplot(d, aes(x,y,color=pointdens))+
    geom_point()+
    theme_ipsum())
                         

  # HISTOGRAM AND DIFF BETWEEN HISTS ---- need to do the diagonal bars where the bins overlap
  # SLA MEAN ----
joined_sla_density <- joined_sla %>%
  gather(key="dataset", value="sla",-y,-x)

(sla_hist <- ggplot(joined_sla_density,aes(x=sla, group=dataset,
                                               fill=dataset)) +
  geom_histogram(bins=100,alpha=0.4)+
  theme_ipsum()+
  scale_fill_discrete(name = "Specific Leaf Area", 
                      labels = c("Cardamom", "Butler")))

  # SLA STDEV ----

joined_slastd_density <- joined_sla_std %>%
  gather(key="dataset", value="sla_std",-x,-y)

(slastd_hist <- ggplot(joined_slastd_density, aes(x=sla_std,group=dataset,
                                                  fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="Specific Leaf Area StDev", 
                        labels=c("Cardamom", "Butler")))


#### PERCENTAGE OVERLAP DENSITY PLOTS ####
# sla mean ----
joined_sla <- joined_sla %>%
  filter(sla!=0, specific.leaf.area!=0)
sla_n <- list(cardamom = joined_sla$sla,
       butler = joined_sla$specific.leaf.area) 
png("./figures/global_density_overlap.png", width = 40, height = 20,
    units = "cm", res = 400)
sla_overl<- my.overlap(sla_n, plot = TRUE)
dev.off()
    
# sla stdev ----
joined_sla_std <- joined_sla_std %>%
  filter(Standard_Deviation!=0, specific.leaf.area!=0)
slastd_n <- list(cardamom_std = joined_sla_std$Standard_Deviation,
                 butler_std = joined_sla_std$specific.leaf.area)
png("./figures/global_std_density_overlap.png", width = 40,
    height = 20, units = "cm", res = 400)
sla_std_overl <- my.overlap(slastd_n, plot = T)
dev.off()


## category Components plot ----

categoryComponentsPlot(cardamom_sla, butler_sla, units = "m^2/kg")

## Cross-tabulate two RasterLayer objects, or mulitiple layers in a RasterStack ----
# or RasterBrick to create a contingency table.
crosstab <- crosstab(cardamom_sla,butler_nc)
crosstab <- as.data.frame(crosstab)
heatmap(crosstab)

cardamom_sla[30,30]
crosstabm <- crosstabm(cardamom_sla, butler_sla, percent = TRUE)
heatmap(crosstabm)

sla_only <- joined_sla %>%
  dplyr::select(-x,-y) 
sla_only <- as.matrix(sla_only)
 
sla_only<- rcorr(sla_only, type="pearson")



## CORRELATION CHART BETWEEN SLA PARAMETERS  ----

png("correlation_chart.png", width=30, height=20, units = "cm", res = 300)
dev.off()

new_joined <- new_joined %>%
  rename("cardamom"=sla, "butler"=specific.leaf.area)
  
png("corr_with_coords_chart.png", width=30, height=20, units = "cm", res = 300)
corr_char_coords <- chart.Correlation(new_joined, histogram=TRUE, pch=19)
dev.off()
#The distribution of each variable is shown on the diagonal.
#On the bottom of the diagonal : the bivariate scatter plots with a fitted line are displayed
#On the top of the diagonal : the value of the correlation plus the significance level as stars
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)





### GLOBAL DATA STATS: r2, rmse, bias ----
# SLA MEAN
global_sla_stat <- joined_sla %>%
  filter(cardamom!=0, butler!=0) %>%
  mutate(mean_c = mean(cardamom),
         mean_b = mean(butler),
         diff_butler = butler-mean_b,
         diff_butler2 = diff_butler^2,
         sum_diff_butler2 = sum(diff_butler2),
         slope_bf = sum((cardamom-mean_c)*(butler-mean_b))/
           sum((cardamom-mean_c)^2),
         b_intercept = mean_b - (slope_bf*mean_c),
         new_b_val = b_intercept + (slope_bf*cardamom),
         dist_mean_new_b = new_b_val - mean_b,
         sqrd_dist_b = dist_mean_new_b^2,
         sum_sqrd_dist_b = sum(sqrd_dist_b),
         sla_r2 = sum_sqrd_dist_b / sum_diff_butler2,
         rmse_av = rmse(butler, cardamom),
         rmse_row = sqrt(se(butler, cardamom)),
         bias = bias(butler, cardamom),
         bias_row = butler-cardamom)

# plotting rmse sla mean 
(sla_rmse_plot <- ggplot(global_sla_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))
(sla_rmse_marginal <- ggMarginal(sla_rmse_plot,type = "density",
                                 color="darkred", size = 5))

# SLA STD

global_slastd_stat <- joined_sla_std %>%
  filter(cardamom_std!=0, butler_std!=0) %>%
  mutate(meanstd_b = mean(butler_std),
         meanstd_c = mean(cardamom_std),
         diff_b_std = butler_std-meanstd_b,
         sqrd_diff_b = diff_b_std^2,
         sum_sqrd_diff_b = sum(sqrd_diff_b),
         slopestd_bf = (sum((cardamom_std-meanstd_c)*(butler_std-meanstd_b))/
                          sum((cardamom_std-meanstd_c)^2)),
         bstd_intercept = meanstd_b-(slopestd_bf*meanstd_c),
         new_bstd_val = bstd_intercept + (slopestd_bf*cardamom_std),
         dist_std_new_b = new_bstd_val - meanstd_b,
         sqrd_dist_std_b = dist_std_new_b^2,
         sum_sqrd_dist_std_b = sum(sqrd_dist_std_b),
         sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b,
         rmse_av = rmse(butler_std, cardamom_std),
         rmse_row = sqrt(se(butler_std, cardamom_std)),
         bias = bias(butler_std, cardamom_std),
         bias_row = butler_std-cardamom_std)

# plotting rmse sla stdev  
(slastd_rmse_plot <- ggplot(global_slastd_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat="identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab(" ")+
    xlab("\nLongitude")+
    labs(color="RMSE (m2.kg-1)")+
    ggtitle("SLA StDev RMSE\n")+
    theme(plot.title= element_text(face = "bold")))
(slastd_rmse_marginal <- ggMarginal(slastd_rmse_plot,type = "density",
                                    color="darkred", size = 5, margins = "y"))

(rmse_panelled <- ggarrange(sla_rmse_plot, slastd_rmse_plot, ncol = 2))
ggsave("./figures/RMSE_panel.png", rmse_panelled, width = 50, height = 20, 
       units="cm", dpi = 300)

## the R2 when im doing a linear regression - NOT BEING USED ATM (all commented)----
#lm_sla <- lm(butler ~ cardamom, data = joined_sla)
#plot(lm_sla)
#summary(lm_sla)$r.squared
# Multiple R-squared:  0.0003106,	Adjusted R-squared:  0.0002266 

#lm_slastd <- lm(butler_std  ~ cardamom_std, data = joined_sla_std)
#plot(lm_slastd)
#summary(lm_slastd) 
# Multiple R-squared: ,	Adjusted R-squared:  0.0002127 




### should i be doing t-test and f-test for sla mean and std respectively? could do 




################################################################################
#                           splitting of the world                             #
#               -- for these also need to calculate the bias!!! ---            #
################################################################################

###############################
#### SPLITTING BY LATITUDE ####
###############################
  # tropics ----
# cardamom:
  # sla mean
trpmat <- matrix(data <- c(-180,180,-23.5,23.5), nrow = 2, ncol = 2,
                 byrow = TRUE)
trpext <- extent(trpmat)
trp_c <- crop(cardamom_sla, trpext)
#plot(trp[[1]]) # tropical lats
tropicsSLA_df <- raster::as.data.frame(trp_c, xy=TRUE)
  # sla std
trpstd<- crop(cardamom_sla_std, trpext)
#plot(trpstd[[1]], asp=NA) # height to fix when (and if) saving it as png+
tropicsSTD_df <- raster::as.data.frame(trpstd, xy=TRUE)

# butler: 
  # sla mean
trp_b <- crop(butler_sla, trpext)
trp_b_df <- raster::as.data.frame(trp_b, xy=TRUE)
  # sla std
trp_b_std <- crop(butler_sla_std, trpext)
trp_b_std_df <- raster::as.data.frame(trp_b_std, xy=TRUE)

# joined tropics cardamom + butler
trpSLA <- left_join(tropicsSLA_df,trp_b_df)
trpSTD <- left_join(tropicsSTD_df, trp_b_std_df)
  # subtropics ----
# cardamom:
  # sla mean 
sbtrpmatN <- matrix(data <- c(-180,23.5,180,35), nrow = 2, ncol = 2)
sbtrpextN <- extent(sbtrpmatN)  
sbtrpN <- crop(cardamom_sla, sbtrpextN)
#plot(sbtrpN[[1]])
sbtrpN_df <- raster::as.data.frame(sbtrpN, xy = TRUE)
sbtrpmatS <- matrix(data <- c(-180,-23.5,180,-35), nrow = 2, ncol = 2)
sbtrpextS <- extent(sbtrpmatS)
sbtrpS <- crop(cardamom_sla, sbtrpextS)
#plot(sbtrpS[[1]])
sbtrpS_df <- raster::as.data.frame(sbtrpS, xy = TRUE)

# merging two datasets from two hemispheres to have them in one dataframe
sbtrp_c <- merge(sbtrpN_df, sbtrpS_df, by=c("x", "y", "sla"), all=TRUE) 

  # sla std 
sbtrp_std_N <- crop(cardamom_sla_std, sbtrpextN)
sbtrpSTD_c_N <- raster::as.data.frame(sbtrp_std_N, xy=TRUE)
sbtrp_std_S <- crop(cardamom_sla_std, sbtrpextS)
sbtrpSTD_c_S <- raster::as.data.frame(sbtrp_std_S, xy=TRUE)
sbtrpSTD_c <- merge(sbtrpSTD_c_N, sbtrpSTD_c_S, 
                    by=c("x","y","Standard_Deviation"), all = TRUE)

# butler:
  # sla mean
sbtrpN_b <- crop(butler_sla, sbtrpextN)  
#plot(sbtrpN_b[[1]])
sbtrpS_b <- crop(butler_sla, sbtrpextS)
#plot(sbtrpS_b[[1]])
sbtrp_bN_df <- raster::as.data.frame(sbtrpN_b, xy=TRUE)
sbtrp_bS_df <- raster::as.data.frame(sbtrpS_b, xy =TRUE)
sbtrp_b <- merge(sbtrp_bN_df, sbtrp_bS_df, by=c("x","y", "specific.leaf.area"),
                 all = TRUE)

  # sla std
sbtrpSTD_b_N <- crop(butler_sla_std, sbtrpextN)
sbtrpSTD_b_S <- crop(butler_sla_std, sbtrpextS)
sbtrpSTD_b_N_df <- raster::as.data.frame(sbtrpSTD_b_N, xy = TRUE)
sbtrpSTD_b_S_df <- raster::as.data.frame(sbtrpSTD_b_S, xy =TRUE)
sbtrpSTD_b <- merge(sbtrpSTD_b_N_df, sbtrpSTD_b_S_df, 
                    by = c("x","y","specific.leaf.area"),
                    all = TRUE)

# subtropics joined cardamom + butler
sbtrp_joined_SLA <- left_join(sbtrp_c, sbtrp_b)
sbtrp_joined_STD <- merge(sbtrpSTD_c, sbtrpSTD_b, by = c("x","y"), all=TRUE)

  # temperate ----
#cardamom:
  # sla mean
tmpmatN <- matrix(data <- c(-180,180,35,66.5), nrow = 2, ncol = 2, 
                 byrow = TRUE)
tmpextN<- extent(tmpmatN)
tmpmatS <- matrix(data <- c(-180, 180, -35,-66.5), nrow = 2, ncol = 2, 
                  byrow=TRUE)
tmpextS <- extent(tmpmatS)
tmpN <- crop(cardamom_sla, tmpextN)
tmpS <- crop(cardamom_sla, tmpextS)
#plot(tmpN[[1]])
#plot(tmpS[[1]])

tmp_c_N_df <- raster::as.data.frame(tmpN, xy=TRUE)
tmp_c_S_df <- raster::as.data.frame(tmpS, xy = TRUE)
tmp_c_df <- merge(tmp_c_N_df, tmp_c_S_df, by = c("x","y","sla"), all = TRUE)

  # sla std 
tmpSTD_c_N <- crop(cardamom_sla_std, tmpextN)
tmpSTD_c_S <- crop(cardamom_sla_std, tmpextS)
tmpSTD_c_N_df <- raster::as.data.frame(tmpSTD_c_N, xy =TRUE)
tmpSTD_c_S_df <- raster::as.data.frame(tmpSTD_c_S, xy =TRUE)
tmpSTD_c_df <- merge(tmpSTD_c_N_df, tmpSTD_c_S_df, 
                     by = c("x","y","Standard_Deviation"), 
                     all = TRUE)
# bulter 
  # sla mean 
tmp_b_N <- crop(butler_sla, tmpextN)
tmp_b_S <- crop(butler_sla, tmpextS)
tmp_b_N_df <- raster::as.data.frame(tmp_b_N, xy = TRUE)
tmp_b_S_df <- raster::as.data.frame(tmp_b_S, xy = TRUE)
tmp_b_df <- merge(tmp_b_N_df, tmp_b_S_df, by = c("x","y","specific.leaf.area"),
                  all = TRUE)

  # sla std 
tmpSTD_b_N <- crop(butler_sla_std, tmpextN)
tmpSTD_b_S <- crop(butler_sla_std, tmpextS)
tmpSTD_b_N_df <- raster::as.data.frame(tmpSTD_b_N, xy = TRUE)
tmpSTD_b_S_df <- raster::as.data.frame(tmpSTD_b_S, xy = TRUE)
tmpSTD_b_df <- merge(tmpSTD_b_N_df, tmpSTD_b_S_df, 
                     by = c("x","y","specific.leaf.area"), all = TRUE)

## joining cardamom and butler for sla mean and stdev of temperate regions
tmp_sla <- left_join(tmp_c_df,tmp_b_df)
tmp_slastd <- left_join(tmpSTD_c_df, tmpSTD_b_df)
  
  # poles ----
  # only data available in N hemisphere
# cardamom
# sla mean
plN <- matrix(data <- c(-180,66.5,180,90), nrow = 2, ncol=2)
plextN <- extent(plN)
plN_c <- crop(cardamom_sla, plextN)
#plot(plN_c[[1]])
plN_c_df <- raster::as.data.frame(plN_c, xy= TRUE)

# sla stdev
plN_slastd_c <- crop(cardamom_sla_std, plextN)
#plot(plN_slastd_c[[1]])
plN_slastd_c_df <- raster::as.data.frame(plN_slastd_c, xy = TRUE)

# butler
# sla mean
plN_b <- crop(butler_sla, plextN)
#plot(plN_b[[1]])
plN_b_df <- raster::as.data.frame(plN_b, xy = TRUE)

# sla stdev
plN_slastd_b <- crop(butler_sla_std, plextN)
#plot(plN_slastd_b[[1]])
plN_slastd_b_df <- raster::as.data.frame(plN_slastd_b, xy = TRUE)

# joining butler and cardamom for poles lat
plN <- left_join(plN_c_df, plN_b_df)
plN_std <- left_join(plN_slastd_c_df, plN_slastd_b_df)


#### CREATING THE FUNCTIONS FOR OVERLAP PLOTS (modif from overlapping package) ####
my.final.plot <- function (x, OV = NULL){
  AREA <- NULL
  for (i1 in 1:(length(x) - 1)) {
    for (i2 in (i1 + 1):(length(x))) {
      A <- data.frame(x = x[[i1]], group = names(x)[i1], 
                      k = paste(names(x)[i1], names(x)[i2], sep = "-", 
                                collapse = ""))
      B <- data.frame(x = x[[i2]], group = names(x)[i2], 
                      k = paste(names(x)[i1], names(x)[i2], sep = "-", 
                                collapse = ""))
      AREA <- rbind(AREA, rbind(A, B))
    }
  }
  if (!is.null(OV)) {
    OV <- data.frame(OV = OV, k = names(OV))
    AREA <- merge(AREA, OV, by = "k")
    AREA$k <- paste0(AREA$k, " (ov. perc. ", round(AREA$OV * 
                                                     100), ")")
  }
  ggplot(AREA, aes(x = x)) + facet_wrap(~k) + 
    geom_density(aes(fill = AREA$group), 
                 alpha = 0.35) + xlab("") + 
    theme(legend.title = element_blank())+
    theme_classic()+
    scale_color_brewer(palette = "Set1")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))
}

my.overlap <- function (x, nbins = 1024, plot = FALSE, partial.plot = FALSE, 
                        boundaries = NULL, ...){
  if (is.null(names(x))) 
    names(x) <- paste("Y", 1:length(x), sep = "")
  dd <- OV <- FUNC <- DD <- xpoints <- COMPTITLE <- NULL
  for (j in 1:length(x)) {
    if (!is.null(boundaries)) {
      Lbound <- lapply(boundaries, FUN = length)
      if ((Lbound$from == 1) & (Lbound$to == 1)) {
        warning("Boundaries were set all equals")
        boundaries$from <- rep(boundaries$from, length(x))
        boundaries$to <- rep(boundaries$to, length(x))
      }
      else {
        if ((Lbound$from != length(x)) | (Lbound$to != 
                                          length(x))) {
          stop("Boundaries not correctly defined")
        }
      }
      from = boundaries$from[j]
      to = boundaries$to[j]
      dj <- density(x[[j]], n = nbins, from = from, to = to, 
                    ...)
    }
    else {
      dj <- density(x[[j]], n = nbins, ...)
    }
    ddd <- data.frame(x = dj$x, y = dj$y, j = names(x)[j])
    FUNC <- c(FUNC, list(with(ddd, approxfun(x, y))))
    dd <- rbind(dd, ddd)
  }
  for (i1 in 1:(length(x) - 1)) {
    for (i2 in (i1 + 1):(length(x))) {
      comptitle <- paste0(names(x)[i1], "-", names(x)[i2])
      dd2 <- data.frame(x = dd$x, y1 = FUNC[[i1]](dd$x), 
                        y2 = FUNC[[i2]](dd$x))
      dd2[is.na(dd2)] <- 0
      dd2$ovy <- apply(dd2[, c("y1", "y2")], 1, min)
      dd2$ally <- apply(dd2[, c("y1", "y2")], 1, max, na.rm = TRUE)
      dd2$dominance <- ifelse(dd2$y1 > dd2$y2, 1, 2)
      dd2$k <- comptitle
      OV <- c(OV, sum(dd2$ovy, na.rm = TRUE)/sum(dd2$ally, 
                                                 na.rm = TRUE))
      dd2 <- dd2[order(dd2$x), ]
      CHANGE <- dd2$x[which(dd2$dominance[2:nrow(dd2)] != 
                              dd2$dominance[1:(nrow(dd2) - 1)])]
      xpoints <- c(xpoints, list(CHANGE))
      if (partial.plot) {
        gg <- ggplot(dd2, aes(x, dd2$y1)) + theme_bw() + 
          geom_vline(xintercept = CHANGE, lty = 2, color = "#cccccc") + 
          geom_line() + geom_line(aes(x, dd2$y2)) + geom_line(aes(x,dd2$ovy), 
                                                              color = "red") + 
          geom_line(aes(x,dd2$ally), color = "blue") + 
          ggtitle(comptitle) + 
          xlab("") + 
          ylab("") + 
          theme(plot.title = element_text(hjust = 0.5))
        print(gg)
      }
      DD <- rbind(DD, dd2)
      COMPTITLE <- c(COMPTITLE, comptitle)
    }
  }
  names(xpoints) <- names(OV) <- COMPTITLE
  if (plot) 
    print(my.final.plot(x, OV))
  return(list(DD = DD, OV = OV, xpoints = xpoints))
}


#### PERCENTAGE OVERLAP DENSITY PLOTS ####
  # sla mean ----
lat.list <- list(trpSLA, sbtrp_joined_SLA, tmp_sla, plN)
lat.list.n <- list()
lat.density.plot <- list()
for (i in 1:length(lat.list)){
  lat.list[[i]] <- lat.list[[i]]%>%
    filter(sla!=0,specific.leaf.area!=0) %>%
    rename("cardamom" = sla, "butler" = specific.leaf.area) 
  list <- list(cardamom=lat.list[[i]][["cardamom"]],
               butler=lat.list[[i]][["butler"]])
  name <- paste(i,"n",sep = "_")
  lat.list.n[[name]] <- list
  for (j in 1:length(lat.list.n)){
    plot_name <- paste("./figures/lat",i,"density_overlap.png",sep = "_")
    png(plot_name, width = 40, height = 20, units = "cm", res = 400)
    lat.density.plot[[j]] <- my.overlap(lat.list.n[[j]], plot = TRUE)
    dev.off()
  }
}

  # sla stdev ----
lat.std.list <- list(trpSTD,sbtrp_joined_STD,tmp_slastd,plN_std)
lat.std.list.n <- list()
lat.density.std.plot <- list()
for (i in 1:length(lat.std.list)){
  lat.std.list[[i]] <- lat.std.list[[i]]%>%
    filter(Standard_Deviation!=0,specific.leaf.area!=0) %>%
    rename("cardamom_std" = Standard_Deviation, 
           "butler_std" = specific.leaf.area) 
  list <- list(cardamom_std=lat.std.list[[i]][["cardamom_std"]],
               butler_std=lat.std.list[[i]][["butler_std"]])
  name <- paste(i,"n",sep = "_")
  lat.std.list.n[[name]] <- list
  for (j in 1:length(lat.std.list.n)){
    plot_name <- paste("./figures/lat_std",i,"density_overlap.png",sep = "_")
    png(plot_name, width = 40, height = 20, units = "cm", res = 400)
    lat.density.std.plot[[j]] <- my.overlap(lat.std.list.n[[j]], plot = TRUE)
    dev.off()
  }
}



### STIPPLING FOR BUTLER SLA MEAN (25th-75th and 25th-95th percentiles) ----
# creating loops for latitude ranges
lat.stp.25 <- list()
lat.stp.75 <-list()
lat.stp.95 <- list()
stp.locs.75 <- list()
stp.locs.95<- list()
diff.stp.locs <- list()
  # sla mean
lat.list <- list(tropics=trp_b,subtropics_N=sbtrpN_b,subtropics_S=sbtrpS_b,
                 temperate_N=tmp_b_N,temperate_S=tmp_b_S,pole_N=plN_b)

for (i in 1:length(lat.list)){
  crop.25 <- crop(cardamom_25th,extent[[i]])
  crop.75 <- crop(cardamom_75th, extent[[i]])
  crop.95 <- crop(cardamom_95th, extent[[i]])
  stp.25 <- raster::mask(crop.25, lat.list[[i]])
  stp.75 <- raster::mask(crop.75, lat.list[[i]])
  stp.95 <- raster::mask(crop.95, lat.list[[i]])
  stp.name25 <- paste("stp25", names(lat.list)[i],sep = ".")
  stp.name75 <- paste("stp75", names(lat.list)[i],sep = ".")
  stp.name95 <- paste("stp95",names(lat.list)[i],sep=".")
  lat.stp.25[[stp.name25]] <- stp.25
  lat.stp.75[[stp.name75]] <- stp.75
  lat.stp.95[[stp.name95]] <- stp.95
  stp.locs.75.calc <- (lat.list[[i]]>=lat.stp.25[[i]])*
    (lat.list[[i]]<= lat.stp.75[[i]])
  stp.locs.95.calc <- (lat.list[[i]]>=lat.stp.25[[i]])*
    (lat.list[[i]]<= lat.stp.95[[i]])
  stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
  stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
  stp.locs.75.p <- stp.locs.75.p[stp.locs.75.p[, "layer"] == 1,]
  stp.locs.95.p <- stp.locs.95.p[stp.locs.95.p[, "layer"] == 1,]
  stp.locs.75n <- paste("stp.locs",names(lat.list)[i],"75",
                        sep = ".")
  stp.locs.95n <- paste("stp.locs",names(lat.list)[i],"95",
                        sep = ".")
  stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
  stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
  diff.stp.locs.calc <- as.data.frame(stp.locs.95[[i]]) %>% 
    setdiff(as.data.frame(stp.locs.75[[i]])) %>% 
    as.matrix 
  diff.stp.locs.n <- paste("diff",names(lat.list)[i],sep = ".")
  diff.stp.locs[[diff.stp.locs.n]] <- diff.stp.locs.calc
  png(paste("./figures/lat.stippling",names(lat.list)[i],".png",
            sep = "."),
      width = 50,height = 25,units = "cm",res = 500)
  plot(lat.list[[i]],col=rev(brewer.pal(10,"RdYlBu")),
       legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                          side=4, font=1, line=2.3),
       main="Butler Sla Mean (stippling)\n")
  points(stp.locs.75[[i]],pch=18,cex=0.8)
  points(diff.stp.locs[[i]],pch=23,cex=0.9,col="darkgreen",bg="green")
  dev.off()
}

  # sla stdev 
lat.std.list <- list(tropics.std=trp_b_std,subtropics.std_N=sbtrpSTD_b_N,
                     subtropics.std_S=sbtrpSTD_b_S,temperate.std_N=tmpSTD_b_N,
                     temperate.std_S=tmpSTD_b_S,pole.std_N=plN_slastd_b) 
for (i in 1:length(lat.std.list)){
  crop.25 <- crop(cardamom_25th,extent[[i]])
  crop.75 <- crop(cardamom_75th, extent[[i]])
  crop.95 <- crop(cardamom_95th, extent[[i]])
  stp.25 <- raster::mask(crop.25, lat.std.list[[i]])
  stp.75 <- raster::mask(crop.75, lat.std.list[[i]])
  stp.95 <- raster::mask(crop.95, lat.std.list[[i]])
  stp.name25 <- paste("stp25", names(lat.std.list)[i],sep = ".")
  stp.name75 <- paste("stp75", names(lat.std.list)[i],sep = ".")
  stp.name95 <- paste("stp95",names(lat.std.list)[i],sep=".")
  lat.stp.25[[stp.name25]] <- stp.25
  lat.stp.75[[stp.name75]] <- stp.75
  lat.stp.95[[stp.name95]] <- stp.95
  stp.locs.75.calc <- (lat.std.list[[i]]>=lat.stp.25[[i]])*
    (lat.std.list[[i]]<= lat.stp.75[[i]])
  stp.locs.95.calc <- (lat.std.list[[i]]>=lat.stp.25[[i]])*
    (lat.std.list[[i]]<= lat.stp.95[[i]])
  stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
  stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
  stp.locs.75.p <- stp.locs.75.p[stp.locs.75.p[, "layer"] == 1,]
  stp.locs.95.p <- stp.locs.95.p[stp.locs.95.p[, "layer"] == 1,]
  stp.locs.75n <- paste("stp.locs",names(lat.std.list)[i],"75",
                        sep = ".")
  stp.locs.95n <- paste("stp.locs",names(lat.std.list)[i],"95",
                        sep = ".")
  stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
  stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
  diff.stp.locs.calc <- as.data.frame(stp.locs.95[[i]]) %>% 
    setdiff(as.data.frame(stp.locs.75[[i]])) %>% 
    as.matrix 
  diff.stp.locs.n <- paste("diff",names(lat.std.list)[i],sep = ".")
  diff.stp.locs[[diff.stp.locs.n]] <- diff.stp.locs.calc
  png(paste("./figures/lat.std.stippling",names(lat.std.list)[i],".png",
            sep = "."),
      width = 50,height = 25,units = "cm",res = 500)
  plot(lat.std.list[[i]],col=rev(brewer.pal(10,"RdYlBu")),
       legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                          side=4, font=1, line=2.3),
       main="Butler Sla Mean (stippling)\n")
  points(stp.locs.75[[i]],pch=18,cex=0.8)
  points(diff.stp.locs[[i]],pch=23,cex=0.9,col="darkgreen",bg="green")
  dev.off()
}




#### latitude STATS ####
lat.df <- list(tropics_sla=trpSLA,subtropics_sla=sbtrp_joined_SLA,
               temperate_sla=tmp_sla, pole_sla=plN)
lat.stats <- list()
for (i in 1:length(lat.df)){
  stats <- lat.df[[i]] %>%
    rename("cardamom"=sla, "butler"=specific.leaf.area)%>%
    filter(cardamom!=0, butler!=0) %>%
    mutate(mean_c = mean(cardamom),
           mean_b = mean(butler),
           diff_butler = butler-mean_b,
           diff_butler2 = diff_butler^2,
           sum_diff_butler2 = sum(diff_butler2),
           slope_bf = sum((cardamom-mean_c)*(butler-mean_b))/
             sum((cardamom-mean_c)^2),
           b_intercept = mean_b - (slope_bf*mean_c),
           new_b_val = b_intercept + (slope_bf*cardamom),
           dist_mean_new_b = new_b_val - mean_b,
           sqrd_dist_b = dist_mean_new_b^2,
           sum_sqrd_dist_b = sum(sqrd_dist_b),
           trps_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2,
           rmse_av = rmse(butler, cardamom),
           rmse_row = sqrt(se(butler, cardamom)),
           bias = bias(butler, cardamom),
           bias_row = butler-cardamom)
  name.lat <- paste("stat",names(lat.df)[i],sep = ".")
  lat.stats[[name.lat]] <- stats
}

# HEATSCATTER: saving figure for sla mean by latitudinal range----
png("./figures/heatsc_latitude_sla.png", width = 40, height = 30, 
    units = "cm", res = 500)
par(mfcol = c(2,2))
# (heatscatter) saving figure for sla stdev by latitudinal range----
png("./figures/heatsc_latitude_std.png", width = 40, height = 30,
    units = "cm", res = 500)
par(mfcol = c(2,2))


#### ANALYSIS TROPICS ####
# heatscatter tropics ----
# sla mean
trp_c_sla_n <- trpSLA$cardamom
trp_b_sla_n <- trpSLA$butler
(heatsc_sla_trp <- heatscatter(trp_c_sla_n, trp_b_sla_n, pch = 19, 
                                cexplot = 0.5, colpal="spectral", 
                                #disco() for all color options / could set to colorblind
                                add.contour=TRUE, main = "Tropics SLA Mean",
                                xlab="", 
                                ylab="Butler"))

# sla stdev 
trp_c_std_n <- trpSTD$cardamom_std
trp_b_std_n <- trpSTD$butler_std

(heatsc_slastd_trp <- heatscatter(trp_c_std_n, trp_b_std_n, pch = 19, 
                               cexplot = 0.5, colpal="spectral", 
                               #disco() for all color options / could set to colorblind
                               add.contour=TRUE, main = "Tropics SLA StDev",
                               xlab="", 
                               ylab="Butler"))

#### ANALYSIS SUBTROPICS ####
# 25th percentile 
sbtrpN_25pc <- crop(cardamom_25th, sbtrpextN)
sbtrpS_25pc <- crop(cardamom_25th, sbtrpextS)
sbtrpN_25pc <- raster::mask(sbtrpN_25pc, sbtrpN)
# 75th percentile
sbtrpN_75pc <- raster::mask(sbtrpN_25pc, sbtrpN)
# 95th percentile 
sbtrpN_95pc <- raster::mask(sbtrpN_25pc, sbtrpN)
sbtrp_stp75pc <- (sbtrpN[[1]] >= sbtrpN_25pc[[1]])*
  (sbtrpN[[1]]<=sbtrpN_75pc[[1]])
sbtrp_stp75pc <- rasterToPoints(sbtrp_stp75pc)
sbtrp_stp75pc <- sbtrp_stp75pc[sbtrp_stp75pc[, "layer"] == 1,]
trp_stp_95 <- (trp_b[[1]]>=trp_25pc[[1]])*
  (trp_b[[1]]<=trp_95pc[[1]])
trp_stp_95 <- rasterToPoints(trp_stp_95)
trp_stp_95 <- trp_stp_95[trp_stp_95[, "layer"] ==1,]
# heatscatter subtrps ----
# sla mean
sbtrp_c_sla_n <- sbtrp_joined_SLA$cardamom
sbtrp_b_sla_n <- sbtrp_joined_SLA$butler

(heatsc_sla_sbtrp <- heatscatter(sbtrp_c_sla_n, sbtrp_b_sla_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Subtropics",
                                  xlab="Cardamom", 
                                  ylab="Butler"))
#sla stdev 
sbtrp_c_slastd_n <- sbtrp_joined_STD$cardamom_std
sbtrp_b_slastd_n <- sbtrp_joined_STD$butler_std

(heatsc_slastd_sbtrp <- heatscatter(sbtrp_c_slastd_n, sbtrp_b_slastd_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Subtropics",
                                  xlab="Cardamom", 
                                  ylab="Butler"))
#### ANALYSIS TEMPERATE ####
# heatscatter ----
  #sla mean
tmp_sla_c_n <- tmp_sla$cardamom
tmp_sla_b_n <- tmp_sla$butler
(heatsc_sla_tmp <- heatscatter(tmp_sla_c_n, tmp_sla_b_n, pch = 19, 
                               cexplot = 0.5, colpal="spectral", 
                               #disco() for all color options / could set to colorblind
                               add.contour=TRUE, main = "Temperate",
                               xlab="", 
                               ylab=""))
  # sla stdev
tmp_slastd_c_n <- tmp_slastd$cardamom_std
tmp_slastd_b_n <- tmp_slastd$butler_std
(heatsc_slastd_tmp <- heatscatter(tmp_slastd_c_n, tmp_slastd_b_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Temperate",
                                  xlab="", 
                                  ylab=""))
#### ANALYSIS POLES ####
# heatscatter ----
  # sla mean 
pl_sla_c_n <- plN$cardamom
pl_sla_b_n <- plN$butler
(heatsc_sla_pl <- heatscatter(pl_sla_c_n, pl_sla_b_n, pch = 19, 
                               cexplot = 0.5, colpal="spectral", 
                               add.contour=TRUE, main = "N pole",
                               xlab="Cardamom", 
                               ylab=""))
# sla stdev
pl_slastd_c_n <- plN_std$cardamom_std
pl_slastd_b_n <- plN_std$butler_std
(heatsc_slastd_pl <- heatscatter(pl_slastd_c_n, pl_slastd_b_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "N Pole",
                                  xlab="Cardamom", 
                                  ylab=""))
# dev.off----
dev.off()


#############################
#### SPLITTING BY BIOME #####
#############################
# OPEN DATASET: ECOREGIONS17 - where i get the biomes from ----
ecoregions17 <- st_read("./DATA/Ecoregions2017/Ecoregions2017.shp")
st_crs(ecoregions17)
st_bbox(ecoregions17)
#      xmin       ymin       xmax       ymax 
#-179.99999  -89.89197  180.00000   83.62313 

(plot_ecoregion17 <- ggplot() + 
                  geom_sf(data = ) + 
                  ggtitle("Ecoregions 2017") + 
                  coord_sf())

unique(ecoregions17[[4]])
ecoregions17_geom <- st_geometry(ecoregions17)
ecoregions17_geom[[1]]

# creation of different raster layers by major biome types ----
boreal_f_taiga <- ecoregions17 %>%
  filter(BIOME_NAME == "Boreal Forests/Taiga")
tundra <- ecoregions17 %>%
  filter(BIOME_NAME == "Tundra")
temp_conif_forest <- ecoregions17 %>%
  filter(BIOME_NAME == "Temperate Conifer Forests")
temp_broad_mix <- ecoregions17 %>%
  filter(BIOME_NAME == "Temperate Broadleaf & Mixed Forests")
trp_sbtrp_dry_broad <- ecoregions17 %>%
  filter(BIOME_NAME == "Tropical & Subtropical Dry Broadleaf Forests")
trp_sbtrp_conif <- ecoregions17 %>%
  filter(BIOME_NAME == "Tropical & Subtropical Coniferous Forests")
trp_sbtrp_moist_broad <- ecoregions17 %>%
  filter(BIOME_NAME == "Tropical & Subtropical Moist Broadleaf Forests")
med_f_w_scr <- ecoregions17 %>%
  filter(BIOME_NAME == "Mediterranean Forests, Woodlands & Scrub")
des_x_scr <- ecoregions17 %>%
  filter(BIOME_NAME == "Deserts & Xeric Shrublands")
temp_grass_sav_shr <- ecoregions17 %>%
  filter(BIOME_NAME == "Temperate Grasslands, Savannas & Shrublands")
mont_grass_shr <- ecoregions17 %>%
  filter(BIOME_NAME == "Montane Grasslands & Shrublands")
mangroves <- ecoregions17 %>%
  filter(BIOME_NAME == "Mangroves")
flo_grass_sav <- ecoregions17 %>%
  filter(BIOME_NAME == "Flooded Grasslands & Savannas")
trp_sbtrp_grass_sav_shr <- ecoregions17 %>%
  filter(BIOME_NAME == "Tropical & Subtropical Grasslands, Savannas & Shrublands")

# 1) mask raster taiga biome  ----
# sla mean
masked_taiga_sla_c <- raster::mask(cardamom_sla, boreal_f_taiga)
masked_taiga_sla_b <- raster::mask(butler_sla, boreal_f_taiga)
#plot(masked_taiga_sla_c[[1]])
#plot(masked_taiga_sla_b[[1]])

# sla stdev 
masked_taiga_slastd_c <- raster::mask(cardamom_sla_std, boreal_f_taiga)
masked_taiga_slastd_b <- raster::mask(butler_sla_std, boreal_f_taiga)
#plot(masked_taiga_slastd_c[[1]])
#plot(masked_taiga_slastd_b[[1]])
# 2) mask raster tundra biome ----
# sla mean 
masked_tundra_sla_c <- raster::mask(cardamom_sla, tundra)
masked_tundra_sla_b <- raster::mask(butler_sla, tundra)
#plot(masked_tundra_sla_c[[1]])
#plot(masked_tundra_sla_b[[1]])

# sla stdev 
masked_tundra_slastd_c <- raster::mask(cardamom_sla_std, tundra)
masked_tundra_slastd_b <- raster::mask(butler_sla_std, tundra)
#plot(masked_tundra_slastd_c[[1]])
#plot(masked_tundra_slastd_b[[1]])

# 3) mask raster temp conif forest biome ----
# sla mean 
mask_temp_conif_sla_c <- raster::mask(cardamom_sla,temp_conif_forest)
mask_temp_conif_sla_b <- raster::mask(butler_sla, temp_conif_forest)
#plot(mask_temp_conif_sla_c[[1]])
#plot(mask_temp_conif_sla_b[[1]])

# sla stdev 
mask_temp_conif_slastd_c <- raster::mask(cardamom_sla_std, temp_conif_forest)
mask_temp_conif_slastd_b <- raster::mask(butler_sla_std, temp_conif_forest)
#plot(mask_temp_conif_slastd_c[[1]])
#plot(mask_temp_conif_slastd_b[[1]])

# 4) mask raster temperate broad and mixed forest biome ----
# sla mean
mask_temp_broad_mix_sla_c <- raster::mask(cardamom_sla, temp_broad_mix)
mask_temp_broad_mix_sla_b <- raster::mask(butler_sla, temp_broad_mix)
#plot(mask_temp_broad_mix_sla_c[[1]])
#plot(mask_temp_broad_mix_sla_b[[1]])

# sla stdev 
mask_temp_broad_mix_slastd_c <- raster::mask(cardamom_sla_std, temp_broad_mix)
mask_temp_broad_mix_slastd_b <- raster::mask(butler_sla_std, temp_broad_mix)
#plot(mask_temp_broad_mix_slastd_c[[1]])
#plot(mask_temp_broad_mix_slastd_b[[1]])

# 5) mask raster tropical and subtropical dry broadleaf biome ----
# sla mean
mask_trp_sbtrp_dry_broad_sla_c <- raster::mask(cardamom_sla, 
                                               trp_sbtrp_dry_broad)
mask_trp_sbtrp_dry_broad_sla_b <-raster::mask(butler_sla, trp_sbtrp_dry_broad)
#plot(mask_trp_sbtrp_dry_broad_sla_c[[1]])
#plot(mask_trp_sbtrp_dry_broad_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_dry_broad_slastd_c <- raster::mask(cardamom_sla_std, 
                                                  trp_sbtrp_dry_broad)
mask_trp_sbtrp_dry_broad_slastd_b <- raster::mask(butler_sla_std, 
                                                  trp_sbtrp_dry_broad)
#plot(mask_trp_sbtrp_dry_broad_slastd_c[[1]])
#plot(mask_trp_sbtrp_dry_broad_slastd_b[[1]])

# 6) mask raster tropical and subtropical conif forest biome ----
# sla mean 
mask_trp_sbtrp_conif_sla_c <- raster::mask(cardamom_sla, trp_sbtrp_conif)
mask_trp_sbtrp_conif_sla_b <- raster::mask(butler_sla, trp_sbtrp_conif)
#plot(mask_trp_sbtrp_conif_sla_c[[1]])
#plot(mask_trp_sbtrp_conif_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_conif_slastd_c <- raster::mask(cardamom_sla_std, trp_sbtrp_conif)
mask_trp_sbtrp_conif_slastd_b <- raster::mask(butler_sla_std, trp_sbtrp_conif)
#plot(mask_trp_sbtrp_conif_slastd_c[[1]])
#plot(mask_trp_sbtrp_conif_slastd_b[[1]])

# 7) mask raster tropical subtropical moist broadleaf biome ----
# sla mean
mask_trp_sbtrp_moist_broad_sla_c <- raster::mask(cardamom_sla, 
                                                 trp_sbtrp_moist_broad)
mask_trp_sbtrp_moist_broad_sla_b <- raster::mask(butler_sla, 
                                                 trp_sbtrp_moist_broad)
#plot(mask_trp_sbtrp_moist_broad_sla_c[[1]])
#plot(mask_trp_sbtrp_moist_broad_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_moist_broad_slastd_c <- raster::mask(cardamom_sla_std, 
                                                    trp_sbtrp_moist_broad)
mask_trp_sbtrp_moist_broad_slastd_b <- raster::mask(butler_sla_std, 
                                                    trp_sbtrp_moist_broad)
#plot(mask_trp_sbtrp_moist_broad_slastd_c[[1]])
#plot(mask_trp_sbtrp_moist_broad_slastd_b[[1]])

### these might be extras ###
# 8) mask raster mediterranean forests, woodlands, scrub biome ----
# sla mean
mask_med_f_w_scr_sla_c <- raster::mask(cardamom_sla, med_f_w_scr)
mask_med_f_w_scr_sla_b <- raster::mask(butler_sla, med_f_w_scr)
#plot(mask_med_f_w_scr_sla_c[[1]])
#plot(mask_med_f_w_scr_sla_b[[1]])

# sla stdev 
mask_med_f_w_scr_slastd_c <- raster::mask(cardamom_sla_std, med_f_w_scr)
mask_med_f_w_scr_slastd_b <- raster::mask(butler_sla_std, med_f_w_scr)
#plot(mask_med_f_w_scr_slastd_c[[1]])
#plot(mask_med_f_w_scr_slastd_b[[1]])

# 9) mask raster desertic and xeric scrubland biome ----
# sla mean
mask_des_x_scr_sla_c <- raster::mask(cardamom_sla, des_x_scr)
mask_des_x_scr_sla_b <- raster:: mask(butler_sla, des_x_scr)
#plot(mask_des_x_scr_sla_c[[1]])
#plot(mask_des_x_scr_sla_b[[1]])

# sla stdev 
mask_des_x_scr_slastd_c <- raster::mask(cardamom_sla_std, des_x_scr)
mask_des_x_scr_slastd_b <- raster::mask(butler_sla_std, des_x_scr)
#plot(mask_des_x_scr_slastd_c[[1]])
#plot(mask_des_x_scr_slastd_b[[1]])

# 10) mask raster temperate grassland, savanna, shrubland biome ----
# sla mean
mask_temp_grass_sav_shr_sla_c <- raster::mask(cardamom_sla,temp_grass_sav_shr)
mask_temp_grass_sav_shr_sla_b <- raster::mask(butler_sla, temp_grass_sav_shr)
#plot(mask_temp_grass_sav_shr_sla_c[[1]])
#plot(mask_temp_grass_sav_shr_sla_b[[1]])

# sla stdev
mask_temp_grass_sav_shr_slastd_c <- raster::mask(cardamom_sla_std, 
                                                 temp_grass_sav_shr)
mask_temp_grass_sav_shr_slastd_b <- raster::mask(butler_sla_std, 
                                                 temp_grass_sav_shr)
#plot(mask_temp_grass_sav_shr_slastd_c[[1]])
#plot(mask_temp_grass_sav_shr_slastd_b[[1]])

# 11) mask raster montane grassland and shrubland biome ----
# sla mean 
mask_mont_grass_shr_sla_c <- raster::mask(cardamom_sla, mont_grass_shr)
mask_mont_grass_shr_sla_b <- raster::mask(butler_sla, mont_grass_shr)
#plot(mask_mont_grass_shr_sla_c[[1]])
#plot(mask_mont_grass_shr_sla_b[[1]])

# sla stdev 
mask_mont_grass_shr_slastd_c <- raster::mask(cardamom_sla_std, 
                                             mont_grass_shr)
mask_mont_grass_shr_slastd_b <- raster::mask(butler_sla_std, 
                                             mont_grass_shr)
#plot(mask_mont_grass_shr_slastd_c[[1]])
#plot(mask_mont_grass_shr_slastd_b[[1]])

# 12) mask raster mangrove biome ----
# sla mean
mask_mangroves_sla_c <- raster::mask(cardamom_sla, mangroves)
mask_mangroves_sla_b <- raster::mask(butler_sla, mangroves)
#plot(mask_mangroves_sla_c[[1]])
#plot(mask_mangroves_sla_b[[1]])

# sla stdev 
mask_mangroves_slastd_c <- raster::mask(cardamom_sla_std, mangroves)
mask_mangroves_slastd_b <- raster::mask(butler_sla_std, mangroves)
#plot(mask_mangroves_slastd_c[[1]])
#plot(mask_mangroves_slastd_b[[1]])

# 13) mask raster flooded grassland and savanna biome ----
# sla mean
mask_flo_grass_sav_sla_c <- raster::mask(cardamom_sla, flo_grass_sav)
mask_flo_grass_sav_sla_b <- raster::mask(butler_sla, flo_grass_sav)
#plot(mask_flo_grass_sav_sla_c[[1]])
#plot(mask_flo_grass_sav_sla_b[[1]])

# sla stdev 
mask_flo_grass_sav_slastd_c <- raster::mask(cardamom_sla_std, flo_grass_sav)
mask_flo_grass_sav_slastd_b <- raster::mask(butler_sla_std, flo_grass_sav)
#plot(mask_flo_grass_sav_slastd_c[[1]])
#plot(mask_flo_grass_sav_slastd_b[[1]])

# 14) mask raster tropical and subtropical grassland, savanna, shrubland biome ----
# sla mean 
mask_trp_sbtrp_grass_sav_shr_sla_c <- raster::mask(cardamom_sla, 
                                                   trp_sbtrp_grass_sav_shr)
mask_trp_sbtrp_grass_sav_shr_sla_b <- raster::mask(butler_sla, 
                                                   trp_sbtrp_grass_sav_shr)
#plot(mask_trp_sbtrp_grass_sav_shr_sla_c[[1]])
#plot(mask_trp_sbtrp_grass_sav_shr_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_grass_sav_shr_slastd_c <- raster::mask(cardamom_sla_std,
                                                      trp_sbtrp_grass_sav_shr)
mask_trp_sbtrp_grass_sav_shr_slastd_b <- raster::mask(butler_sla_std,
                                                      trp_sbtrp_grass_sav_shr)
#plot(mask_trp_sbtrp_grass_sav_shr_slastd_c[[1]])
#plot(mask_trp_sbtrp_grass_sav_shr_slastd_b[[1]])






#### VISUAL AND STAT ANALYSIS BY BIOME ####
# first thing - turning all masked rasterlayers to dataframes ----
mask.to.df <- function(x){
  new.list <- list()
  for (i in 1:length(x)){
    df <- raster::as.data.frame(x[[i]], xy =TRUE)
    name.df <- paste("df",names(x)[i],sep = ".")
    new.list[[name.df]] <- df
  }
  new.list
}
  # sla mean - butler
biome.butler.raster.list <- list(taiga.b=masked_taiga_sla_b,
                          tundra.b= masked_tundra_sla_b,
                          des.b = mask_des_x_scr_sla_b,
                          trp.g.b = mask_trp_sbtrp_grass_sav_shr_sla_b,
                          temp.c.b=mask_temp_conif_sla_b,
                          med.b=mask_med_f_w_scr_sla_b,
                          mang.b=mask_mangroves_sla_b,
                          temp.bm.b=mask_temp_broad_mix_sla_b,
                          trp.c.b=mask_trp_sbtrp_conif_sla_b,
                          flo.g.b=mask_flo_grass_sav_sla_b,
                          mont.g.b=mask_mont_grass_shr_sla_b,
                          trp.db.b=mask_trp_sbtrp_dry_broad_sla_b,
                          trp.mb.b=mask_trp_sbtrp_moist_broad_sla_b,
                          temp.g.b=mask_temp_grass_sav_shr_sla_b)

biome.butler.df <- mask.to.df(biome.butler.raster.list)

  # sla mean - cardamom
biome.cardamom.raster.list <- list(taiga.c=masked_taiga_sla_c,
                                 tundra.c= masked_tundra_sla_c,
                                 des.c = mask_des_x_scr_sla_c,
                                 trp.g.c = mask_trp_sbtrp_grass_sav_shr_sla_c,
                                 temp.c.c=mask_temp_conif_sla_c,
                                 med.c=mask_med_f_w_scr_sla_c,
                                 mang.c=mask_mangroves_sla_c,
                                 temp.bm.c=mask_temp_broad_mix_sla_c,
                                 trp.c.c=mask_trp_sbtrp_conif_sla_c,
                                 flo.g.c=mask_flo_grass_sav_sla_c,
                                 mont.g.c=mask_mont_grass_shr_sla_c,
                                 trp.db.c=mask_trp_sbtrp_dry_broad_sla_c,
                                 trp.mb.c=mask_trp_sbtrp_moist_broad_sla_c,
                                 temp.g.c=mask_temp_grass_sav_shr_sla_c)

biome.cardamom.df <- mask.to.df(biome.cardamom.raster.list)

  # sla stdev - butler 
biome.cardamom.raster.std <- list(taiga.c=masked_taiga_slastd_c,
                                   tundra.c= masked_tundra_slastd_c,
                                   des.c = mask_des_x_scr_slastd_c,
                                   trp.g.c = mask_trp_sbtrp_grass_sav_shr_slastd_c,
                                   temp.c.c=mask_temp_conif_slastd_c,
                                   med.c=mask_med_f_w_scr_slastd_c,
                                   mang.c=mask_mangroves_slastd_c,
                                   temp.bm.c=mask_temp_broad_mix_slastd_c,
                                   trp.c.c=mask_trp_sbtrp_conif_slastd_c,
                                   flo.g.c=mask_flo_grass_sav_slastd_c,
                                   mont.g.c=mask_mont_grass_shr_slastd_c,
                                   trp.db.c=mask_trp_sbtrp_dry_broad_slastd_c,
                                   trp.mb.c=mask_trp_sbtrp_moist_broad_slastd_c,
                                   temp.g.c=mask_temp_grass_sav_shr_slastd_c)

biome.cardamom.df.std <- mask.to.df(biome.cardamom.raster.std)

  # sla stdev - cardamom 
biome.butler.raster.std <- list(taiga.b=masked_taiga_slastd_b,
                                tundra.b= masked_tundra_slastd_b,
                                des.b = mask_des_x_scr_slastd_b,
                                trp.g.b = mask_trp_sbtrp_grass_sav_shr_slastd_b,
                                temp.c.b=mask_temp_conif_slastd_b,
                                med.b=mask_med_f_w_scr_slastd_b,
                                mang.b=mask_mangroves_slastd_b,
                                temp.bm.b=mask_temp_broad_mix_slastd_b,
                                trp.c.b=mask_trp_sbtrp_conif_slastd_b,
                                flo.g.b=mask_flo_grass_sav_slastd_b,
                                mont.g.b=mask_mont_grass_shr_slastd_b,
                                trp.db.b=mask_trp_sbtrp_dry_broad_slastd_b,
                                trp.mb.b=mask_trp_sbtrp_moist_broad_slastd_b,
                                temp.g.b=mask_temp_grass_sav_shr_slastd_b)
biome.butler.df.std <- mask.to.df(biome.butler.raster.std)

# joining dataframes cardamom + butler ----
join.f <- function(x,y,na.omit=TRUE){
  new.list <- list()
  join <- mapply(left_join, x, y,SIMPLIFY = FALSE)
  join <- lapply(join, na.omit)
  name.df <- names(x)[i]
  new.list[[name.df]] <- join
}

# sla mean
j.biome.sla <- join.f(biome.cardamom.df,biome.butler.df, na.omit = TRUE)

# sla stdev 
j.biome.slastd <- join.f(biome.cardamom.df.std,biome.butler.df.std)


#### PERCENTAGE OVERLAP DENSITY PLOTS ####
  # sla mean ----
biome.list <- list(taiga_sla=j_taiga_sla, 
                   tundra_sla=j_tundra_sla, tmp_c_sla=j_tmp_c_sla,
                   tmp_b_m_sla=j_tmp_b_m_sla, 
                   trp_sbtrp_d_b_sla=j_trp_sbtrp_d_b_sla,
                   trp_sbtrp_c_sla=j_trp_sbtrp_c_sla, 
                   trp_sbtrp_m_br_sla=j_trp_sbtrp_m_br_sla,
                   med_f_sla=j_med_f_sla, des_x_s_sla=j_des_x_s_sla, 
                   temp_g_s_sh_sla=j_temp_g_s_sh_sla,
                   mont_g_shr_sla=j_mont_g_shr_sla,
                   mangr_sla=j_mangr_sla, flo_g_sav_sla=j_flo_g_sav_sla,
                   trpsbtrp_g_sav_shr_sla=j_trpsbtrp_g_sav_shr_sla)
biome.list.n <- list()
density.plot <- list()

for (i in 1:length(biome.list)){
  biome.list[[i]] <- biome.list[[i]]%>%
    filter(sla!=0,specific.leaf.area!=0) %>%
    rename("cardamom" = sla, "butler" = specific.leaf.area) 
  list <- list(cardamom=biome.list[[i]][["cardamom"]],
               butler=biome.list[[i]][["butler"]])
  name <- paste(names(biome.list)[i],"n",sep = "_")
  biome.list.n[[name]] <- list
  for (j in 1:length(biome.list.n)){
    plot_name <- paste("./figures/biome",names(biome.list)[i], # this is the way to do it
                       "density_overlap.png",sep = "_")
    png(plot_name, width = 40, height = 20, units = "cm", res = 400)
    density.plot[[j]] <- my.overlap(biome.list.n[[j]], plot = TRUE)
   # print(density.plot[[j]])
    dev.off()
  }
}

  # sla stdev ----
biome.std.list <- list(taiga_slastd=j_taiga_slastd, 
                       tundra_slastd=j_tundra_slastd,
                       tmp_c_slastd=j_tmp_c_slastd, 
                       tmp_b_m_slastd=j_tmp_b_m_slastd,
                       trp_sbtrp_d_b_slastd=j_trp_sbtrp_d_b_slastd, 
                       trp_sbtrp_c_slastd=j_trp_sbtrp_c_slastd,
                       trp_sbtrp_m_br_slastd=j_trp_sbtrp_m_br_slastd, 
                       med_f_slastd=j_med_f_slastd,
                       des_x_s_slastd=j_des_x_s_slastd, 
                       temp_g_s_sh_slastd=j_temp_g_s_sh_slastd,
                       mont_g_shr_slastd=j_mont_g_shr_slastd, 
                       mangr_slastd=j_mangr_slastd,
                       flo_g_sav_slastd=j_flo_g_sav_slastd, 
                       trpsbtrp_g_sav_shr_slastd=j_trpsbtrp_g_sav_shr_slastd)
biome.std.list.n <- list()
density.std.plot <- list()
for (i in 1:length(biome.std.list)){
  biome.std.list[[i]] <- biome.std.list[[i]]%>%
    filter(Standard_Deviation!=0,specific.leaf.area!=0) %>%
    rename("cardamom" = Standard_Deviation, "butler" = specific.leaf.area) 
  list <- list(cardamom=biome.std.list[[i]][["cardamom"]],
               butler=biome.std.list[[i]][["butler"]])
  name <- paste(names(biome.std.list)[i],"n",sep = "_")
  biome.std.list.n[[name]] <- list
  for (j in 1:length(biome.std.list.n)){
    # density.plot.name <- paste("biome",j,"density_plot", sep = "_")
    plot_name <- paste("./figures/biome_std",names(biome.std.list)[i],
                       "density_overlap.png",sep = "_")
    png(plot_name, width = 40, height = 20, units = "cm", res = 400)
    density.std.plot[[j]] <- my.overlap(biome.std.list.n[[j]], plot = TRUE)
    # print(density.plot[[j]])
    dev.off()
  }
}



#### STIPPLING BY BIOME ####
# list of masks by biome only for butler 

# list of masks by biome only for butler (std)
biome.std.raster.list <- list(m.std.taiga.b = masked_taiga_slastd_b,
                            m.std.tundra.b = masked_tundra_slastd_b,
                            m.std.des.b = mask_des_x_scr_slastd_b,
                            m.std.trp.g.b = mask_trp_sbtrp_grass_sav_shr_slastd_b,
                            m.std.temp.c.b = mask_temp_conif_slastd_b,
                            m.std.med.b = mask_med_f_w_scr_slastd_b,
                            m.std.mangr.b = mask_mangroves_slastd_b,
                            m.std.temp.bm.b = mask_temp_broad_mix_slastd_b,
                            m.std.trp.c.b = mask_trp_sbtrp_conif_slastd_b,
                            m.std.flo.g.b = mask_flo_grass_sav_slastd_b,
                            m.std.mont.g.b = mask_mont_grass_shr_slastd_b,
                            m.std.trp.db.b = mask_trp_sbtrp_dry_broad_slastd_b,
                            m.std.trp.mb.b = mask_trp_sbtrp_moist_broad_slastd_b,
                            m.std.temp.g.b = mask_temp_grass_sav_shr_slastd_b)

## creation of function to do the automatisation of stippling by biome
stippling.biomes <- function(x){
  biome.stp.25 <- list()
  biome.stp.75 <-list()
  biome.stp.95 <- list()
  stp.locs.75 <- list()
  stp.locs.95<- list()
  diff.stp.locs <- list()
  for (i in 1:length(x)){
    stp.25 <- raster::mask(cardamom_25th, x[[i]])
    stp.75 <- raster::mask(cardamom_75th, x[[i]])
    stp.95 <- raster::mask(cardamom_95th, x[[i]])
    stp.name25 <- paste("stp25", names(x)[i],sep = ".")
    stp.name75 <- paste("stp75", names(x)[i],sep = ".")
    stp.name95 <- paste("stp95",names(x)[i],sep=".")
    biome.stp.25[[stp.name25]] <- stp.25
    biome.stp.75[[stp.name75]] <- stp.75
   biome.stp.95[[stp.name95]] <- stp.95
   stp.locs.75.calc <- (x[[i]]>=biome.stp.25[[i]])*2
    (x[[i]]<= biome.stp.75[[i]])
  stp.locs.95.calc <- (x[[i]]>=biome.stp.25[[i]])*
    (x[[i]]<= biome.stp.95[[i]])
  stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
  stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
  stp.locs.75.p <- stp.locs.75.p[stp.locs.75.p[, "layer"] == 1,]
  stp.locs.95.p <- stp.locs.95.p[stp.locs.95.p[, "layer"] == 1,]
  stp.locs.75n <- paste("stp.locs",names(x)[i],"75",sep = ".")
  stp.locs.95n <- paste("stp.locs",names(x)[i],"95",sep = ".")
  stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
  stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
  diff.stp.locs.calc <- as.data.frame(stp.locs.95[[i]]) %>% 
    setdiff(as.data.frame(stp.locs.75[[i]])) %>% 
    as.matrix 
  diff.stp.locs.n <- paste("diff",names(x)[i],sep = ".")
    diff.stp.locs[[diff.stp.locs.n]] <- diff.stp.locs.calc
  png(paste("./figures/stippling",names(x)[i],".png",sep = "."),
      width = 50,height = 25,units = "cm",res = 500)
  plot(x[[i]],asp=NA,col=rev(brewer.pal(10,"RdYlBu")),
       legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                          side=4, font=1, line=2.3),
       main="Butler Sla Mean (stippling)\n")
  points(stp.locs.75[[i]],pch=18,cex=0.8)
  points(diff.stp.locs[[i]],pch=23,cex=0.9,col="darkgreen",bg="green")}
  dev.off()
}

### note ###
# there is something wrong when i try to apply this function to stdev, 
# for some reason it seems that the format of one particular biome becomes 
# different from others, thus doesnt work and loops stops there once it reaches 
# that biome - m.std.trp.db.b (the tropics-subtropics dry broadleaf biome)
# if have time, can check better and try to fix it
# for the time being, i've copied the code to match the std, without the 
# function, and without the code for finding diff between 25-75pc and 25-95pc
# since stdev doesnt have points that go beyond the 75, as i had found out from 
# the observation of global stippling 
###---###

# stippling biomes for sla mean
stippling.biomes(biome.raster.list)

# stippling biomes for sla stdev
biome.stp.25 <- list()
biome.stp.75 <-list()
biome.stp.95 <- list()
stp.locs.75 <- list()
stp.locs.95<- list()
diff.stp.locs <- list()
for (i in 1:length(biome.std.raster.list)){
    stp.25 <- raster::mask(cardamom_25th, biome.std.raster.list[[i]])
    stp.75 <- raster::mask(cardamom_75th, biome.std.raster.list[[i]])
    stp.95 <- raster::mask(cardamom_95th, biome.std.raster.list[[i]])
    stp.name25 <- paste("stp25", names(biome.std.raster.list)[i],sep = ".")
    stp.name75 <- paste("stp75", names(biome.std.raster.list)[i],sep = ".")
    stp.name95 <- paste("stp95",names(biome.std.raster.list)[i],sep=".")
    biome.stp.25[[stp.name25]] <- stp.25
    biome.stp.75[[stp.name75]] <- stp.75
    biome.stp.95[[stp.name95]] <- stp.95
    stp.locs.75.calc <- (biome.std.raster.list[[i]]>=biome.stp.25[[i]])*
      (biome.std.raster.list[[i]]<= biome.stp.75[[i]])
    stp.locs.95.calc <- (biome.std.raster.list[[i]]>=biome.stp.25[[i]])*
      (biome.std.raster.list[[i]]<= biome.stp.95[[i]])
    stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
    stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
    stp.locs.75.p <- stp.locs.75.p[stp.locs.75.p[, "layer"] == 1,]
    stp.locs.95.p <- stp.locs.95.p[stp.locs.95.p[, "layer"] == 1,]
    stp.locs.75n <- paste("stp.locs",names(biome.std.raster.list)[i],"75",
                          sep = ".")
    stp.locs.95n <- paste("stp.locs",names(biome.std.raster.list)[i],"95",
                          sep = ".")
    stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
    stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
    png(paste("./figures/stippling",names(biome.std.raster.list)[i],".png",
              sep = "."),
        width = 50,height = 25,units = "cm",res = 500)
    plot(biome.std.raster.list[[i]],asp=NA,col=rev(brewer.pal(10,"RdYlBu")),
         legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                            side=4, font=1, line=2.3),
         main="Butler Sla Mean (stippling)\n")
    points(stp.locs.75[[i]],pch=18,cex=0.8)
    dev.off()
}






# HEATSCATTER SLA MEAN MAJ BIOMES ----
png("./figures/panel_heatsc_biome_sla.png", width = 50, height = 30,
    units = "cm", res = 500)
par(mfcol=c(2,4))
# 1) taiga ----
taiga_c_n <- j_taiga_sla$sla
taiga_b_n <- j_taiga_sla$specific.leaf.area
(heatsc_taiga_sla <- heatscatter(taiga_c_n, taiga_b_n, 
                                 pch = 19, cexplot = 0.5, colpal="spectral", 
                                 add.contour=TRUE, main = "Taiga",
                                 xlab=" ", 
                                 ylab="Butler"))
  
# 2) tundra ----
tundra_c_n <- j_tundra_sla$sla
tundra_b_n <- j_tundra_sla$specific.leaf.area
(heatsc_taiga_sla <- heatscatter(tundra_c_n, tundra_b_n, 
                                 pch = 19, cexplot = 0.5, colpal="spectral", 
                                 add.contour=TRUE, main = "Tundra",
                                 xlab="Cardamom", 
                                 ylab="Butler"))

# 3) temp conif forest ----
temp_con_c_n <- j_tmp_c_sla$sla
temp_con_b_n <- j_tmp_c_sla$specific.leaf.area
(heatsc_temp_con <- heatscatter(temp_con_c_n, temp_con_b_n, 
                                 pch = 19, cexplot = 0.5, colpal="spectral", 
                                 add.contour=TRUE, main = "Temperate coniferous",
                                 xlab=" ", 
                                 ylab=" "))

# 4) temp broad mix forest ----
temp_broad_mix_c_n <- j_tmp_b_m_sla$sla
temp_broad_mix_b_n <- j_tmp_b_m_sla$specific.leaf.area
(heatsc_temp_con <- heatscatter(temp_broad_mix_c_n, temp_broad_mix_b_n, 
                                pch = 19, cexplot = 0.5, colpal="spectral", 
                                add.contour=TRUE, main = "Temperate broad/mixed",
                                xlab="Cardamom", 
                                ylab=" "))

# 5) tropical and subtropical dry broadleaf ----
trpsbtrp_dry_broad_c_n <- j_trp_sbtrp_d_b_sla$sla
trpsbtrp_dry_broad_b_n <- j_trp_sbtrp_d_b_sla$specific.leaf.area
(heatsc_trpsbtrp_d_broad <- heatscatter(trpsbtrp_dry_broad_c_n, 
                                        trpsbtrp_dry_broad_b_n, 
                                pch = 19, cexplot = 0.5, colpal="spectral", 
                                add.contour=TRUE, 
                                main = "Tropical/subtropical dry broadleaf",
                                xlab=" ", 
                                ylab=" "))

# 6) tropical and subtropical conif forest ----
trpsbtrp_con_c_n <- j_trp_sbtrp_c_sla$sla
trpsbtrp_con_b_n <- j_trp_sbtrp_c_sla$specific.leaf.area
(heatsc_trpsbtrp_con <- heatscatter(trpsbtrp_con_c_n, 
                                    trpsbtrp_con_b_n, 
                                pch = 19, cexplot = 0.5, colpal="spectral", 
                                add.contour=TRUE, main = "Tropical/subtropical coniferous",
                                xlab="Cardamom", 
                                ylab=" "))

# 7) tropical subtropical moist broadleaf ----
trpsbtrp_m_broad_c_n <- j_trp_sbtrp_m_br_sla$sla
trpsbtrp_m_broad_b_n <- j_trp_sbtrp_m_br_sla$specific.leaf.area
(heatsc_trpsbtrp_m_broad <- heatscatter(trpsbtrp_m_broad_c_n, 
                                        trpsbtrp_m_broad_b_n, 
                                    pch = 19, cexplot = 0.5, colpal="spectral", 
                                    add.contour=TRUE, 
                                    main = "Tropical/subtropical moist broadleaf",
                                    xlab=" ", 
                                    ylab=" "))

# 8) mediterranean forests, woodlands, scrub ----
med_f_w_scr_c_n <- j_med_f_sla$sla
med_f_w_scr_b_n <- j_med_f_sla$specific.leaf.area
(heatsc_med_f_w_scr <- heatscatter(med_f_w_scr_c_n, 
                                   med_f_w_scr_b_n, 
                                        pch = 19, cexplot = 0.5, colpal="spectral", 
                                        add.contour=TRUE, 
                                   main = "Mediterranean\nwoodland and scrubland",
                                        xlab="Cardamom", 
                                        ylab=" "))

# dev.off ----
dev.off()



# HEATSCATTER SLA STDEV MAJ BIOMES ----
png("./figures/panel_heatsc_slastd_biomes.png", width = 50,
    height = 30, units = "cm", res = 500)
par(mfcol = c(2,4))
# 1) Taiga ----
taiga_slastd_c_n <- j_taiga_slastd$Standard_Deviation
taiga_slastd_b_n <- j_taiga_slastd$specific.leaf.area
(heatsc_taiga_slastd <- heatscatter(taiga_slastd_c_n, taiga_slastd_b_n, 
                                    pch = 19, cexplot = 0.5, colpal="spectral", 
                                    add.contour=TRUE, main = "Taiga",
                                    xlab="", 
                                    ylab="Butler"))
# 2) Tundra ----
tundra_c_slastd_n <- j_tundra_slastd$Standard_Deviation
tundra_b_slastd_n <- j_tundra_slastd$specific.leaf.area
(heatsc_tundra_slastd <- heatscatter(tundra_c_slastd_n, tundra_b_slastd_n, 
                                    pch = 19, cexplot = 0.5, colpal="spectral", 
                                    add.contour=TRUE, main = "Tundra",
                                    xlab="Cardamom", 
                                    ylab="Butler"))

# 3) temp conif forest ----
temp_con_std_c_n <- j_tmp_c_slastd$Standard_Deviation
temp_con_std_b_n <- j_tmp_c_slastd$specific.leaf.area
(heatsc_temp_con_slastd <- heatscatter(temp_con_std_c_n, 
                                       temp_con_std_b_n, 
                                     pch = 19, cexplot = 0.5, colpal="spectral", 
                                     add.contour=TRUE, main = "Temperate coniferous",
                                     xlab="", 
                                     ylab=""))

# 4) temp broad mix forest ----
temp_broad_mix_std_c_n <- j_tmp_b_m_slastd$Standard_Deviation
temp_broad_mix_std_b_n <- j_tmp_b_m_slastd$specific.leaf.area
(heatsc_temp_b_m_slastd <- heatscatter(temp_broad_mix_std_c_n, 
                                       temp_broad_mix_std_b_n, 
                                       pch = 19, cexplot = 0.5, 
                                       colpal="spectral", 
                                       add.contour=TRUE, 
                                       main = "Temperate broadleaf mixed",
                                       xlab="Cardamom", 
                                       ylab=""))

# 5) tropical and subtropical dry broadleaf ----
trpsbtrp_d_broad_std_c_n <- j_trp_sbtrp_d_b_slastd$Standard_Deviation
trpsbtrp_d_broad_std_b_n <- j_trp_sbtrp_d_b_slastd$specific.leaf.area
(heatsc_trpsbtrp_d_b_slastd <- heatscatter(trpsbtrp_d_broad_std_c_n, 
                                           trpsbtrp_d_broad_std_b_n, 
                                       pch = 19, cexplot = 0.5, 
                                       colpal="spectral", 
                                       add.contour=TRUE, 
                                       main = "Tropical/subtropical dry broadleaf",
                                       xlab="", 
                                       ylab=""))

# 6) tropical and subtropical conif forest ----
trpsbtrp_con_std_c_n <- j_trp_sbtrp_c_slastd$Standard_Deviation
trpsbtrp_con_std_b_n <- j_trp_sbtrp_c_slastd$specific.leaf.area
(heatsc_trpsbtrp_con_slastd <- heatscatter(trpsbtrp_con_std_c_n, 
                                           trpsbtrp_con_std_b_n, 
                                           pch = 19, cexplot = 0.5, 
                                           colpal="spectral", 
                                           add.contour=TRUE, 
                                           main = "Tropical/subtropical coniferous",
                                           xlab="Cardamom", 
                                           ylab=""))

# 7) tropical subtropical moist broadleaf ----
trpsbtrp_m_broad_std_c_n <- j_trp_sbtrp_m_br_slastd$Standard_Deviation
trpsbtrp_m_broad_std_b_n <- j_trp_sbtrp_m_br_slastd$specific.leaf.area
(heatsc_trpsbtrp_con_slastd <- heatscatter(trpsbtrp_m_broad_std_c_n, 
                                           trpsbtrp_m_broad_std_b_n, 
                                           pch = 19, cexplot = 0.5, 
                                           colpal="spectral", 
                                           add.contour=TRUE, 
                                           main = "Tropical/subtropical moist broadleaf",
                                           xlab="", 
                                           ylab=""))

# 8) mediterranean forests, woodlands, scrub ----
med_f_w_scr_std_c_n <- j_med_f_slastd$Standard_Deviation
med_f_w_scr_std_b_n <- j_med_f_slastd$specific.leaf.area
(heatsc_med_f_w_scr_slastd <- heatscatter(med_f_w_scr_std_c_n, 
                                          med_f_w_scr_std_b_n, 
                                           pch = 19, cexplot = 0.5, 
                                           colpal="spectral", 
                                           add.contour=TRUE, 
                                           main = "Mediterranean woodland and scrubland",
                                           xlab="Cardamom", 
                                           ylab=""))
# dev.off ----
dev.off()




# STATS MAJ BIOMES: r2, rmse, bias ----

# first thing - renaming sla and specific.leaf.area parameters to cardamom and
# butler respectively

  # for sla mean
for (i in 1:length(j.biome.sla)){
  colnames(j.biome.sla[[i]]) <- sub("sla","cardamom",colnames(j.biome.sla[[i]]))
  colnames(j.biome.sla[[i]]) <- sub("specific.leaf.area","butler",
                                    colnames(j.biome.sla[[i]]))
}
  # for sla stdev
for (i in 1:length(j.biome.slastd)){
  colnames(j.biome.slastd[[i]]) <- sub("Standard_Deviation","cardamom",
                                       colnames(j.biome.slastd[[i]]))
  colnames(j.biome.slastd[[i]]) <- sub("specific.leaf.area","butler",
                                       colnames(j.biome.slastd[[i]]))
}

# creating function to perform calculations of stats for each biome
stats.f <- function(x){
  lapply(x,function(df){
    df$c_mean <- mean(df$cardamom)      
    df$b_mean <- mean(df$butler)
    df$diff_butler <- df$butler-df$b_mean
    df$diff_butler2<- df$diff_butler^2
    df$sum_diff_butler2 <-sum(df$diff_butler2)
    df$slope_bf <-sum((df$cardamom-df$c_mean)*(df$butler-df$b_mean))/
      sum((df$cardamom-df$c_mean)^2)
    df$b_intercept <- df$b_mean - (df$slope_bf*df$c_mean)
    df$new_b_val <- df$b_intercept + (df$slope_bf*df$cardamom)
    df$dist_mean_new_b <- df$new_b_val - df$b_mean
    df$sqrd_dist_b <- df$dist_mean_new_b^2
    df$sum_sqrd_dist_b <- sum(df$sqrd_dist_b)
    df$sla_r2 <- df$sum_sqrd_dist_b / df$sum_diff_butler2
    df$bias_av <- bias(df$butler, df$cardamom)
    df$bias_row <- df$butler - df$cardamom
    df$rmse_av <- rmse(df$butler, df$cardamom)
    df$rmse_row <- sqrt(se(df$butler, df$cardamom))
    df
  })
}

  # sla mean 
biome.sla.stats <- stats.f(j.biome.sla)

  # sla stdev
biome.slastd.stats <- stats.f(j.biome.slastd)

###############################################
# stat analysis by biome in diff continents        #
###############################################

#### splitting biomes by continent ####



#### MAKE A TABLE WITH R2 AND RMSE AVERAGES for all the different subdivisions ####
Area <- c("global", "tropics", "subtropics",
         "temperate", "pole (N)", 
         "tundra", "taiga", "temperate coniferous",
         "temperate broadleaf/mixed",
         "tropical and subtropical dry broadleaf",
         "tropical and subtropical coniferous",
         "tropical and subtropical moist broadleaf",
         "mediterranean forest, woodland and scrubland")
R2 <- c(unique(c(global_sla_stat$trps_sla_r2,
                 trps_sla_stat$trps_sla_r2, 
                 sbtrp_sla_stat$sbtrp_sla_r2,
                 tmp_sla_stat$tmp_sla_r2,
                 pl_stat$pl_sla_r2, tundra_sla_stat$sla_r2, 
                 taiga_sla_stat$sla_r2, temp_con_sla_stat$sla_r2,
                 temp_b_m_sla_stat$sla_r2, trpsbtrp_d_b_sla_stat$sla_r2,
                 trpsbtrp_con_sla_stat$sla_r2, trpsbtrp_m_b_sla_stat$sla_r2,
                 med_f_w_scr_sla_stat$sla_r2)))
RMSE <- c(unique(c(global_sla_stat$rmse_av, trps_sla_stat$rmse_av,
                   sbtrp_sla_stat$rmse_av, tmp_sla_stat$rmse_av,
                   pl_stat$rmse_av, tundra_sla_stat$rmse_av,
                   taiga_sla_stat$rmse_av, temp_con_sla_stat$rmse_av,
                   temp_b_m_sla_stat$rmse_av, trpsbtrp_d_b_sla_stat$rmse_av,
                   trpsbtrp_con_sla_stat$rmse_av, trpsbtrp_m_b_sla_stat$rmse_av,
                   med_f_w_scr_sla_stat$rmse_av)))
R2_std <-   c(unique(c(global_slastd_stat$sla_std_r2, 
                       trps_std_stat$trps_sla_std_r2, 
                       sbtrp_std_stat$sbtrp_sla_std_r2,
                       tmp_std_stat$tmp_sla_std_r2,
                       pl_std_stat$pl_sla_std_r2, tundra_slastd_stat$sla_std_r2,
                       taiga_slastd_stat$sla_std_r2, 
                       temp_con_slastd_stat$sla_std_r2,
                       temp_b_m_slastd_stat$sla_std_r2, 
                       trpsbtrp_d_b_slastd_stat$sla_std_r2,
                       trpsbtrp_con_slastd_stat$sla_std_r2,
                       trpsbtrp_m_b_slastd_stat$sla_std_r2, 
                       med_f_w_scr_slastd_stat$sla_std_r2)))
RMSE_std <- c(unique(c(global_slastd_stat$rmse_av, 
                       trps_std_stat$rmse_av, sbtrp_std_stat$rmse_av,
                       tmp_std_stat$rmse_av,
                       pl_std_stat$rmse_av, tundra_slastd_stat$rmse_av,
                       taiga_slastd_stat$rmse_av, temp_con_slastd_stat$rmse_av,
                       temp_b_m_slastd_stat$rmse_av, 
                       trpsbtrp_d_b_slastd_stat$rmse_av, 
                       trpsbtrp_con_slastd_stat$rmse_av, 
                       trpsbtrp_m_b_slastd_stat$rmse_av,
                       med_f_w_scr_slastd_stat$rmse_av)))
bias <- c(unique(c(global_sla_stat$bias, trps_sla_stat$bias,
                   sbtrp_sla_stat$bias, tmp_sla_stat$bias,
                   pl_stat$bias, tundra_sla_stat$bias, taiga_sla_stat$bias,
                   temp_con_sla_stat$bias, temp_b_m_sla_stat$bias,
                   trpsbtrp_d_b_sla_stat$bias, trpsbtrp_con_sla_stat$bias,
                   trpsbtrp_m_b_sla_stat$bias, med_f_w_scr_sla_stat$bias)))
bias_std <- c(unique(c(global_slastd_stat$bias, trps_std_stat$bias,
                       sbtrp_std_stat$bias, tmp_std_stat$bias,
                       pl_std_stat$bias, tundra_slastd_stat$bias,
                       taiga_slastd_stat$bias, temp_con_slastd_stat$bias,
                       temp_b_m_slastd_stat$bias, trpsbtrp_d_b_slastd_stat$bias,
                       trpsbtrp_con_slastd_stat$bias,
                       trpsbtrp_m_b_slastd_stat$bias, 
                       med_f_w_scr_slastd_stat$bias)))
stat_results <- as.data.frame(cbind(c(stat_results, Area)))
stat_results <- stat_results %>%
  rename("area" = V1) %>%
  filter(area != 0) %>%
  mutate(R2 = R2, R2_std = R2_std, RMSE = RMSE, RMSE_std = RMSE_std,
         bias = bias, bias_std = bias_std)

(stat_results_table <- stat_results %>%
  dplyr::rename("Sla Mean (r2)" = R2, "Sla StDev (r2)" = R2_std, 
                "Sla Mean (rmse)" = RMSE,
                "Sla StDev (rmse)" = RMSE_std,
                "Area" = area, "Sla Mean (bias)" = bias,
                "Sla StDev (bias)" = bias_std) %>%
  kable(digits = 30, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "scale_down"),
    full_width = F,
                position = "center", font_size = 10) %>%
  add_header_above(c(" ", "R2" = 2, "RMSE" = 2, "Bias" = 2), bold = T) %>%
  kableExtra::group_rows("Latitudinal gradient", 2,5) %>%
  kableExtra::group_rows("Biome",6,13) %>%
  as_image(stat_results_table, file= "./figures/table_stat_results.png", 
           width = 4, dpi = 500))






# im not sure im understanding how to fucking calculate bias and rmse for fucks sake







