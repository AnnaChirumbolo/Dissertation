###############################################################################
#                 BSc Hons Dissertation                                       #
#                 Assessing degree of consistency between                     # 
#                 functional leaf trait estimates from diff methods           #
###############################################################################

## things to do 
  
## STIPPLING!!!!!!!!

# scatter plot - make heatmap --| ASK CODING CLUB / MADE A HEATSCATTER (BETTER REPR) / ASK IF I COULD DO IT WITH THE MAP OF THE WORLD?

# for results (to put in table and have the figures showing)
# matrix of map, see how they differ from each other by each grid point 
# calculation of bias / RMSE / R2 - visual and tabular repr.

# MAKING MASKS
  # splitting the globe, by latitudinal bands to observe differences by biome   
  # splitting the globe, between different evergreen/deciduous (FACTOR SEASONALITY) - to carry out same observations at that level of comparison 

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


### opening netcdf, to data frame ----
cardamom_sla <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc", 
                       varname="sla")
cardamom_sla_std <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc", 
                           varname="Standard_Deviation")

cardamom_sla_df <- raster::as.data.frame(cardamom_sla, xy = TRUE)
cardamom_sla_std_df <- raster::as.data.frame(cardamom_sla_std, xy=TRUE)

butler_sla <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", 
                     varname="sla")
butler_sla_std <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", 
                         varname ="sla_std")
butler_sla_df <- raster::as.data.frame(butler_sla, xy = TRUE) 
butler_sla_std_df <- raster::as.data.frame(butler_sla_std, xy=TRUE)

### visualising sla from both datasets DIFF SCALES ----
png("./figures/plot_sla_DIFFSCALES.png", width = 50, height = 20, units = "cm", res = 200)
par(mfrow=c(1,2), oma = c(0,3,8,0) + 0.1, mar = c(7,0,2,8) + 0.1, new=FALSE)
plot(cardamom_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")), 
     main="Cardamom\n")
#raster::image(cardamom_sla[[1]], asp=NA, col= rev(brewer.pal(10,"RdBu")))
plot(butler_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")),
     legend.args=list(text='\nSpecific Leaf Area (m2.kg-1)', side=4, font=2, line=2.3),
     main="Butler\n")
#grid.text("Specific Leaf Area (m2.kg-1)", x=unit(0.95, "npc"), y=unit(0.50, "npc"), rot=-90)
dev.off()

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
setMinMax(cardamom_sla_std[[1]]-butler_sla_std[[1]])
breakpoints <- c(-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45, 50,
                     55,60,65)
colors <- c("#27408B", "#36648B", "#4876FF", "#8DEEEE", "#FFFFFF", 
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

### Joining datasets ----

joined_sla <- left_join(cardamom_sla_df, butler_sla_df) 
joined_sla_std <- left_join(cardamom_sla_std_df, butler_sla_std_df)

joined_sla <- joined_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area)
joined_sla_std <- joined_sla_std %>%
  rename("cardamom_std" = Standard_Deviation, "butler_std" = specific.leaf.area)

  # new datasets with removed NAs/x-y coords from the joined datasets 
joined_sla_noNA <- joined_sla %>%
  filter(cardamom!=0&butler!=0)
joined_sla_nocoord <- joined_sla %>%
  dplyr::select(-x,-y)
joined_sla_nocoordNA <- joined_sla_nocoord %>%
  filter(cardamom!=0 & butler != 0)

joined_sla_std_noNA <- joined_sla_std %>%
  filter(cardamom_std !=0 & butler_std != 0)
joined_sla_std_nocoord <- joined_sla_std %>%
  dplyr::select(-x,-y)
joined_sla_std_nocoordNA <- joined_sla_nocoord %>%
  filter(cardamom!=0 & butler != 0)


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

  # HEATSCATTERS ----
  # creating numeric vectors of sla mean and sla stdev to inpu in heatscatter ----
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

  # calc diff between HISTOGRAM areas? PANELLED PLOTS ----
  # 1) SLA:
hist_slastd_data <- as.data.table(ggplot_build(slastd_hist)$data[1])
hist_slastd_data <- hist_slastd_data %>%
  select(count, xmin, xmax, group) 
# atm for some reason, error "unused groups (count,xmin,xmax,group)" dont know why
slastd_hist_c <- hist_slastd_data[group==1]
slastd_hist_b <- hist_slastd_data[group==2]
diff_slastd_hist <- merge(slastd_hist_c, slastd_hist_b, by=c("xmin","xmax"), 
                          suffixes = c(".c",".b"), allow.cartesian=TRUE)
    #Error in vecseq(f__, len__, if (allow.cartesian || notjoin || !anyDuplicated(f__,  : 
    #Join results in 3094 rows; more than 1024 = nrow(x)+nrow(i). 
    #Check for duplicate key values in i each of which join to the same group in x over and over again. If that's ok, try by=.EACHI to run j for each group to avoid the large allocation. If you are sure you wish to proceed, rerun with allow.cartesian=TRUE. Otherwise, please search for this error message in the FAQ, Wiki, Stack Overflow and data.table issue tracker for advice.

diff_slastd_hist <- diff_slastd_hist[,Difference:=count.c-count.b]
setnames(diff_slastd_hist, old = c("count.c","count.b"),new = c("Cardamom",
                                                                "Butler"))
diff_slastd_melt <- melt(diff_slastd_hist, id.vars = c("xmin","xmax"), 
                         measure.vars = c("Cardamom","Difference","Butler"))

(diff_histstd_plot <- ggplot(diff_slastd_melt, aes(xmin=xmin, xmax=xmax, 
                                                   ymax=value,ymin=0, 
                                                   group=variable, 
                                         fill=variable, color=variable, 
                                         alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    theme(legend.title = element_blank()))  ## figure out how to show the overlap 
# the whole length it says the difference between cardamom-butler!
ggsave("./figures/Difference_hist.png", diff_dens_plot, width = 30, height = 20,
       units = "cm", dpi = 300)
    
  # save the file of difference SLA STDEV to csv
diff_hist <- write.csv(diff_dens, "difference_hist_slastd.csv")

  # 2) SLA STDEV:
hist_slastd_data <- as.data.table(ggplot_build(slastd_hist)$data[1])
hist_slastd_data <- hist_slastd_data %>%
  select(count, xmin, xmax, group)
str(hist_slastd_data)
slastd_hist_c <- hist_slastd_data[group==1]
slastd_hist_b <- hist_slastd_data[group==2]
diff_slastd_hist <- merge(slastd_hist_c, slastd_hist_b, by=c("xmin","xmax"), suffixes = c(".c",".b"), allow.cartesian=TRUE)
#Error in vecseq(f__, len__, if (allow.cartesian || notjoin || !anyDuplicated(f__,  : 
#Join results in 3094 rows; more than 1024 = nrow(x)+nrow(i). Check for duplicate key values in i each of which join to the same group in x over and over again. If that's ok, try by=.EACHI to run j for each group to avoid the large allocation. If you are sure you wish to proceed, rerun with allow.cartesian=TRUE. Otherwise, please search for this error message in the FAQ, Wiki, Stack Overflow and data.table issue tracker for advice.
diff_slastd_hist <- diff_slastd_hist[,Difference:=count.c-count.b]
setnames(diff_slastd_hist, old = c("count.c","count.b"),new = c("Cardamom",
                                                                "Butler"))
diff_slastd_hist1 <- melt(diff_slastd_hist, id.vars = c("xmin","xmax"), 
                          measure.vars = c("Cardamom","Difference","Butler"))

(diff_histstd_plot <- ggplot(diff_slastd_hist1, aes(xmin=xmin, xmax=xmax, ymax=value,ymin=0, group=variable, 
                                              fill=variable, color=variable, alpha = 0.9))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    theme(legend.title = element_blank())) # the whole length it says the difference between cardamom-butler!
ggsave("./figures/Difference_hist_slastd.png", diff_histstd_plot, width = 30, height = 20, units = "cm", 
       dpi = 300)

# save the file of difference SLA to csv
diff_hist <- write.csv(diff_dens, "difference_hist.csv")



### scatterplot: by ID (unique lat x lon point) ----

(scatter_degree <- ggplot()+
   geom_point(new_joined, mapping = aes(degree,sla, colour="Cardamom"), 
            #  colour= "#8B0000",
              alpha=0.7)+
   geom_point(new_joined, mapping = aes(degree,specific.leaf.area, colour="Butler"), 
             # colour="#EEB422", 
              alpha=0.7)+
   scale_color_manual(name="Specific Leaf Area", values= c("#8B0000", "#EEB422"))+
   theme(legend.position =c(0.9,0.8), 
         panel.background = element_rect(fill = "white", colour = "black"))+
   xlab("\nID (single lat/lon point)")+
   ylab("Specific Leaf Area (m^2/kg)\n"))

ggsave("scatter_by1x1degree.jpg", plot = last_plot(), width = 20, height = 10,
       dpi = 300)

### scatterplot by single lat and lon, panelled together ----

(scatter_bylon <- ggplot()+
   geom_point(new_joined, mapping = aes(x,sla, colour="Cardamom"), alpha=0.7)+
   geom_point(new_joined, mapping = aes(x,specific.leaf.area, colour="Butler"),
              alpha=0.7)+
   scale_color_manual(name="Specific Leaf Area", values= c("#8B0000", "#EEB422"))+
   theme(legend.position =c(0.9,0.8), 
         panel.background = element_rect(fill = "white", colour = "black"))+
   scale_y_continuous(sec.axis = ~ .)+
   xlab("\n Longitude (West to East)")+
   ylab("Specific Leaf Area (m^2/kg)\n"))


(scatter_bylat <- ggplot()+
    geom_point(new_joined, mapping = aes(sla,y.y,colour="Cardamom"), alpha=0.7)+
    geom_point(new_joined, mapping = aes(specific.leaf.area,y.y, colour="Butler"),
               alpha=0.7)+
    scale_color_manual(name="Specific Leaf Area", values= c("#8B0000", "#EEB422"))+
    theme(legend.position =c(0.9,0.8), 
          panel.background = element_rect(fill = "white", colour = "black"))+
    xlab("\nSpecific Leaf Area (m^2/kg)")+
    ylab("Latitude (North to South)\n"))

library(gridExtra)
panelled_scatter <- grid.arrange(scatter_bylat, scatter_bylon, ncol=2)

ggsave("panel_scatter_bycoord.jpg", plot = panelled_scatter, width = 20, height = 10,
       dpi = 300)


### heatmap ----

new_joined <- left_join(cardamom_df, butler_nc) %>%
  filter(sla!=0 & specific.leaf.area!=0)

#new_joined <- new_joined %>%
 # select(-degree) 

new_joined <- new_joined %>%
  rename("cardamom"=sla,"butler"=specific.leaf.area)%>%
  gather(key=dataset,value=sla,-x,-y)


joined_wide <- as.matrix(joined_wide)
rownames(joined_wide) <- joined_wide[, 1]
joined_wide <- joined_wide[, -1];
str(joined_wide)
heatmap.2(joined_wide,Colv = NA, Rowv = NA)


### COMPARING RESULT ESTIMATES ---- 

  # BIAS - BETTER ONLY FOR REGIONAL ESTIMATIONS - for later----

bias <- Metrics::bias(cardamom_sla, butler_sla) ## let's try it?? ASK CC?
bias # gives overall rel bias, and because there's lots of NAs it gives NA as result?
butler_sla <- butler_nc$sla  
butler_sla <- as.data.frame(butler_sla) 
butler_sla
cardamom_sla <- cardamom_df$sla
cardamom_sla <- as.data.frame(cardamom_sla) 



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




### R SQUARED CALCULATIONS AND CSV OUTPUT ----
# Step by step calculation of R2 value sla mean/std correlation ----
# sla mean ----
sla_r2 <- joined_sla
sla_r2 <- sla_r2 %>%
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
         sla_r2 = sum_sqrd_dist_b / sum_diff_butler2)
# resulting R2: 0.0003105955

# sla stdev ---- 

slastd_r2 <- joined_sla_std
slastd_r2 <- slastd_r2 %>%
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
         sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b)
# result for sla stdev r2: 0.0002966 

# joined sla r2 results for sla mean and sla stdev ----

r2 <- as.data.frame(c(0.0003105955, 0.0002966)) # R2 results - global
r2 <- r2 %>%
  rename("r2" = 1) %>%
  mutate(parameter = c("SLA mean", "SLA StDev"))

#write.csv(rsq_results, "R2_results.csv")

# RMSE- NEED TO SET AT SAME SCALE????----
# SLA MEAN
sla_rmse <- joined_sla %>%
  mutate(rmse = sqrt((cardamom-butler)^2)) %>%
  dplyr::select(-cardamom & -butler) %>%
  filter(rmse!=0)

# average RMSE value for MEAN SLA 8.335043
sla_rmse_all <- joined_sla %>%
  filter(cardamom !=0, butler!=0) %>%
  summarise(rmse= sqrt(mean((cardamom-butler)^2)))

# plotting rmse sla mean only 
(sla_rmse_plot <- ggplot(sla_rmse, aes(x,y,color=rmse))+
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

slastd_rmse <- joined_sla_std %>%
  mutate(rmse = sqrt((cardamom_std - butler_std)^2)) %>%
  dplyr::select(-cardamom_std, -butler_std) %>%
  filter(rmse!=0)

# average RMSE value for STDEV SLA 
slastd_rmse_all <- joined_sla_std %>%
  filter(cardamom_std!=0, butler_std!=0) %>%
  summarise(rmse= sqrt(mean((cardamom_std-butler_std)^2)))

# plotting rmse sla stdev only 
(slastd_rmse_plot <- ggplot(slastd_rmse, aes(x,y,color=rmse))+
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

# general RMSE formula
# sqrt(mean((m - o)^2)), where m is the "ref model" and o is the "observed model"

rmse_all <- rbind(sla_rmse_all, slastd_rmse_all)


## trial for other way to calculate rmse ----
# calculation of the residual - difference between estimated and actual val of the model
# how much a model disagrees with the actual data
# sla mean rmse ----

sla_rmse <- sla_r2 %>%
  mutate(residuals = butler - new_b_val, 
         resid2 = residuals^2,
         sum_resid2 = sum(resid2),
         mean_sum_resid2 = sum_resid2/(11910-1),
         rmse = sqrt(mean_sum_resid2))
# sla mean rmse: 4.274008

sla_rmse_byrow <- sla_r2 %>%
  mutate(residuals = butler - new_b_val,
         resid2 = residuals^2, 
      #   mean_resid2 = resid2 / (11910-1), a bit unsure about this? whether to divide it by the dof if there's actually no mean, because the mean it's itself in each single data point
         rmse_byrow = sqrt(resid2))

# plotting sla mean rmse ----
(sla_rmse_plot1 <- ggplot(sla_rmse_byrow, aes(x,y,color=rmse_byrow))+
   geom_jitter(stat = "identity")+
   theme_classic()+
   scale_color_gradient(low = "yellow", high = "red4",
                        limits=c(0,20))+
   ylab("Latitude\n")+
   xlab("\nLongitude")+
   labs(color=" ")+
   ggtitle("Mean SLA RMSE\n")+
   theme(plot.title = element_text(face = "bold")))

# sla stdev rmse ----
slastd_rmse <- slastd_r2 %>%
  mutate(resid_std = butler_std - new_bstd_val,
         resid_std2 = resid_std^2, 
         sum_resid_std2 = sum(resid_std2),
         mean_sum_resid_std2 = sum_resid_std2/(11910-1),
         rmse_std = sqrt(mean_sum_resid_std2))
# sla std rmse: 2.727371

# sla stdev rmse by row ----

slastd_rmse_byrow <- slastd_r2 %>%
  mutate(resid_std = butler_std - new_bstd_val, 
         resid_std2 = resid_std^2,
         rmse_std_byrow = sqrt(resid_std2))

(slastd_rmse_plot1 <- ggplot(slastd_rmse_byrow, aes(x,y,color=rmse_std_byrow))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "red4",
                         limits = c(0,20))+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold")))

## panel rmse plots and save ----
(rmse_panelled <- ggarrange(sla_rmse_plot1, slastd_rmse_plot1, ncol = 2))
ggsave("./figures/panelled_rmse.png", rmse_panelled, width = 50, height = 20,
       units = "cm", dpi = 300)

## joining rmse results ----
rmse <- as.data.frame(c(4.274008,2.727371))
rmse <- rmse %>%
  rename("rmse"=1)

### global: table with r2 and rmse ----

stat_world <- cbind(r2,rmse)
formattable(stat_world) # to create the table

## the R2 when im doing a linear regression - NOT BEING USED ATM ----
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
trp <- crop(cardamom_sla, trpext)
plot(trp[[1]]) # tropical lats
tropicsSLA_df <- raster::as.data.frame(trp, xy=TRUE)
raster::zoom(trp)
  # sla std
trpstd<- crop(cardamom_sla_std, trpext)
plot(trpstd[[1]], asp=NA) # height to fix when (and if) saving it as png+
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
trpSLA <- trpSLA %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area)
trpSTD <- left_join(tropicsSTD_df, trp_b_std_df)
trpSTD <- trpSTD %>%
  rename("cardamom_std" =Standard_Deviation, "butler_std" = specific.leaf.area)

  # subtropics ----
# cardamom:
  # sla mean 
sbtrpmatN <- matrix(data <- c(-180,23.5,180,35), nrow = 2, ncol = 2)
sbtrpextN <- extent(sbtrpmatN)  
sbtrpN <- crop(cardamom_sla, sbtrpextN)
plot(sbtrpN[[1]])
sbtrpN_df <- raster::as.data.frame(sbtrpN, xy = TRUE)
sbtrpmatS <- matrix(data <- c(-180,-23.5,180,-35), nrow = 2, ncol = 2)
sbtrpextS <- extent(sbtrpmatS)
sbtrpS <- crop(cardamom_sla, sbtrpextS)
plot(sbtrpS[[1]])
sbtrpS_df <- raster::as.data.frame(sbtrpS, xy = TRUE)

# merging two datasets from two hemispheres to have them in one dataframe
sbtrp_c <- merge(sbtrpN_df, sbtrpS_df, by=c("x", "y", "sla"), all=TRUE) 
# this works 
sbtrp_c <- sbtrp_c %>% rename("cardamom" = sla)
  # plotting the two strips of latitude of subtropics
(trial <- ggplot(sbtrpSTD_b, aes(x,y,color=butler_std))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("SUBTROPICS Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))

  # sla std 
sbtrp_std_N <- crop(cardamom_sla_std, sbtrpextN)
sbtrpSTD_c_N <- raster::as.data.frame(sbtrp_std_N, xy=TRUE)
sbtrp_std_S <- crop(cardamom_sla_std, sbtrpextS)
sbtrpSTD_c_S <- raster::as.data.frame(sbtrp_std_S, xy=TRUE)
sbtrpSTD_c <- merge(sbtrpSTD_c_N, sbtrpSTD_c_S, 
                    by=c("x","y","Standard_Deviation"), all = TRUE)
sbtrpSTD_c <- sbtrpSTD_c %>%
  rename("cardamom_std" = Standard_Deviation )

# butler:
  # sla mean
sbtrpN_b <- crop(butler_sla, sbtrpextN)  
plot(sbtrpN_b[[1]])
sbtrpS_b <- crop(butler_sla, sbtrpextS)
plot(sbtrpS_b[[1]])
sbtrp_bN_df <- raster::as.data.frame(sbtrpN_b, xy=TRUE)
sbtrp_bS_df <- raster::as.data.frame(sbtrpS_b, xy =TRUE)
sbtrp_b <- merge(sbtrp_bN_df, sbtrp_bS_df, by=c("x","y", "specific.leaf.area"),
                 all = TRUE)
sbtrp_b <- sbtrp_b %>%
  rename("butler"=specific.leaf.area)

  # sla std
sbtrpSTD_b_N <- crop(butler_sla_std, sbtrpextN)
sbtrpSTD_b_S <- crop(butler_sla_std, sbtrpextS)
sbtrpSTD_b_N_df <- raster::as.data.frame(sbtrpSTD_b_N, xy = TRUE)
sbtrpSTD_b_S_df <- raster::as.data.frame(sbtrpSTD_b_S, xy =TRUE)
sbtrpSTD_b <- merge(sbtrpSTD_b_N_df, sbtrpSTD_b_S_df, 
                    by = c("x","y","specific.leaf.area"),
                    all = TRUE)
sbtrpSTD_b <- sbtrpSTD_b %>%
  rename("butler_std"= specific.leaf.area)

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
tmpextS <- extent(tmpmat2)
tmpN <- crop(cardamom_sla, tmpextN)
tmpS <- crop(cardamom_sla, tmpextS)
plot(tmpN[[1]])
plot(tmpS[[1]])

tmp_c_N_df <- raster::as.data.frame(tmpN, xy=TRUE)
tmp_c_S_df <- raster::as.data.frame(tmpS, xy = TRUE)
tmp_c_df <- merge(tmp_c_N_df, tmp_c_S_df, by = c("x","y","sla"), all = TRUE)
tmp_c_df <- tmp_c_df %>% rename("cardamom"= sla)
  
  # sla std 
tmpSTD_c_N <- crop(cardamom_sla_std, tmpextN)
tmpSTD_c_S <- crop(cardamom_sla_std, tmpextS)
tmpSTD_c_N_df <- raster::as.data.frame(tmpSTD_c_N, xy =TRUE)
tmpSTD_c_S_df <- raster::as.data.frame(tmpSTD_c_S, xy =TRUE)
tmpSTD_c_df <- merge(tmpSTD_c_N_df, tmpSTD_c_S_df, 
                     by = c("x","y","Standard_Deviation"), 
                     all = TRUE)
tmpSTD_c_df <- tmpSTD_c_df %>% rename("cardamom_std" = Standard_Deviation)

# bulter 
  # sla mean 
tmp_b_N <- crop(butler_sla, tmpextN)
tmp_b_S <- crop(butler_sla, tmpextS)
tmp_b_N_df <- raster::as.data.frame(tmp_b_N, xy = TRUE)
tmp_b_S_df <- raster::as.data.frame(tmp_b_S, xy = TRUE)
tmp_b_df <- merge(tmp_b_N_df, tmp_b_S_df, by = c("x","y","specific.leaf.area"),
                  all = TRUE)
tmp_b_df <- tmp_b_df %>%
  rename("butler" = specific.leaf.area)

  # sla std 
tmpSTD_b_N <- crop(butler_sla_std, tmpextN)
tmpSTD_b_S <- crop(butler_sla_std, tmpextS)
tmpSTD_b_N_df <- raster::as.data.frame(tmpSTD_b_N, xy = TRUE)
tmpSTD_b_S_df <- raster::as.data.frame(tmpSTD_b_S, xy = TRUE)
tmpSTD_b_df <- merge(tmpSTD_b_N_df, tmpSTD_b_S_df, 
                     by = c("x","y","specific.leaf.area"), all = TRUE)
tmpSTD_b_df <- tmpSTD_b_df %>%
  rename("butler_std" = specific.leaf.area)

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
plot(plN_c[[1]])
plN_c_df <- raster::as.data.frame(plN_c, xy= TRUE)
plN_c_df <- plN_c_df %>% rename("cardamom" = sla)

# sla stdev
plN_slastd_c <- crop(cardamom_sla_std, plextN)
plot(plN_slastd_c[[1]])
plN_slastd_c_df <- raster::as.data.frame(plN_slastd_c, xy = TRUE)
plN_slastd_c_df <- plN_slastd_c_df %>% 
  rename("cardamom_std" = Standard_Deviation)

# butler
# sla mean
plN_b <- crop(butler_sla, plextN)
plot(plN_b[[1]])
plN_b_df <- raster::as.data.frame(plN_b, xy = TRUE)
plN_b_df <- plN_b_df %>% rename("butler"=specific.leaf.area)

# sla stdev
plN_slastd_b <- crop(butler_sla_std, plextN)
plot(plN_slastd_b[[1]])
plN_slastd_b_df <- raster::as.data.frame(plN_slastd_b, xy = TRUE)
plN_slastd_b_df <- plN_slastd_b_df %>% 
  rename("butler_std"= specific.leaf.area)

# joining butler and cardamom for poles lat
plN <- left_join(plN_c_df, plN_b_df)
plN_std <- left_join(plN_slastd_c_df, plN_slastd_b_df)



#### ANALYSIS TROPICS ####
# heatscatter tropics ----
# sla mean
trp_c_sla_n <- trpSLA$cardamom
trp_b_sla_n <- trpSLA$butler

png("./figures/heatscatter_trps_sla.png", width = 30, height = 15, units = "cm",
    res = 300)
(heatsc_sla_trp <- heatscatter(trp_c_sla_n, trp_b_sla_n, pch = 19, 
                                cexplot = 0.5, colpal="spectral", 
                                #disco() for all color options / could set to colorblind
                                add.contour=TRUE, main = "Tropics SLA Mean\n",
                                xlab="\nCardamom", 
                                ylab="\nButler"))
dev.off()

#corr_trp_sla <- chart.Correlation(trpSLA)

# sla stdev 
trp_c_std_n <- trpSTD$cardamom_std
trp_b_std_n <- trpSTD$butler_std

png("./figures/heatscatter_trps_slastd.png", width = 30, height = 15, 
    units = "cm", res = 300)
(heatsc_slastd_trp <- heatscatter(trp_c_std_n, trp_b_std_n, pch = 19, 
                               cexplot = 0.5, colpal="spectral", 
                               #disco() for all color options / could set to colorblind
                               add.contour=TRUE, main = "Tropics SLA StDev\n",
                               xlab="\nCardamom", 
                               ylab="\nButler"))
dev.off()

## stats: R2 ----
# sla mean
trps_r2 <- trpSLA
trps_r2 <- trps_r2 %>%
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
         trps_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2)
# Multiple R-squared:  0.1271

# sla stdev
trps_std_r2 <- trpSTD
trps_std_r2 <- trps_std_r2 %>%
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
         trps_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b)
# R2 stdev tropics: 0.0380433

## stats: RMSE ----
# sla mean 
  # rmse for each data point
trp_sla_rmse <- trpSLA %>%
  mutate(rmse = sqrt((cardamom-butler)^2)) %>%
  dplyr::select(-cardamom & -butler) %>%
  filter(rmse!=0)

  # plotting rmse for each data point 
(trps_sla_rmse_plot <- ggplot(trp_sla_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TROPICS Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))
ggsave("./figures/trps_sla_rmse.png", trps_sla_rmse_plot, width = 30, 
       height = 10, units = "cm", dpi = 300)

  # average rmse sla mean: 9.464195
trps_sla_rmse_av <- trpSLA %>%
  filter(cardamom !=0, butler!=0) %>%
  summarise(rmse= sqrt(mean((cardamom-butler)^2)))

# sla stdev 
  # rmse for each point 
trp_slastd_rmse <- trpSTD %>%
  mutate(rmse = sqrt((cardamom_std-butler_std)^2)) %>%
  dplyr::select(-cardamom_std & -butler_std) %>%
  filter(rmse!=0)
# plotting rmse for each data point 
(trp_slastd_rmse_plot <- ggplot(trp_slastd_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TROPICS SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold")))
ggsave("./figures/trps_slastd_rmse.png", trp_slastd_rmse_plot, width = 30, 
       height = 10, units = "cm", dpi = 300)

# average rmse sla stdev: 40.99694
trps_slastd_rmse_av <- trpSTD %>%
  filter(cardamom_std !=0, butler_std!=0) %>%
  summarise(rmse= sqrt(mean((cardamom_std-butler_std)^2)))


#### ANALYSIS SUBTROPICS ####
# heatscatter subtrps ----
# sla mean
sbtrp_c_sla_n <- sbtrp_joined_SLA$cardamom
sbtrp_b_sla_n <- sbtrp_joined_SLA$butler

(heatsc_sla_sbtrp <- heatscatter(sbtrp_c_sla_n, sbtrp_b_sla_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Subtropics SLA Mean\n",
                                  xlab="\nCardamom", 
                                  ylab="\nButler"))
#sla stdev 
sbtrp_c_slastd_n <- sbtrp_joined_STD$cardamom_std
sbtrp_b_slastd_n <- sbtrp_joined_STD$butler_std

(heatsc_slastd_sbtrp <- heatscatter(sbtrp_c_slastd_n, sbtrp_b_slastd_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Subtropics SLA StDev\n",
                                  xlab="\nCardamom", 
                                  ylab="\nButler"))
# need to save them panelled 
## stats: R2 ----
# sla mean r2
sbtrp_r2 <- sbtrp_joined_SLA %>%
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
         sbtrp_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2)
# r2 result for mean sla in subtropics: 0.1372143

# sla stdev r2 
sbtrp_std_r2 <- sbtrp_joined_STD
sbtrp_std_r2 <- sbtrp_std_r2 %>%
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
         sbtrp_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b)
# r2 result for sla stdev in subtropics: 0.01140631
lm_trial <- lm(butler_std  ~ cardamom_std, data = sbtrp_joined_STD)
summary(lm_trial) # same result as manual calculation 

## stats: RMSE ----
# rmse sla mean 
  # for each single point
sbtrp_sla_rmse <- sbtrp_joined_SLA %>%
  mutate(rmse = sqrt((cardamom-butler)^2)) %>%
  dplyr::select(-cardamom & -butler) %>%
  filter(rmse!=0)
  # plotting rmse for each data point
(sbtrp_sla_rmse_plot <- ggplot(sbtrp_sla_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("SUBTROPICS Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold"))) # this is still showing me only the northern hemisphere...
ggsave("./figures/sbtrp_sla_rmse.png", sbtrp_sla_rmse_plot, width = 30, 
       height = 15, units = "cm", dpi = 300)

# average rmse sla mean subtropics: 7.340305
subtrp_rmse_av_sla <- sbtrp_joined_SLA %>%
  filter(cardamom !=0, butler!=0) %>%
  summarise(rmse= sqrt(mean((cardamom-butler)^2)))

# sla stdev 
  # rmse each point
sbtrp_slastd_rmse <- sbtrp_joined_STD %>%
  mutate(rmse = sqrt((cardamom_std-butler_std)^2)) %>%
  dplyr::select(-cardamom_std & -butler_std) %>%
  filter(rmse!=0)
# plotting rmse for each data point 
(sbtrp_slastd_rmse_plot <- ggplot(sbtrp_slastd_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("SUBTROPICS SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold")))
ggsave("./figures/sbtrp_slastd_rmse.png", sbtrp_slastd_rmse_plot, width = 30, 
       height = 10, units = "cm", dpi = 300)

# average rmse sla stdev: 40.03824
sbtrp_slastd_rmse_av <- sbtrp_joined_STD %>%
  filter(cardamom_std !=0, butler_std!=0) %>%
  summarise(rmse= sqrt(mean((cardamom_std-butler_std)^2)))


#### ANALYSIS TEMPERATE ####
# heatscatter ----
  #sla mean
tmp_sla_c_n <- tmp_sla$cardamom
tmp_sla_b_n <- tmp_sla$butler
(heatsc_sla_tmp <- heatscatter(tmp_sla_c_n, tmp_sla_b_n, pch = 19, 
                               cexplot = 0.5, colpal="spectral", 
                               #disco() for all color options / could set to colorblind
                               add.contour=TRUE, main = "Temperate SLA Mean\n",
                               xlab="\nCardamom", 
                               ylab="\nButler"))
  # sla stdev
tmp_slastd_c_n <- tmp_slastd$cardamom_std
tmp_slastd_b_n <- tmp_slastd$butler_std
(heatsc_slastd_tmp <- heatscatter(tmp_slastd_c_n, tmp_slastd_b_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Temperate SLA StDev\n",
                                  xlab="\nCardamom", 
                                  ylab="\nButler"))
## stats: R2 ----
# sla mean
tmp_r2 <- tmp_sla %>%
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
         tmp_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2)
# result of r2 for temperate sla mean: 0.00640377

# sla stdev
tmp_std_r2 <- tmp_slastd
tmp_std_r2 <- tmp_std_r2 %>%
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
         tmp_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b)
# results for r2 in temperate for sla stdev: 0.0001560524

## stats: RMSE ----
# sla mean 
  # each data point
tmp_sla_rmse <- tmp_sla %>%
  mutate(rmse = sqrt((cardamom-butler)^2)) %>%
  dplyr::select(-cardamom & -butler) %>%
  filter(rmse!=0)
# plotting rmse for each data point
(tmp_sla_rmse_plot <- ggplot(tmp_sla_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TEMPERATE Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))
ggsave("./figures/tmp_sla_rmse.png", tmp_sla_rmse_plot, width = 30, 
       height = 15, units = "cm", dpi = 300)
  # RMSE average data points: 7.215715
tmp_rmse_av_sla <- tmp_sla %>%
  filter(cardamom !=0, butler!=0) %>%
  summarise(rmse= sqrt(mean((cardamom-butler)^2)))
# sla stdev 
  # each data point 
tmp_slastd_rmse <- tmp_slastd %>%
  mutate(rmse = sqrt((cardamom_std-butler_std)^2)) %>%
  dplyr::select(-cardamom_std & -butler_std) %>%
  filter(rmse!=0)
# plotting rmse for each data point
(tmp_slastd_rmse_plot <- ggplot(tmp_slastd_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TEMPERATE SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold"))) # this is still showing me only the northern hemisphere...
ggsave("./figures/tmp_slastd_rmse.png", tmp_slastd_rmse_plot, width = 30, 
       height = 15, units = "cm", dpi = 300)
  # RMSE average data points: 29.67961
tmp_rmse_av_slastd <- tmp_slastd %>%
  filter(cardamom_std !=0, butler_std!=0) %>%
  summarise(rmse= sqrt(mean((cardamom_std-butler_std)^2)))


#### ANALYSIS POLES ####
## stats: R2 ----
# sla mean
pl_r2 <- plN %>%
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
         pl_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2)
# result of r2 for temperate sla mean: 0.02279148

# sla stdev
pl_std_r2 <- plN_std %>%
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
         pl_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b)
# results for r2 in temperate for sla stdev: 0.001805205
## stats: RMSE ----
# sla mean 
  # each data point
pl_sla_rmse <- plN %>%
  mutate(rmse = sqrt((cardamom-butler)^2)) %>%
  dplyr::select(-cardamom & -butler) %>%
  filter(rmse!=0)
# plotting rmse for each data point
(pl_sla_rmse_plot <- ggplot(pl_sla_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("N POLE Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold"))) # this is still showing me only the northern hemisphere...
ggsave("./figures/Npl_sla_rmse.png", pl_sla_rmse_plot, width = 20, 
       height = 10, units = "cm", dpi = 300)
# RMSE average data points: 10.76312
pl_rmse_av_sla <- plN %>%
  filter(cardamom !=0, butler!=0) %>%
  summarise(rmse= sqrt(mean((cardamom-butler)^2)))

# sla stdev 
  # each data point 
pl_slastd_rmse <- plN_std %>%
  mutate(rmse = sqrt((cardamom_std-butler_std)^2)) %>%
  dplyr::select(-cardamom_std & -butler_std) %>%
  filter(rmse!=0)
  # plotting rmse for each data point
(pl_slastd_rmse_plot <- ggplot(pl_slastd_rmse, aes(x,y,color=rmse))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("N POLE SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold"))) # this is still showing me only the northern hemisphere...
ggsave("./figures/Npl_slastd_rmse.png", pl_slastd_rmse_plot, width = 20, 
       height = 10, units = "cm", dpi = 300)

# RMSE average data points: 25.24611
pl_rmse_av_slastd <- plN_std %>%
  filter(cardamom_std !=0, butler_std!=0) %>%
  summarise(rmse= sqrt(mean((cardamom_std-butler_std)^2)))



#### ADDITIONAL STATS: BIAS #### not sure how to to it... 
# trial - tropics ----
trps_r2 <- trps_r2 %>%
  mutate(bias_point = new_b_val - butler,
         sum_error_b = sum(bias_point),
         bias_av = sum_error_b / 3403) # very close to 0 --> 1.137954e-16

(trps_bias_sla <- ggplot(trps_r2, aes(x,y,color=bias_point))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_viridis(direction = 1)+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("Trial bias tropics sla mean\n")+
    theme(plot.title = element_text(face = "bold")))
# need to put all results in a single dataframe - incl global and by lat 
   # not sure if this is right


cold_deciduous_forest <-'./DATA/cold.deciduous.forest_p_1km_s0..0cm_2000..2017_v0.1.tif' 
raster_decidious_forest =raster(cold_deciduous_forest)
plot(imported_raster[[1]])

cold_evergreen_needleleaf <- "./DATA/cold.evergreen.needleleaf.forest_p_1km_s0..0cm_2000..2017_v0.1.tif"
raster_evergreen_needleleaf <- raster(cold_evergreen_needleleaf)
plot(raster_evergreen_needleleaf[[1]])

#### OTHER ADDITIONAL STATS: T-TEST AND F-TEST? ####

#############################
#### SPLITTING BY BIOME #####
#############################
# OPEN DATASET: ECOREGIONS17 - where i get the biome masks from ----
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

# mask the raster taiga biome  ----
# sla mean
masked_taiga_sla_c <- raster::mask(cardamom_sla, boreal_f_taiga)
masked_taiga_sla_b <- raster::mask(butler_sla, boreal_f_taiga)
plot(masked_taiga_sla_c[[1]])
plot(masked_taiga_sla_b[[1]])

# sla stdev 
masked_taiga_slastd_c <- raster::mask(cardamom_sla_std, boreal_f_taiga)
masked_taiga_slastd_b <- raster::mask(butler_sla_std, boreal_f_taiga)
plot(masked_taiga_slastd_c[[1]])
plot(masked_taiga_slastd_b[[1]])

# mask the raster tundra biome ----
# sla mean 
masked_tundra_sla_c <- raster::mask(cardamom_sla, tundra)
masked_tundra_sla_b <- raster::mask(butler_sla, tundra)
plot(masked_tundra_sla_c[[1]])
plot(masked_tundra_sla_b[[1]])

# sla stdev 
masked_tundra_slastd_c <- raster::mask(cardamom_sla_std, tundra)
masked_tundra_slastd_b <- raster::mask(butler_sla_std, tundra)
plot(masked_tundra_slastd_c[[1]])
plot(masked_tundra_slastd_b[[1]])

# mask the raster temp conif forest biome ----
# sla mean 
mask_temp_conif_sla_c <- mask(cardamom_sla,temp_conif_forest)
mask_temp_conif_sla_b <- mask(butler_sla, temp_conif_forest)
plot(mask_temp_conif_sla_c[[1]])
plot(mask_temp_conif_sla_b[[1]])

# sla stdev 
mask_temp_conif_slastd_c <- raster::mask(cardamom_sla_std, temp_conif_forest)
mask_temp_conif_slastd_b <- raster::mask(butler_sla_std, temp_conif_forest)
plot(mask_temp_conif_slastd_c[[1]])
plot(mask_temp_conif_slastd_b[[1]])

# mask the raster temperate broad and mixed forest biome ----
# sla mean
mask_temp_broad_mix_sla_c <- raster::mask(cardamom_sla, temp_broad_mix)
mask_temp_broad_mix_sla_b <- raster::mask(butler_sla, temp_broad_mix)
plot(mask_temp_broad_mix_sla_c[[1]])
plot(mask_temp_broad_mix_sla_b[[1]])

# sla stdev 
mask_temp_broad_mix_slastd_c <- raster::mask(cardamom_sla_std, temp_broad_mix)
mask_temp_broad_mix_slastd_b <- raster::mask(butler_sla_std, temp_broad_mix)
plot(mask_temp_broad_mix_slastd_c[[1]])
plot(mask_temp_broad_mix_slastd_b[[1]])

# mask the raster tropical and subtropical dry broadleaf biome ----
# sla mean
mask_trp_sbtrp_dry_broad_sla_c <- raster::mask(cardamom_sla, 
                                               trp_sbtrp_dry_broad)
mask_trp_sbtrp_dry_broad_sla_b <-raster::mask(butler_sla, trp_sbtrp_dry_broad)
plot(mask_trp_sbtrp_dry_broad_sla_c[[1]])
plot(mask_trp_sbtrp_dry_broad_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_dry_broad_slastd_c <- raster::mask(cardamom_sla_std, 
                                                  trp_sbtrp_dry_broad)
mask_trp_sbtrp_dry_broad_slastd_b <- raster::mask(butler_sla_std, 
                                                  trp_sbtrp_dry_broad)
plot(mask_trp_sbtrp_dry_broad_slastd_c[[1]])
plot(mask_trp_sbtrp_dry_broad_slastd_b[[1]])

# mask the raster tropical and subtropical conif forest ----
# sla mean 
mask_trp_sbtrp_conif_sla_c <- raster::mask(cardamom_sla, trp_sbtrp_conif)
mask_trp_sbtrp_conif_sla_b <- raster::mask(butler_sla, trp_sbtrp_conif)
plot(mask_trp_sbtrp_conif_sla_c[[1]])
plot(mask_trp_sbtrp_conif_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_conif_slastd_c <- raster::mask(cardamom_sla_std, trp_sbtrp_conif)
mask_trp_sbtrp_conif_slastd_b <- raster::mask(butler_sla_std, trp_sbtrp_conif)
plot(mask_trp_sbtrp_conif_slastd_c[[1]])
plot(mask_trp_sbtrp_conif_slastd_b[[1]])

# mask the raster tropical subtropical moist broadleaf biome ----
# sla mean
mask_trp_sbtrp_moist_broad_sla_c <- raster::mask(cardamom_sla, 
                                                 trp_sbtrp_moist_broad)
mask_trp_sbtrp_moist_broad_sla_b <- raster::mask(butler_sla, 
                                                 trp_sbtrp_moist_broad)
plot(mask_trp_sbtrp_moist_broad_sla_c[[1]])
plot(mask_trp_sbtrp_moist_broad_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_moist_broad_slastd_c <- raster::mask(cardamom_sla_std, 
                                                    trp_sbtrp_moist_broad)
mask_trp_sbtrp_moist_broad_slastd_b <- raster::mask(butler_sla_std, 
                                                    trp_sbtrp_moist_broad)
plot(mask_trp_sbtrp_moist_broad_slastd_c[[1]])
plot(mask_trp_sbtrp_moist_broad_slastd_b[[1]])

# mask the raster mediterranean forests, woodlands, scrub biome ----
# sla mean
mask_med_f_w_scr_sla_c <- raster::mask(cardamom_sla, med_f_w_scr)
mask_med_f_w_scr_sla_b <- raster::mask(butler_sla, med_f_w_scr)
plot(mask_med_f_w_scr_sla_c[[1]])
plot(mask_med_f_w_scr_sla_b[[1]])

# sla stdev 
mask_med_f_w_scr_slastd_c <- raster::mask(cardamom_sla_std, med_f_w_scr)
mask_med_f_w_scr_slastd_b <- raster::mask(butler_sla_std, med_f_w_scr)
plot(mask_med_f_w_scr_slastd_c[[1]])
plot(mask_med_f_w_scr_slastd_b[[1]])

# mask the raster desertic and xeric scrubland biome ----
# sla mean
mask_des_x_scr_sla_c <- raster::mask(cardamom_sla, des_x_scr)
mask_des_x_scr_sla_b <- raster:: mask(butler_sla, des_x_scr)
plot(mask_des_x_scr_sla_c[[1]])
plot(mask_des_x_scr_sla_b[[1]])

# sla stdev 
mask_des_x_scr_slastd_c <- raster::mask(cardamom_sla_std, des_x_scr)
mask_des_x_scr_slastd_b <- raster::mask(butler_sla_std, des_x_scr)
plot(mask_des_x_scr_slastd_c[[1]])
plot(mask_des_x_scr_slastd_b[[1]])

# mask the raster temperate grassland, savanna, shrubland biome ----
# sla mean
mask_temp_grass_sav_shr_sla_c <- raster::mask(cardamom_sla,temp_grass_sav_shr)
mask_temp_grass_sav_shr_sla_b <- raster::mask(butler_sla, temp_grass_sav_shr)
plot(mask_temp_grass_sav_shr_sla_c[[1]])
plot(mask_temp_grass_sav_shr_sla_b[[1]])

# sla stdev
mask_temp_grass_sav_shr_slastd_c <- raster::mask(cardamom_sla_std, 
                                                 temp_grass_sav_shr)
mask_temp_grass_sav_shr_slastd_b <- raster::mask(butler_sla_std, 
                                                 temp_grass_sav_shr)
plot(mask_temp_grass_sav_shr_slastd_c[[1]])
plot(mask_temp_grass_sav_shr_slastd_b[[1]])

# mask the raster montane grassland and shrubland biome ----
# sla mean 
mask_mont_grass_shr_sla_c <- raster::mask(cardamom_sla, mont_grass_shr)
mask_mont_grass_shr_sla_b <- raster::mask(butler_sla, mont_grass_shr)
plot(mask_mont_grass_shr_sla_c[[1]])
plot(mask_mont_grass_shr_sla_b[[1]])

# sla stdev 
mask_mont_grass_shr_slastd_c <- raster::mask(cardamom_sla_std, 
                                             mont_grass_shr)
mask_mont_grass_shr_slastd_b <- raster::mask(butler_sla_std, 
                                             mont_grass_shr)
plot(mask_mont_grass_shr_slastd_c[[1]])
plot(mask_mont_grass_shr_slastd_b[[1]])

# mask the raster mangrove biome ----
# sla mean
mask_mangroves_sla_c <- raster::mask(cardamom_sla, mangroves)
mask_mangroves_sla_b <- raster::mask(butler_sla, mangroves)
plot(mask_mangroves_sla_c[[1]])
plot(mask_mangroves_sla_b[[1]])

# sla stdev 
mask_mangroves_slastd_c <- raster::mask(cardamom_sla_std, mangroves)
mask_mangroves_slastd_b <- raster::mask(butler_sla_std, mangroves)
plot(mask_mangroves_slastd_c[[1]])
plot(mask_mangroves_slastd_b[[1]])

# mask raster flooded grassland and savanna biome ----
# sla mean
mask_flo_grass_sav_sla_c <- raster::mask(cardamom_sla, flo_grass_sav)
mask_flo_grass_sav_sla_b <- raster::mask(butler_sla, flo_grass_sav)
plot(mask_flo_grass_sav_sla_c[[1]])
plot(mask_flo_grass_sav_sla_b[[1]])

# sla stdev 
mask_flo_grass_sav_slastd_c <- raster::mask(cardamom_sla_std, flo_grass_sav)
mask_flo_grass_sav_slastd_b <- raster::mask(butler_sla_std, flo_grass_sav)
plot(mask_flo_grass_sav_slastd_c[[1]])
plot(mask_flo_grass_sav_slastd_b[[1]])

# mask raster tropical and subtropical grassland, savanna, shrubland biome ----
# sla mean 
mask_trp_sbtrp_grass_sav_shr_sla_c <- raster::mask(cardamom_sla, 
                                                   trp_sbtrp_grass_sav_shr)
mask_trp_sbtrp_grass_sav_shr_sla_b <- raster::mask(butler_sla, 
                                                   trp_sbtrp_grass_sav_shr)
plot(mask_trp_sbtrp_grass_sav_shr_sla_c[[1]])
plot(mask_trp_sbtrp_grass_sav_shr_sla_b[[1]])

# sla stdev 
mask_trp_sbtrp_grass_sav_shr_slastd_c <- raster::mask(cardamom_sla_std,
                                                      trp_sbtrp_grass_sav_shr)
mask_trp_sbtrp_grass_sav_shr_slastd_b <- raster::mask(butler_sla_std,
                                                      trp_sbtrp_grass_sav_shr)
plot(mask_trp_sbtrp_grass_sav_shr_slastd_c[[1]])
plot(mask_trp_sbtrp_grass_sav_shr_slastd_b[[1]])







