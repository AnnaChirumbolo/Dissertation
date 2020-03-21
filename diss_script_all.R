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
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)),
         bias = bias(butler, new_b_val))

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
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)),
         bias = bias(butler_std, new_bstd_val))

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
trp <- crop(cardamom_sla, trpext)
plot(trp[[1]]) # tropical lats
tropicsSLA_df <- raster::as.data.frame(trp, xy=TRUE)
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



# (heatscatter) saving figure for sla mean by latitudinal range----
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

## stats: r2, rmse, bias ----
# sla mean
trps_sla_stat <- trpSLA %>%
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
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)),
         bias = bias(butler, new_b_val))
# Multiple R-squared:  0.1271
# plotting rmse for each data point: sla mean
(trps_sla_rmse_plot <- ggplot(trps_sla_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TROPICS Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))

# sla stdev
trps_std_stat <- trpSTD %>%
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
         trps_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b,
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)),
         bias = bias(butler_std, new_bstd_val))
# R2 stdev tropics: 0.0380433

# plotting rmse for each data point: stdev
(trp_slastd_rmse_plot <- ggplot(trps_std_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TROPICS SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold")))


#### ANALYSIS SUBTROPICS ####
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
## stats: r2, rmse, bias ----
# sla mean 
sbtrp_sla_stat <- sbtrp_joined_SLA %>%
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
         sbtrp_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2,
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)),
         bias = bias(butler, new_b_val))
# r2 result for mean sla in subtropics: 0.1372143
# plotting rmse for each data point sla mean
(sbtrp_sla_rmse_plot <- ggplot(sbtrp_sla_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("SUBTROPICS Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))

# sla stdev 
sbtrp_std_stat <- sbtrp_joined_STD %>%
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
         sbtrp_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b,
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)),
         bias = bias(butler_std, new_bstd_val))
# r2 result for sla stdev in subtropics: 0.01140631

# plotting rmse for each data point stdev
(sbtrp_slastd_rmse_plot <- ggplot(sbtrp_std_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("SUBTROPICS SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold")))



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
## stats: R2 ----
# sla mean
tmp_sla_stat <- tmp_sla %>%
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
         tmp_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2,
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)),
         bias = bias(butler, new_b_val))
# result of r2 for temperate sla mean: 0.00640377
# plotting rmse for each data point: sla mean
(tmp_sla_rmse_plot <- ggplot(tmp_sla_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TEMPERATE Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))

# sla stdev
tmp_std_stat <- tmp_slastd %>%
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
         tmp_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b,
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)),
         bias = bias(butler_std, new_bstd_val))
# results for r2 in temperate for sla stdev: 0.0001560524

# plotting rmse for each data point
(tmp_slastd_rmse_plot <- ggplot(tmp_std_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("TEMPERATE SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold")))



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
## stats r2, rmse, bias ----
# sla mean
pl_stat <- plN %>%
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
         pl_sla_r2 = sum_sqrd_dist_b / sum_diff_butler2,
         rmse_av = Metrics::rmse(butler, new_b_val), # this is the mean of the mean values
         rmse_row = sqrt(Metrics::se(butler, new_b_val)), # square root of squared error - each row is already a mean value
         bias = Metrics::bias(butler, new_b_val)) 
# plotting rmse for each data point: sla mean
(pl_sla_rmse_plot <- ggplot(pl_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("N POLE Mean SLA RMSE\n")+
    theme(plot.title = element_text(face = "bold")))

# sla stdev
pl_std_stat <- plN_std %>%
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
         pl_sla_std_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b,
         rmse_av = Metrics::rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(Metrics::se(butler_std, new_bstd_val)),
         bias = bias(butler_std, new_bstd_val))

  # plotting rmse for each data point: stdev 
(pl_slastd_rmse_plot <- ggplot(pl_std_stat, aes(x,y,color=rmse_row))+
    geom_jitter(stat = "identity")+
    theme_classic()+
    scale_color_gradient(low = "yellow", high = "darkred")+
    ylab("Latitude\n")+
    xlab("\nLongitude")+
    labs(color=" ")+
    ggtitle("N POLE SLA StDev RMSE\n")+
    theme(plot.title = element_text(face = "bold")))

# dev.off----
dev.off()


#### OTHER ADDITIONAL STATS: T-TEST AND F-TEST? ####

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
plot(masked_taiga_sla_c[[1]])
plot(masked_taiga_sla_b[[1]])

# sla stdev 
masked_taiga_slastd_c <- raster::mask(cardamom_sla_std, boreal_f_taiga)
masked_taiga_slastd_b <- raster::mask(butler_sla_std, boreal_f_taiga)
plot(masked_taiga_slastd_c[[1]])
plot(masked_taiga_slastd_b[[1]])

# 2) mask raster tundra biome ----
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

# 3) mask raster temp conif forest biome ----
# sla mean 
mask_temp_conif_sla_c <- raster::mask(cardamom_sla,temp_conif_forest)
mask_temp_conif_sla_b <- raster::mask(butler_sla, temp_conif_forest)
plot(mask_temp_conif_sla_c[[1]])
plot(mask_temp_conif_sla_b[[1]])

# sla stdev 
mask_temp_conif_slastd_c <- raster::mask(cardamom_sla_std, temp_conif_forest)
mask_temp_conif_slastd_b <- raster::mask(butler_sla_std, temp_conif_forest)
plot(mask_temp_conif_slastd_c[[1]])
plot(mask_temp_conif_slastd_b[[1]])

# 4) mask raster temperate broad and mixed forest biome ----
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

# 5) mask raster tropical and subtropical dry broadleaf biome ----
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

# 6) mask raster tropical and subtropical conif forest biome ----
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

# 7) mask raster tropical subtropical moist broadleaf biome ----
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

### these might be extras ###
# 8) mask raster mediterranean forests, woodlands, scrub biome ----
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

# 9) mask raster desertic and xeric scrubland biome ----
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

# 10) mask raster temperate grassland, savanna, shrubland biome ----
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

# 11) mask raster montane grassland and shrubland biome ----
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

# 12) mask raster mangrove biome ----
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

# 13) mask raster flooded grassland and savanna biome ----
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

# 14) mask raster tropical and subtropical grassland, savanna, shrubland biome ----
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




#### VISUAL AND STAT ANALYSIS BY BIOME ####
# first thing - turning all masked rasterlayers to dataframes ----
  # sla mean
# 1) taiga
taiga_sla_c <- raster::as.data.frame(masked_taiga_sla_c, xy = TRUE)
taiga_sla_b <- raster::as.data.frame(masked_taiga_sla_b, xy = TRUE)
# 2) tundra forest boreal 
tundra_sla_c <- raster::as.data.frame(masked_tundra_sla_c, xy = TRUE)
tundra_sla_b <- raster::as.data.frame(masked_tundra_sla_b, xy = TRUE)
# 3) temperate coniferous forest 
tmp_conif_sla_c <- raster::as.data.frame(mask_temp_conif_sla_c, xy = TRUE)
tmp_conif_sla_b <- raster::as.data.frame(mask_temp_conif_sla_b, xy = TRUE)
# 4) temperate broad and mixed forest 
tmp_broad_mix_sla_c <- raster::as.data.frame(mask_temp_broad_mix_sla_c,
                                             xy = TRUE)
tmp_broad_mix_sla_b <- raster::as.data.frame(mask_temp_broad_mix_sla_b,
                                             xy = TRUE)
# 5) tropical and subtropical dry broadleaf
trpsbtrp_dry_broad_sla_c <- raster::as.data.frame(mask_trp_sbtrp_dry_broad_sla_c,
                                                   xy = TRUE)
trpsbtrp_dry_broad_sla_b <- raster::as.data.frame(mask_trp_sbtrp_dry_broad_sla_b,
                                                  xy = TRUE)
# 6) tropical and subtropical conif forest
trpsbtrp_conif_sla_c <- raster::as.data.frame(mask_trp_sbtrp_conif_sla_c,
                                              xy = TRUE)
trpsbtrp_conif_sla_b <- raster::as.data.frame(mask_trp_sbtrp_conif_sla_b,
                                              xy = TRUE)
# 7) tropical subtropical moist broadleaf biome
trpsbtrp_moist_broad_sla_c <- raster::as.data.frame(mask_trp_sbtrp_moist_broad_sla_c,
                                                    xy = TRUE)
trpsbtrp_moist_broad_sla_b <- raster::as.data.frame(mask_trp_sbtrp_moist_broad_sla_b,
                                                    xy = TRUE)
# 8) mediterranean forests, woodlands, scrub biome
med_f_w_scr_sla_c <- raster::as.data.frame(mask_med_f_w_scr_sla_c, 
                                           xy = TRUE)
med_f_w_scr_sla_b <- raster::as.data.frame(mask_med_f_w_scr_sla_b,
                                           xy = TRUE)
# 9) desertic and xeric scrubland
des_x_scr_sla_c <- raster::as.data.frame(mask_des_x_scr_sla_c, xy = TRUE)
des_x_scr_sla_b <- raster::as.data.frame(mask_des_x_scr_sla_b, xy = TRUE)
# 10) temperate grassland, savanna, shrubland
temp_grass_sav_shr_sla_c <- raster::as.data.frame(mask_temp_grass_sav_shr_sla_c,
                                                  xy = TRUE)
temp_grass_sav_shr_sla_b <- raster::as.data.frame(mask_temp_grass_sav_shr_sla_b,
                                                  xy = TRUE)
# 11) montane grassland and shrubland
mont_grass_shr_sla_c <- raster::as.data.frame(mask_mont_grass_shr_sla_c, 
                                              xy = TRUE)
mont_grass_shr_sla_b <- raster::as.data.frame(mask_mont_grass_shr_sla_b,
                                              xy = TRUE)
# 12) mangroves
mangr_sla_c <- raster::as.data.frame(mask_mangroves_sla_c, xy = TRUE)
mangr_sla_b <- raster::as.data.frame(mask_mangroves_sla_b, xy = TRUE)
# 13) flooded grassland and savanna
flo_grass_sav_sla_c <- raster::as.data.frame(mask_flo_grass_sav_sla_c,
                                             xy = TRUE)
flo_grass_sav_sla_b <- raster::as.data.frame(mask_flo_grass_sav_sla_b, 
                                             xy = TRUE)
# 14) tropical and subtropical grassland, savanna, shrubland
trpsbtrp_grass_sav_shr_sla_c <- raster::as.data.frame(mask_trp_sbtrp_grass_sav_shr_sla_c,
                                                       xy = TRUE)
trpsbtrp_grass_sav_shr_sla_b <- raster::as.data.frame(mask_trp_sbtrp_grass_sav_shr_sla_b,
                                                      xy = TRUE)

# sla stdev
# 1) taiga
taiga_slastd_c <- raster::as.data.frame(masked_taiga_slastd_c, xy = TRUE)
taiga_slastd_b <- raster::as.data.frame(masked_taiga_slastd_b, xy = TRUE)
# 2) tundra forest boreal 
tundra_slastd_c <- raster::as.data.frame(masked_tundra_slastd_c, xy = TRUE)
tundra_slastd_b <- raster::as.data.frame(masked_tundra_slastd_b, xy = TRUE)
# 3) temperate coniferous forest 
tmp_conif_slastd_c <- raster::as.data.frame(mask_temp_conif_slastd_c, 
                                            xy = TRUE)
tmp_conif_slastd_b <- raster::as.data.frame(mask_temp_conif_slastd_b,
                                            xy = TRUE)
# 4) temperate broad and mixed forest 
tmp_broad_mix_slastd_c <- raster::as.data.frame(mask_temp_broad_mix_slastd_c,
                                             xy = TRUE)
tmp_broad_mix_slastd_b <- raster::as.data.frame(mask_temp_broad_mix_slastd_b,
                                             xy = TRUE)
# 5) tropical and subtropical dry broadleaf
trpsbtrp_dry_broad_slastd_c <- raster::as.data.frame(mask_trp_sbtrp_dry_broad_slastd_c,
                                                  xy = TRUE)
trpsbtrp_dry_broad_slastd_b <- raster::as.data.frame(mask_trp_sbtrp_dry_broad_slastd_b,
                                                  xy = TRUE)
# 6) tropical and subtropical conif forest
trpsbtrp_conif_slastd_c <- raster::as.data.frame(mask_trp_sbtrp_conif_slastd_c,
                                              xy = TRUE)
trpsbtrp_conif_slastd_b <- raster::as.data.frame(mask_trp_sbtrp_conif_slastd_b,
                                              xy = TRUE)
# 7) tropical subtropical moist broadleaf biome
trpsbtrp_moist_broad_slastd_c <- raster::as.data.frame(mask_trp_sbtrp_moist_broad_slastd_c,
                                                    xy = TRUE)
trpsbtrp_moist_broad_slastd_b <- raster::as.data.frame(mask_trp_sbtrp_moist_broad_slastd_b,
                                                    xy = TRUE)
# 8) mediterranean forests, woodlands, scrub biome
med_f_w_scr_slastd_c <- raster::as.data.frame(mask_med_f_w_scr_slastd_c, 
                                           xy = TRUE)
med_f_w_scr_slastd_b <- raster::as.data.frame(mask_med_f_w_scr_slastd_b,
                                           xy = TRUE)
# 9) desertic and xeric scrubland
des_x_scr_slastd_c <- raster::as.data.frame(mask_des_x_scr_slastd_c, xy = TRUE)
des_x_scr_slastd_b <- raster::as.data.frame(mask_des_x_scr_slastd_b, xy = TRUE)
# 10) temperate grassland, savanna, shrubland
temp_grass_sav_shr_slastd_c <- raster::as.data.frame(mask_temp_grass_sav_shr_slastd_c,
                                                  xy = TRUE)
temp_grass_sav_shr_slastd_b <- raster::as.data.frame(mask_temp_grass_sav_shr_slastd_b,
                                                  xy = TRUE)
# 11) montane grassland and shrubland
mont_grass_shr_slastd_c <- raster::as.data.frame(mask_mont_grass_shr_slastd_c, 
                                              xy = TRUE)
mont_grass_shr_slastd_b <- raster::as.data.frame(mask_mont_grass_shr_slastd_b,
                                              xy = TRUE)
# 12) mangroves
mangr_slastd_c <- raster::as.data.frame(mask_mangroves_slastd_c, xy = TRUE)
mangr_slastd_b <- raster::as.data.frame(mask_mangroves_slastd_b, xy = TRUE)
# 13) flooded grassland and savanna
flo_grass_sav_slastd_c <- raster::as.data.frame(mask_flo_grass_sav_slastd_c,
                                             xy = TRUE)
flo_grass_sav_slastd_b <- raster::as.data.frame(mask_flo_grass_sav_slastd_b, 
                                             xy = TRUE)
# 14) tropical and subtropical grassland, savanna, shrubland
trpsbtrp_grass_sav_shr_slastd_c <- raster::as.data.frame(mask_trp_sbtrp_grass_sav_shr_slastd_c,
                                                      xy = TRUE)
trpsbtrp_grass_sav_shr_slastd_b <- raster::as.data.frame(mask_trp_sbtrp_grass_sav_shr_slastd_b,
                                                      xy = TRUE)


# joining dataframes cardamom + butler ----
  # sla mean
j_taiga_sla <- left_join(taiga_sla_c, taiga_sla_b) 
j_tundra_sla <- left_join(tundra_sla_c, tundra_sla_b) 
j_tmp_c_sla <- left_join(tmp_conif_sla_c, tmp_conif_sla_b) 
j_tmp_b_m_sla <- left_join(tmp_broad_mix_sla_c, tmp_broad_mix_sla_b) 
j_trp_sbtrp_d_b_sla <- left_join(trpsbtrp_dry_broad_sla_c, 
                                 trpsbtrp_dry_broad_sla_b) 
j_trp_sbtrp_c_sla <- left_join(trpsbtrp_conif_sla_c, 
                               trpsbtrp_conif_sla_b) 
j_trp_sbtrp_m_br_sla <- left_join(trpsbtrp_moist_broad_sla_c, 
                                  trpsbtrp_moist_broad_sla_b) 
j_med_f_sla <- left_join(med_f_w_scr_sla_c, med_f_w_scr_sla_b) 
j_des_x_s_sla <- left_join(des_x_scr_sla_c, des_x_scr_sla_b) 
j_temp_g_s_sh_sla <- left_join(temp_grass_sav_shr_sla_c, 
                               temp_grass_sav_shr_sla_b) 
j_mont_g_shr_sla <- left_join(mont_grass_shr_sla_c, 
                              mont_grass_shr_sla_b)
j_mangr_sla <- left_join(mangr_sla_c, mangr_sla_b)
j_flo_g_sav_sla <- left_join(flo_grass_sav_sla_c, flo_grass_sav_sla_b) 
  # joined flooded grass and savanna
j_trpsbtrp_g_sav_shr_sla <- left_join(trpsbtrp_grass_sav_shr_sla_c,
                                  trpsbtrp_grass_sav_shr_sla_b) 

## --- ##

  # sla stdev
j_taiga_slastd <- left_join(taiga_slastd_c, taiga_slastd_b)
j_tundra_slastd <- left_join(tundra_slastd_c, tundra_slastd_b)
j_tmp_c_slastd <- left_join(tmp_conif_slastd_c, tmp_conif_slastd_b)
j_tmp_b_m_slastd <- left_join(tmp_broad_mix_slastd_c,
                           tmp_broad_mix_slastd_b) 
j_trp_sbtrp_d_b_slastd <- left_join(trpsbtrp_dry_broad_slastd_c,
                                    trpsbtrp_dry_broad_slastd_b)
j_trp_sbtrp_c_slastd <- left_join(trpsbtrp_conif_slastd_c,
                                  trpsbtrp_conif_slastd_b)
j_trp_sbtrp_m_br_slastd <- left_join(trpsbtrp_moist_broad_slastd_c,
                                     trpsbtrp_moist_broad_slastd_b)
j_med_f_slastd <- left_join(med_f_w_scr_slastd_c,
                            med_f_w_scr_slastd_b)
j_des_x_s_slastd <- left_join(des_x_scr_slastd_c, des_x_scr_slastd_b)
j_temp_g_s_sh_slastd <- left_join(temp_grass_sav_shr_slastd_c, 
                                  temp_grass_sav_shr_slastd_b)
j_mont_g_shr_slastd <- left_join(mont_grass_shr_slastd_c, 
                                 mont_grass_shr_slastd_b)
j_mangr_slastd <- left_join(mangr_slastd_c, mangr_slastd_b)
j_flo_g_sav_slastd <- left_join(flo_grass_sav_slastd_c,
                                flo_grass_sav_slastd_b)
j_trpsbtrp_g_sav_shr_slastd <- left_join(trpsbtrp_grass_sav_shr_slastd_c,
                                         trpsbtrp_grass_sav_shr_slastd_b)

+
  
# DIFFERENCE HISTOGRAMS ----
## major biomes of interest ##
# 1) taiga ----
      # sla mean----
j_taiga_sla_h <- j_taiga_sla %>%
  gather(key = "dataset",value="sla", -x,-y)
  # - histogram for cardamom and butler in taiga biome
(taiga_sla_hist <- ggplot(j_taiga_sla_h, aes(x=sla,group=dataset,
                                                  fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (Taiga)", 
                        labels=c("Cardamom", "Butler")))

taiga_sla_hist_data <- as.data.table(ggplot_build(taiga_sla_hist)$data[1])
taiga_sla_hist_data <- taiga_sla_hist_data %>%
  dplyr::select(count, xmin, xmax, group) 
# atm for some reason, error "unused groups (count,xmin,xmax,group)" dont know why
taiga_sla_hist_data_c <- taiga_sla_hist_data[group==1]
taiga_sla_hist_data_b <- taiga_sla_hist_data[group==2]
taiga_sla_hist_diff <- merge(taiga_sla_hist_data_c, taiga_sla_hist_data_b, 
                             by=c("xmin","xmax"), 
                          suffixes = c(".c",".b"), allow.cartesian=TRUE)
taiga_sla_hist_diff <- taiga_sla_hist_diff[,Difference:=count.c-count.b]
setnames(taiga_sla_hist_diff, old = c("count.c","count.b"),new = c("Cardamom",
                                                                "Butler"))
taiga_sla_hist_diff_melt <- melt(taiga_sla_hist_diff, id.vars = c("xmin","xmax"), 
                         measure.vars = c("Cardamom","Difference","Butler"))

(taiga_sla_hist_diff_plot <- ggplot(taiga_sla_hist_diff_melt, 
                                    aes(xmin=xmin, xmax=xmax,
                                        ymax=value,ymin=0,
                                        group=variable,
                                        fill=variable, 
                                        color=variable,
                                        alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
#   xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Taiga\n") +
    theme(legend.title = element_blank())+
    xlim(0, 65)+
    ylim(-150,400))  

      # sla stdev ----
j_taiga_slastd_h <- j_taiga_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
# - histogram for cardamom and butler in taiga biome
(taiga_slastd_hist <- ggplot(j_taiga_slastd_h, aes(x=sla_std,group=dataset,
                                             fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (Taiga)"))

taiga_slastd_hist_data <- as.data.table(ggplot_build(taiga_slastd_hist)$data[1])
taiga_slastd_hist_data <- taiga_slastd_hist_data %>%
  dplyr::select(count, xmin, xmax, group) 
# atm for some reason, error "unused groups (count,xmin,xmax,group)" dont know why
taiga_slastd_hist_data_c <- taiga_slastd_hist_data[group==2]
taiga_slastd_hist_data_b <- taiga_slastd_hist_data[group==1]
taiga_slastd_hist_diff <- merge(taiga_slastd_hist_data_c, 
                                taiga_slastd_hist_data_b, 
                             by=c("xmin","xmax"), 
                             suffixes = c(".c",".b"), allow.cartesian=TRUE)
taiga_slastd_hist_diff <- taiga_slastd_hist_diff[,Difference:=count.c-count.b]
setnames(taiga_slastd_hist_diff, old = c("count.c","count.b"),new = c("Cardamom",
                                                                   "Butler"))
taiga_slastd_hist_diff_melt <- melt(taiga_slastd_hist_diff, 
                                    id.vars = c("xmin","xmax"), 
                                 measure.vars = c("Cardamom","Difference",
                                                  "Butler"))

(taiga_slastd_hist_diff_plot <- ggplot(taiga_slastd_hist_diff_melt, 
                                    aes(xmin=xmin, xmax=xmax,
                                        ymax=value,ymin=0,
                                        group=variable,
                                        fill=variable, 
                                        color=variable,
                                        alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
  #  xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Taiga\n") +
    theme(legend.title = element_blank())+
    xlim(0,75)+
    ylim(-500,500))  

# 2) tundra ---- 
     # sla mean ----
j_tundra_sla_hist <- j_tundra_sla %>%
  gather(dataset, sla, -x,-y)
(tundra_sla_hist <- ggplot(j_tundra_sla_hist, aes(x=sla,group=dataset,
                                             fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (Taiga)", 
                        labels=c("Cardamom", "Butler")))
tundra_sla_hist_data <- as.data.table(ggplot_build(tundra_sla_hist)$data[1])
tundra_sla_hist_data <- tundra_sla_hist_data %>%
  dplyr::select(count, xmin, xmax, group) 
tundra_sla_hist_data_c <- tundra_sla_hist_data[group==1]
tundra_sla_hist_data_b <- tundra_sla_hist_data[group==2]
tundra_sla_hist_diff <- merge(tundra_sla_hist_data_c, tundra_sla_hist_data_b, 
                             by=c("xmin","xmax"), 
                             suffixes = c(".c",".b"), allow.cartesian=TRUE)
tundra_sla_hist_diff <- tundra_sla_hist_diff[,Difference:=count.c-count.b]
setnames(tundra_sla_hist_diff, old = c("count.c","count.b"),new = c("Cardamom",
                                                                   "Butler"))
tundra_sla_hist_diff_melt <- melt(tundra_sla_hist_diff, id.vars = c("xmin","xmax"), 
                                 measure.vars = c("Cardamom","Difference","Butler"))

(tundra_sla_hist_diff_plot <- ggplot(tundra_sla_hist_diff_melt, 
                                    aes(xmin=xmin, xmax=xmax,
                                        ymax=value,ymin=0,
                                        group=variable,
                                        fill=variable, 
                                        color=variable,
                                        alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
  #  xlab("\nSpecific Leaf Area (m2.kg-1)")+
  #  ylab("Count\n")+
    ggtitle("Tundra\n")+
    xlim(0, 65)+
    ylim(-150,400))  

     # sla stdev ----
tundra_slastd_h <- j_tundra_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(tundra_slastd_h_plot <- ggplot(tundra_slastd_h, aes(x=sla_std,group=dataset,
                                                   fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (Tundra)"))
tundra_slastd_h_data <- as.data.table(ggplot_build(tundra_slastd_h_plot)$data[1])
tundra_slastd_h_data <- tundra_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
tundra_slastd_h_data_c <- tundra_slastd_h_data[group==2]
tundra_slastd_h_data_b <- tundra_slastd_h_data[group==1]
tundra_slastd_h_diff <- merge(tundra_slastd_h_data_c, 
                              tundra_slastd_h_data_b, 
                                by=c("xmin","xmax"), 
                                suffixes = c(".c",".b"), allow.cartesian=TRUE)
tundra_slastd_h_diff <- tundra_slastd_h_diff[,Difference:=count.c-count.b]
setnames(tundra_slastd_h_diff, old = c("count.c","count.b"),new = c("Cardamom",
                                                                      "Butler"))
tundra_slastd_h_melt <- melt(tundra_slastd_h_diff, 
                                    id.vars = c("xmin","xmax"), 
                                    measure.vars = c("Cardamom","Difference",
                                                     "Butler"))

(tundra_slastd_h_diff <- ggplot(tundra_slastd_h_melt, 
                                       aes(xmin=xmin, xmax=xmax,
                                           ymax=value,ymin=0,
                                           group=variable,
                                           fill=variable, 
                                           color=variable,
                                           alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    #  xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    #  ylab("Count\n")+
    ggtitle("Tundra\n") +
    theme(legend.title = element_blank())+
    xlim(0,75)+
    ylim(-500,500)) 

# 3) temp conif forest biome ----
      # sla mean ---- 
temp_con_sla_hist <- j_tmp_c_sla %>%
  gather(dataset, sla, -x,-y)
(temp_con_sla_hist_plot <- ggplot(temp_con_sla_hist, aes(x=sla,group=dataset,
                                                  fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (Temperate Coniferous)"))
temp_con_sla_hist_data <- as.data.table(ggplot_build(temp_con_sla_hist_plot)$data[1])
temp_con_sla_hist_data <- temp_con_sla_hist_data %>%
  dplyr::select(count, xmin, xmax, group)
temp_con_sla_hist_data_c <- temp_con_sla_hist_data[group==1]
temp_con_sla_hist_data_b <- temp_con_sla_hist_data[group==2]
temp_con_sla_hist_merge <- merge(temp_con_sla_hist_data_c,
                                 temp_con_sla_hist_data_b, 
                              by=c("xmin","xmax"), 
                              suffixes = c(".c",".b"), allow.cartesian=TRUE)
temp_con_sla_hist_merge <- temp_con_sla_hist_merge[,Difference:=count.c-count.b]
setnames(temp_con_sla_hist_merge, old = c("count.c","count.b"),new = c("Cardamom",
                                                                    "Butler"))
temp_con_sla_hist_merge <- melt(temp_con_sla_hist_merge, 
                                id.vars = c("xmin","xmax"), 
                                  measure.vars = c("Cardamom",
                                                   "Difference","Butler"))

(temp_con_sla_hist_diff_plot <- ggplot(temp_con_sla_hist_merge, 
                                     aes(xmin=xmin, xmax=xmax,
                                         ymax=value,ymin=0,
                                         group=variable,
                                         fill=variable, 
                                         color=variable,
                                         alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
  #  xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Temperate coniferous forest\n") +
    theme(legend.title = element_blank())+
    xlim(0, 33)+
    ylim(-15,25))  

      # sla stdev ----
temp_c_slastd_h <- j_temp_g_s_sh_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(temp_c_slastd_h_plot <- ggplot(temp_c_slastd_h, aes(x=sla_std,group=dataset,
                                                     fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
temp_c_slastd_h_data <- as.data.table(ggplot_build(temp_c_slastd_h_plot)$data[1])
temp_c_slastd_h_data <- temp_c_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
temp_c_slastd_h_data_c <- temp_c_slastd_h_data[group==2]
temp_c_slastd_h_data_b <- temp_c_slastd_h_data[group==1]
temp_c_slastd_h_diff <- merge(temp_c_slastd_h_data_c, 
                              temp_c_slastd_h_data_b, 
                              by=c("xmin","xmax"), 
                              suffixes = c(".c",".b"), allow.cartesian=TRUE)
temp_c_slastd_h_diff <- temp_c_slastd_h_diff[,Difference:=count.c-count.b]
setnames(temp_c_slastd_h_diff, old = c("count.c","count.b"),new = c("Cardamom",
                                                                    "Butler"))
temp_c_slastd_h_melt <- melt(temp_c_slastd_h_diff, 
                             id.vars = c("xmin","xmax"), 
                             measure.vars = c("Cardamom","Difference",
                                              "Butler"))

(temp_c_slastd_h_diff <- ggplot(temp_c_slastd_h_melt, 
                                aes(xmin=xmin, xmax=xmax,
                                    ymax=value,ymin=0,
                                    group=variable,
                                    fill=variable, 
                                    color=variable,
                                    alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    #xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Temperate coniferous forest\n") +
    theme(legend.title = element_blank())+
    xlim(0,75)+
    ylim(-500,500)) 

# 4) temp broad mix forest biome ----
      # sla mean ----
tmp_b_m_sla_h <- j_tmp_b_m_sla %>%
  gather(dataset, sla, -x, -y)
(tmp_b_m_sla_h_plot <- ggplot(tmp_b_m_sla_h, aes(x=sla,
                                                     group=dataset,
                                                     fill=dataset))+
                            geom_histogram(bins = 100, alpha=0.4)+
                            theme_ipsum()+
                            scale_fill_discrete(name="SLA Mean (temp broad/mix)", 
                                                labels=c("Cardamom", "Butler")))
tmp_b_m_sla_h_data <- as.data.table(ggplot_build(tmp_b_m_sla_h_plot)$data[1])
tmp_b_m_sla_h_data <- dplyr::select(tmp_b_m_sla_h_data, 
                               count, xmin, xmax, ymin, ymax, group)
tmp_b_m_sla_h_data_c <- tmp_b_m_sla_h_data[group==1]
tmp_b_m_sla_h_data_b <- tmp_b_m_sla_h_data[group==2]
tmp_b_m_sla_h_merge <- merge(tmp_b_m_sla_h_data_c,
                             tmp_b_m_sla_h_data_b, 
                                by=c("xmin","xmax"), 
                                suffixes = c(".c",".b"), allow.cartesian=TRUE)
tmp_b_m_sla_h_merge <- tmp_b_m_sla_h_merge[,Difference:=count.c-count.b]
setnames(tmp_b_m_sla_h_merge, old = c("count.c","count.b"),new = c("Cardamom",
                                                                      "Butler"))
tmp_b_m_sla_h_merge <- melt(tmp_b_m_sla_h_merge, 
                                    id.vars = c("xmin","xmax"), 
                                    measure.vars = c("Cardamom","Difference",
                                                     "Butler"))

(tmp_b_m_sla_h_diff <- ggplot(tmp_b_m_sla_h_merge, 
                                       aes(xmin=xmin, xmax=xmax,
                                           ymax=value,ymin=0,
                                           group=variable,
                                           fill=variable, 
                                           color=variable,
                                           alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Temperate broadleaf/mixed forest\n") +
    theme(legend.title = element_blank())+
    xlim(0, 65)+
    ylim(-150,400))  
    
      # sla stdev ----
temp_broad_mix_slastd_h <- j_tmp_b_m_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(temp_broad_mix_slastd_h_plot <- ggplot(temp_broad_mix_slastd_h, 
                                        aes(x=sla_std,group=dataset,
                                                     fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
temp_broad_mix_slastd_h_data <- as.data.table(ggplot_build(temp_broad_mix_slastd_h_plot)$data[1])
temp_broad_mix_slastd_h_data <- temp_broad_mix_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
temp_broad_mix_slastd_h_data_c <- temp_broad_mix_slastd_h_data[group==2]
temp_broad_mix_slastd_h_data_b <- temp_broad_mix_slastd_h_data[group==1]
temp_broad_mix_slastd_h_merge <- merge(temp_broad_mix_slastd_h_data_c, 
                                       temp_broad_mix_slastd_h_data_b, 
                              by=c("xmin","xmax"), 
                              suffixes = c(".c",".b"), allow.cartesian=TRUE)
temp_broad_mix_slastd_h_merge <- temp_broad_mix_slastd_h_merge[,Difference:=count.c-count.b]
setnames(temp_broad_mix_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
temp_broad_mix_slastd_h_merge <- melt(temp_broad_mix_slastd_h_merge, 
                             id.vars = c("xmin","xmax"), 
                             measure.vars = c("Cardamom","Difference",
                                              "Butler"))
(temp_broad_mix_slastd_h_diff <- ggplot(temp_broad_mix_slastd_h_merge, 
                                aes(xmin=xmin, xmax=xmax,
                                    ymax=value,ymin=0,
                                    group=variable,
                                    fill=variable, 
                                    color=variable,
                                    alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Temperate broadleaf/mixed forest\n") +
    theme(legend.title = element_blank())+
    xlim(0,75)+
    ylim(-500,500)) 

# 5) tropical and subtropical dry broadleaf biome ----
      # sla mean ----
trpsbtrp_d_broad_h <- j_trp_sbtrp_d_b_sla %>%
  gather(dataset, sla, -x, -y)
(trpsbtrp_d_broad_h_plot <- ggplot(trpsbtrp_d_broad_h, 
                                   aes(sla, group = dataset, 
                                       fill = dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (tropical and subtropical dry broadleaf)",
                        labels=c("Cardamom", "Butler")))
trpsbtrp_d_broad_h_data <- as.data.table(ggplot_build(trpsbtrp_d_broad_h_plot)$data[1])
trpsbtrp_d_broad_h_data <- dplyr::select(trpsbtrp_d_broad_h_data,
                                         count, xmin, xmax, ymin, ymax, group)
trpsbtrp_d_broad_h_data_c <- trpsbtrp_d_broad_h_data[group == 1]
trpsbtrp_d_broad_h_data_b <- trpsbtrp_d_broad_h_data[group == 2]
trpsbtrp_d_broad_h_merge <- merge(trpsbtrp_d_broad_h_data_c,
                                  trpsbtrp_d_broad_h_data_b, 
                             by=c("xmin","xmax"), 
                             suffixes = c(".c",".b"), allow.cartesian=TRUE)
trpsbtrp_d_broad_h_merge <- trpsbtrp_d_broad_h_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_d_broad_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
trpsbtrp_d_broad_h_merge <- melt(trpsbtrp_d_broad_h_merge, 
                            id.vars = c("xmin","xmax"), 
                            measure.vars = c("Cardamom","Difference",
                                             "Butler"))
(trpsbtrp_d_broad_h_diff <- ggplot(trpsbtrp_d_broad_h_merge, 
                              aes(xmin=xmin, xmax=xmax,
                                  ymax=value,ymin=0,
                                  group=variable,
                                  fill=variable, 
                                  color=variable,
                                  alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
   # xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Tropical and subtropical dry broadleaf forest\n") +
    theme(legend.title = element_blank())+
    xlim(0, 33)+
    ylim(-15,25))  

      # sla stdev ----
trpsbtrp_dry_broad_slastd_h <- j_trp_sbtrp_d_b_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(trpsbtrp_dry_broad_slastd_h_plot <- ggplot(trpsbtrp_dry_broad_slastd_h, 
                                        aes(x=sla_std,group=dataset,
                                            fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
trpsbtrp_dry_broad_slastd_h_data <- as.data.table(ggplot_build(trpsbtrp_dry_broad_slastd_h_plot)$data[1])
trpsbtrp_dry_broad_slastd_h_data <- trpsbtrp_dry_broad_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
trpsbtrp_dry_broad_slastd_h_data_c <- trpsbtrp_dry_broad_slastd_h_data[group==2]
trpsbtrp_dry_broad_slastd_h_data_b <- trpsbtrp_dry_broad_slastd_h_data[group==1]
trpsbtrp_dry_broad_slastd_h_merge <- merge(trpsbtrp_dry_broad_slastd_h_data_c, 
                                           trpsbtrp_dry_broad_slastd_h_data_b, 
                                       by=c("xmin","xmax"), 
                                       suffixes = c(".c",".b"), allow.cartesian=TRUE)
trpsbtrp_dry_broad_slastd_h_merge <- trpsbtrp_dry_broad_slastd_h_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_dry_broad_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
trpsbtrp_dry_broad_slastd_h_merge <- melt(trpsbtrp_dry_broad_slastd_h_merge, 
                                      id.vars = c("xmin","xmax"), 
                                      measure.vars = c("Cardamom","Difference",
                                                       "Butler"))
(trpsbtrp_dry_broad_slastd_h_diff <- ggplot(trpsbtrp_dry_broad_slastd_h_merge, 
                                        aes(xmin=xmin, xmax=xmax,
                                            ymax=value,ymin=0,
                                            group=variable,
                                            fill=variable, 
                                            color=variable,
                                            alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    #  xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Tropical and subtropical dry broadleaf forest\n") +
    theme(legend.title = element_blank())+
    xlim(0,65)+
    ylim(-50,50)) 

# 6) tropical and subtropical conif forest biome ----
      # sla mean ----
trpsbtrp_conif_h <- j_trp_sbtrp_c_sla %>%
  gather(dataset, sla, -x, -y)
(trpsbtrp_conif_h_plot <- ggplot(trpsbtrp_conif_h, aes(sla, 
                                                       group = dataset,
                                                       fill = dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (tropical and subtropical coniferous biome)",
                        labels=c("Cardamom", "Butler")))
trpsbtrp_conif_h_data <- as.data.table(ggplot_build(trpsbtrp_conif_h_plot)$data[1])
trpsbtrp_conif_h_data <- dplyr::select(trpsbtrp_conif_h_data,
                                       count, xmin,xmax,ymin,ymax,group)
trpsbtrp_conif_h_data_c <- trpsbtrp_conif_h_data[group == 1]
trpsbtrp_conif_h_data_b <- trpsbtrp_conif_h_data[group == 2]
trpsbtrp_conif_h_data_merge <- merge(trpsbtrp_conif_h_data_c,
                                     trpsbtrp_conif_h_data_b, 
                                  by=c("xmin","xmax"), 
                                  suffixes = c(".c",".b"), 
                                  allow.cartesian=TRUE)
trpsbtrp_conif_h_data_merge <- trpsbtrp_conif_h_data_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_conif_h_data_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
trpsbtrp_conif_h_data_merge <- melt(trpsbtrp_conif_h_data_merge, 
                                 id.vars = c("xmin","xmax"), 
                                 measure.vars = c("Cardamom","Difference",
                                                  "Butler"))
(trpsbtrp_conif_h_diff <- ggplot(trpsbtrp_conif_h_data_merge, 
                                   aes(xmin=xmin, xmax=xmax,
                                       ymax=value,ymin=0,
                                       group=variable,
                                       fill=variable, 
                                       color=variable,
                                       alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Tropical and subtropical coniferous forest\n") +
    theme(legend.title = element_blank())+
    xlim(0, 33)+
    ylim(-15,25))  

      # sla stdev ----
trpsbtrp_conif_slastd_h <- j_trp_sbtrp_c_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(trpsbtrp_conif_slastd_h_plot <- ggplot(trpsbtrp_conif_slastd_h, 
                                            aes(x=sla_std,group=dataset,
                                                fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
trpsbtrp_conif_slastd_h_data <- as.data.table(ggplot_build(trpsbtrp_conif_slastd_h_plot)$data[1])
trpsbtrp_conif_slastd_h_data <- trpsbtrp_conif_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
trpsbtrp_conif_slastd_h_data_c <- trpsbtrp_conif_slastd_h_data[group==2]
trpsbtrp_conif_slastd_h_data_b <- trpsbtrp_conif_slastd_h_data[group==1]
trpsbtrp_conif_slastd_h_merge <- merge(trpsbtrp_conif_slastd_h_data_c, 
                                       trpsbtrp_conif_slastd_h_data_b, 
                                           by=c("xmin","xmax"), 
                                           suffixes = c(".c",".b"), allow.cartesian=TRUE)
trpsbtrp_conif_slastd_h_merge <- trpsbtrp_conif_slastd_h_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_conif_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
trpsbtrp_conif_slastd_h_merge <- melt(trpsbtrp_conif_slastd_h_merge, 
                                          id.vars = c("xmin","xmax"), 
                                          measure.vars = c("Cardamom","Difference",
                                                           "Butler"))
(trpsbtrp_conif_slastd_h_diff <- ggplot(trpsbtrp_conif_slastd_h_merge, 
                                            aes(xmin=xmin, xmax=xmax,
                                                ymax=value,ymin=0,
                                                group=variable,
                                                fill=variable, 
                                                color=variable,
                                                alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    #  xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    #  ylab("Count\n")+
    ggtitle("Tropical and subtropical coniferous forest\n") +
    theme(legend.title = element_blank())+
    xlim(0,65)+
    ylim(-50,50)) 
# 7) tropical subtropical moist broadleaf biome ----
      # sla mean ----
trpsbtrp_m_broad_h <- j_trp_sbtrp_m_br_sla %>%
  gather(dataset, sla, -x, -y )
(trpsbtrp_m_broad_h_plot <- ggplot(trpsbtrp_m_broad_h, aes(sla, 
                                                           group = dataset,
                                                           fill = dataset))+
    geom_histogram(bins = 50, alpha=0.4)+
    theme_ipsum()+
    scale_fill_discrete(name="SLA Mean (tropical and subtropical moist broadleaf forest biome)",
                        labels=c("Cardamom", "Butler")))
trpsbtrp_m_broad_h_data <- as.data.table(ggplot_build(trpsbtrp_m_broad_h_plot)$data[1])
trpsbtrp_m_broad_h_data_c <- trpsbtrp_m_broad_h_data[group==1]
trpsbtrp_m_broad_h_data_b <- trpsbtrp_m_broad_h_data[group==2]
trpsbtrp_m_broad_h_merge <- merge(trpsbtrp_m_broad_h_data_c,
                                  trpsbtrp_m_broad_h_data_b, 
                                  by=c("xmin","xmax"), 
                                  suffixes = c(".c",".b"), 
                                  allow.cartesian=TRUE)
trpsbtrp_m_broad_h_merge <- trpsbtrp_m_broad_h_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_m_broad_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
trpsbtrp_m_broad_h_merge <- melt(trpsbtrp_m_broad_h_merge, 
                                    id.vars = c("xmin","xmax"), 
                                    measure.vars = c("Cardamom","Difference",
                                                     "Butler"))
(trpsbtrp_m_broad_h_diff <- ggplot(trpsbtrp_m_broad_h_merge, 
                                 aes(xmin=xmin, xmax=xmax,
                                     ymax=value,ymin=0,
                                     group=variable,
                                     fill=variable, 
                                     color=variable,
                                     alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Tropical and subtropical moist broadleaf forest\n") +
    theme(legend.title = element_blank())+
    xlim(0, 65)+
    ylim(-150,400))  

      # sla stdev ----
trpsbtrp_m_broad_slastd_h <- j_trp_sbtrp_m_br_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(trpsbtrp_m_broad_slastd_h_plot <- ggplot(trpsbtrp_m_broad_slastd_h, 
                                        aes(x=sla_std,group=dataset,
                                            fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
trpsbtrp_m_broad_slastd_h_data <- as.data.table(ggplot_build(trpsbtrp_m_broad_slastd_h_plot)$data[1])
trpsbtrp_m_broad_slastd_h_data <- trpsbtrp_m_broad_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
trpsbtrp_m_broad_slastd_h_data_c <- trpsbtrp_m_broad_slastd_h_data[group==2]
trpsbtrp_m_broad_slastd_h_data_b <- trpsbtrp_m_broad_slastd_h_data[group==1]
trpsbtrp_m_broad_slastd_h_merge <- merge(trpsbtrp_m_broad_slastd_h_data_c, 
                                         trpsbtrp_m_broad_slastd_h_data_b, 
                                       by=c("xmin","xmax"), 
                                       suffixes = c(".c",".b"), allow.cartesian=TRUE)
trpsbtrp_m_broad_slastd_h_merge <- trpsbtrp_m_broad_slastd_h_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_m_broad_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
trpsbtrp_m_broad_slastd_h_merge <- melt(trpsbtrp_m_broad_slastd_h_merge, 
                                      id.vars = c("xmin","xmax"), 
                                      measure.vars = c("Cardamom","Difference",
                                                       "Butler"))
(trpsbtrp_m_broad_slastd_h_diff <- ggplot(trpsbtrp_m_broad_slastd_h_merge, 
                                        aes(xmin=xmin, xmax=xmax,
                                            ymax=value,ymin=0,
                                            group=variable,
                                            fill=variable, 
                                            color=variable,
                                            alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
   # ylab("Count\n")+
    ggtitle("Tropical and subtropical moist broadleaf forest\n") +
    theme(legend.title = element_blank())+
    xlim(0,75)+
    ylim(-500,500)) 

# 8) mediterranean forests, woodlands, scrub ----
      # sla mean ----
med_f_w_scr_h <- j_med_f_sla %>%
  gather(dataset, sla, -x, -y)
(med_f_w_scr_h_plot <- ggplot(med_f_w_scr_h, aes(sla, group = dataset, 
                                                 fill = dataset))+
    geom_histogram(bins = 100, alpha = 0.4)+
    theme_ipsum())
med_f_w_scr_h_plot_data <- as.data.table(ggplot_build(med_f_w_scr_h_plot)$data[1])
med_f_w_scr_h_plot_data <- dplyr::select(med_f_w_scr_h_plot_data, 
                                         count, xmin, xmax, ymin, ymax, group)
med_f_w_scr_h_plot_data_c <- med_f_w_scr_h_plot_data[group == 1]
med_f_w_scr_h_plot_data_b <- med_f_w_scr_h_plot_data[group == 2]
med_f_w_scr_h_plot_merge <- merge(med_f_w_scr_h_plot_data_c,
                                  med_f_w_scr_h_plot_data_b, 
                                  by=c("xmin","xmax"), 
                                  suffixes = c(".c",".b"), 
                                  allow.cartesian=TRUE)
med_f_w_scr_h_plot_merge <- med_f_w_scr_h_plot_merge[,Difference:=count.c-count.b]
setnames(med_f_w_scr_h_plot_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
med_f_w_scr_h_plot_merge <- melt(med_f_w_scr_h_plot_merge, 
                                 id.vars = c("xmin","xmax"), 
                                 measure.vars = c("Cardamom","Difference",
                                                  "Butler"))
(med_f_w_scr_h_plot_diff <- ggplot(med_f_w_scr_h_plot_merge, 
                                   aes(xmin=xmin, xmax=xmax,
                                       ymax=value,ymin=0,
                                       group=variable,
                                       fill=variable, 
                                       color=variable,
                                       alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area (m2.kg-1)")+
 #  ylab("Count\n")+
    ggtitle("Mediterranean forest, woodland and scrubland\n") +
    theme(legend.title = element_blank())+
    xlim(0,33)+
    ylim(-15,25)) 

      # sla stdev ----
med_f_w_scr_slastd_h <- j_med_f_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(med_f_w_scr_slastd_h_plot <- ggplot(med_f_w_scr_slastd_h, 
                                          aes(x=sla_std,group=dataset,
                                              fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
med_f_w_scr_slastd_h_data <- as.data.table(ggplot_build(med_f_w_scr_slastd_h_plot)$data[1])
med_f_w_scr_slastd_h_data <- med_f_w_scr_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
med_f_w_scr_slastd_h_data_c <- med_f_w_scr_slastd_h_data[group==2]
med_f_w_scr_slastd_h_data_b <- med_f_w_scr_slastd_h_data[group==1]
med_f_w_scr_slastd_h_merge <- merge(med_f_w_scr_slastd_h_data_c, 
                                    med_f_w_scr_slastd_h_data_b, 
                                         by=c("xmin","xmax"), 
                                         suffixes = c(".c",".b"), allow.cartesian=TRUE)
med_f_w_scr_slastd_h_merge <- med_f_w_scr_slastd_h_merge[,Difference:=count.c-count.b]
setnames(med_f_w_scr_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
med_f_w_scr_slastd_h_merge <- melt(med_f_w_scr_slastd_h_merge, 
                                        id.vars = c("xmin","xmax"), 
                                        measure.vars = c("Cardamom","Difference",
                                                         "Butler"))
(med_f_w_scr_slastd_h_diff <- ggplot(med_f_w_scr_slastd_h_merge, 
                                          aes(xmin=xmin, xmax=xmax,
                                              ymax=value,ymin=0,
                                              group=variable,
                                              fill=variable, 
                                              color=variable,
                                              alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    #ylab("Count\n")+
    ggtitle("Mediterranean forests, woodlands and scrubland\n") +
    theme(legend.title = element_blank())+
    xlim(0,65)+
    ylim(-50,50)) 

## other biomes ##
# 9) desertic and xeric scrubland ---- 
      # sla mean ---- 
des_x_scr_h <- j_des_x_s_sla %>%
  gather(dataset, sla, -x, -y)
(des_x_scr_h_plot <- ggplot(des_x_scr_h, aes(sla, group = dataset, 
                                             fill = dataset))+
    geom_histogram(bins = 100, alpha = 0.4)+
    theme_ipsum())
des_x_scr_h_plot_data <- as.data.table(ggplot_build(des_x_scr_h_plot)$data[1])
des_x_scr_h_plot_data <- dplyr::select(des_x_scr_h_plot_data, 
                                       count,xmin,xmax,ymin,ymax,group)
des_x_scr_h_plot_data_c <- des_x_scr_h_plot_data[group == 1]
des_x_scr_h_plot_data_b <- des_x_scr_h_plot_data[group == 2]
des_x_scr_h_plot_merge <- merge(des_x_scr_h_plot_data_c, 
                                des_x_scr_h_plot_data_b, 
                                by=c("xmin","xmax"), 
                                suffixes = c(".c",".b"), 
                                allow.cartesian=TRUE)
des_x_scr_h_plot_merge <- des_x_scr_h_plot_merge[,Difference:=count.c-count.b]
setnames(des_x_scr_h_plot_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
des_x_scr_h_plot_merge <- melt(des_x_scr_h_plot_merge, 
                                 id.vars = c("xmin","xmax"), 
                                 measure.vars = c("Cardamom","Difference",
                                                  "Butler"))
(des_x_scr_h_plot_merge_diff <- ggplot(des_x_scr_h_plot_merge, 
                                   aes(xmin=xmin, xmax=xmax,
                                       ymax=value,ymin=0,
                                       group=variable,
                                       fill=variable, 
                                       color=variable,
                                       alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
 #   xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Desert and xeric scrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-80,150)+
    xlim(0,40))

      # sla stdev ----
des_x_scr_slastd_h <- j_des_x_s_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(des_x_scr_slastd_h_plot <- ggplot(des_x_scr_slastd_h, 
                                     aes(x=sla_std,group=dataset,
                                         fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
des_x_scr_slastd_h_data <- as.data.table(ggplot_build(des_x_scr_slastd_h_plot)$data[1])
des_x_scr_slastd_h_data <- des_x_scr_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
des_x_scr_slastd_h_data_c <- des_x_scr_slastd_h_data[group==2]
des_x_scr_slastd_h_data_b <- des_x_scr_slastd_h_data[group==1]
des_x_scr_slastd_h_merge <- merge(des_x_scr_slastd_h_data_c, 
                                  des_x_scr_slastd_h_data_b, 
                                    by=c("xmin","xmax"), 
                                    suffixes = c(".c",".b"), allow.cartesian=TRUE)
des_x_scr_slastd_h_merge <- des_x_scr_slastd_h_merge[,Difference:=count.c-count.b]
setnames(des_x_scr_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
des_x_scr_slastd_h_merge <- melt(des_x_scr_slastd_h_merge, 
                                   id.vars = c("xmin","xmax"), 
                                   measure.vars = c("Cardamom","Difference",
                                                    "Butler"))
(des_x_scr_slastd_h_diff <- ggplot(des_x_scr_slastd_h_merge, 
                                     aes(xmin=xmin, xmax=xmax,
                                         ymax=value,ymin=0,
                                         group=variable,
                                         fill=variable, 
                                         color=variable,
                                         alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    #  xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Desert and xeric scrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-200, 200)+
    xlim(0, 70)) 

# 10) temperate grassland, savanna, shrubland ----
      # sla mean ----
temp_grass_sav_shr_h <- j_temp_g_s_sh_sla %>%
  gather(dataset, sla, -x, -y)
(temp_grass_sav_shr_h_plot <- ggplot(temp_grass_sav_shr_h, 
                                     aes(sla, group = dataset, 
                                         fill = dataset))+
    geom_histogram(bins = 100, alpha = 0.4)+
    theme_ipsum())
temp_grass_sav_shr_h_data <- as.data.table(ggplot_build(temp_grass_sav_shr_h_plot)$data[1])
temp_grass_sav_shr_h_data <- dplyr::select(temp_grass_sav_shr_h_data, 
                                           count,xmin,xmax,ymin,ymax,group)
temp_grass_sav_shr_h_data_c <- temp_grass_sav_shr_h_data[group == 1]
temp_grass_sav_shr_h_data_b <- temp_grass_sav_shr_h_data[group == 2]
temp_grass_sav_shr_h_merge <- merge(temp_grass_sav_shr_h_data_c,
                                    temp_grass_sav_shr_h_data_b, 
                                    by=c("xmin","xmax"), 
                                    suffixes = c(".c",".b"), 
                                    allow.cartesian=TRUE)
temp_grass_sav_shr_h_merge <- temp_grass_sav_shr_h_merge[,Difference:=count.c-count.b]
setnames(temp_grass_sav_shr_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
temp_grass_sav_shr_h_merge <- melt(temp_grass_sav_shr_h_merge, 
                               id.vars = c("xmin","xmax"), 
                               measure.vars = c("Cardamom","Difference",
                                                "Butler"))
(temp_grass_sav_shr_h_diff <- ggplot(temp_grass_sav_shr_h_merge, 
                                       aes(xmin=xmin, xmax=xmax,
                                           ymax=value,ymin=0,
                                           group=variable,
                                           fill=variable, 
                                           color=variable,
                                           alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
  #  xlab("\nSpecific Leaf Area (m2.kg-1)")+
   # ylab("Count\n")+
    ggtitle("Temperate grassland, savanna and shrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-80,150)+
    xlim(0,40))

      # sla stdev ----
temp_grass_sav_shr_slastd_h <- j_temp_g_s_sh_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(temp_grass_sav_shr_slastd_h_plot <- ggplot(temp_grass_sav_shr_slastd_h, 
                                   aes(x=sla_std,group=dataset,
                                       fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
temp_grass_sav_shr_slastd_h_data <- as.data.table(ggplot_build(temp_grass_sav_shr_slastd_h_plot)$data[1])
temp_grass_sav_shr_slastd_h_data <- temp_grass_sav_shr_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
temp_grass_sav_shr_slastd_h_data_c <- temp_grass_sav_shr_slastd_h_data[group==2]
temp_grass_sav_shr_slastd_h_data_b <- temp_grass_sav_shr_slastd_h_data[group==1]
temp_grass_sav_shr_slastd_h_merge <- merge(temp_grass_sav_shr_slastd_h_data_c, 
                                           temp_grass_sav_shr_slastd_h_data_b, 
                                  by=c("xmin","xmax"), 
                                  suffixes = c(".c",".b"), allow.cartesian=TRUE)
temp_grass_sav_shr_slastd_h_merge <- temp_grass_sav_shr_slastd_h_merge[,Difference:=count.c-count.b]
setnames(temp_grass_sav_shr_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
temp_grass_sav_shr_slastd_h_merge <- melt(temp_grass_sav_shr_slastd_h_merge, 
                                 id.vars = c("xmin","xmax"), 
                                 measure.vars = c("Cardamom","Difference",
                                                  "Butler"))
(temp_grass_sav_shr_slastd_h_diff <- ggplot(temp_grass_sav_shr_slastd_h_merge, 
                                   aes(xmin=xmin, xmax=xmax,
                                       ymax=value,ymin=0,
                                       group=variable,
                                       fill=variable, 
                                       color=variable,
                                       alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    #  ylab("Count\n")+
    ggtitle("Temperate grassland, savanna and shrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-200, 200)+
    xlim(0, 70)) 

# 11) montane grassland and shrubland biome ----
      # sla mean ----
mont_grass_shr_h <- j_mont_g_shr_sla %>%
  gather(dataset, sla, -x, -y)
mont_grass_shr_h_plot <- ggplot(mont_grass_shr_h, 
                                 aes(sla, group = dataset,
                                     fill = dataset))+
  geom_histogram(bins = 100, alpha = 0.4)
mont_grass_shr_h_plot_data <- as.data.table(ggplot_build(mont_grass_shr_h_plot)$data[1])
mont_grass_shr_h_plot_data <- dplyr::select(mont_grass_shr_h_plot_data, 
                                            count,xmin,xmax,ymin,ymax,group)
mont_grass_shr_h_plot_data_c <- mont_grass_shr_h_plot_data[group == 1]  
mont_grass_shr_h_plot_data_b <- mont_grass_shr_h_plot_data[group == 2]
mont_grass_shr_h_plot_merge <- merge(mont_grass_shr_h_plot_data_c,
                                     mont_grass_shr_h_plot_data_b, 
                                     by=c("xmin","xmax"), 
                                     suffixes = c(".c",".b"), 
                                     allow.cartesian=TRUE)
mont_grass_shr_h_plot_merge <- mont_grass_shr_h_plot_merge[,Difference:=count.c-count.b]
setnames(mont_grass_shr_h_plot_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
mont_grass_shr_h_plot_merge <- melt(mont_grass_shr_h_plot_merge, 
                                   id.vars = c("xmin","xmax"), 
                                   measure.vars = c("Cardamom","Difference",
                                                    "Butler"))
(mont_grass_shr_h_diff <- ggplot(mont_grass_shr_h_plot_merge, 
                                     aes(xmin=xmin, xmax=xmax,
                                         ymax=value,ymin=0,
                                         group=variable,
                                         fill=variable, 
                                         color=variable,
                                         alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
#    xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Montane grassland and shrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-25,50)+
    xlim(0,25))

      # sla stdev ----
mont_grass_shr_slastd_h <- j_mont_g_shr_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(mont_grass_shr_slastd_h_plot <- ggplot(mont_grass_shr_slastd_h, 
                                            aes(x=sla_std,group=dataset,
                                                fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
mont_grass_shr_slastd_h_data <- as.data.table(ggplot_build(mont_grass_shr_slastd_h_plot)$data[1])
mont_grass_shr_slastd_h_data <- mont_grass_shr_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
mont_grass_shr_slastd_h_data_c <- mont_grass_shr_slastd_h_data[group==2]
mont_grass_shr_slastd_h_data_b <- mont_grass_shr_slastd_h_data[group==1]
mont_grass_shr_slastd_h_merge <- merge(mont_grass_shr_slastd_h_data_c, 
                                       mont_grass_shr_slastd_h_data_b, 
                                           by=c("xmin","xmax"), 
                                           suffixes = c(".c",".b"), allow.cartesian=TRUE)
mont_grass_shr_slastd_h_merge <- mont_grass_shr_slastd_h_merge[,Difference:=count.c-count.b]
setnames(mont_grass_shr_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
mont_grass_shr_slastd_h_merge <- melt(mont_grass_shr_slastd_h_merge, 
                                          id.vars = c("xmin","xmax"), 
                                          measure.vars = c("Cardamom","Difference",
                                                           "Butler"))
(mont_grass_shr_slastd_h_diff <- ggplot(mont_grass_shr_slastd_h_merge, 
                                            aes(xmin=xmin, xmax=xmax,
                                                ymax=value,ymin=0,
                                                group=variable,
                                                fill=variable, 
                                                color=variable,
                                                alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    #  xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Montane grassland and scrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-70, 70)+
    xlim(0, 65)) 

# 12) mangrove ----
      # sla mean ----
mangr_h <- j_mangr_sla %>%
  gather(dataset, sla, -x, -y)
(mangr_h_plot <- ggplot(mangr_h, aes(sla, group = dataset, 
                                     fill = dataset))+
    geom_histogram(bins = 100, alpha = 0.4)+
    theme_ipsum())
mangr_h_plot_data <- as.data.table(ggplot_build(mangr_h_plot)$data[1])
mangr_h_plot_data <- dplyr::select(mangr_h_plot_data, 
                                   count, xmin, xmax, ymin, ymax, group)
mangr_h_plot_data_c <- mangr_h_plot_data[group == 1]
mangr_h_plot_data_b <- mangr_h_plot_data[group == 2]
mangr_h_plot_merge <- merge(mangr_h_plot_data_c, 
                            mangr_h_plot_data_b, 
                            by=c("xmin","xmax"), 
                            suffixes = c(".c",".b"), 
                            allow.cartesian=TRUE)
mangr_h_plot_merge <- mangr_h_plot_merge[,Difference:=count.c-count.b]
setnames(mangr_h_plot_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
mangr_h_plot_merge <- melt(mangr_h_plot_merge, 
                                   id.vars = c("xmin","xmax"), 
                                   measure.vars = c("Cardamom","Difference",
                                                    "Butler"))
(mangr_h_plot_diff <- ggplot(mangr_h_plot_merge, 
                                     aes(xmin=xmin, xmax=xmax,
                                         ymax=value,ymin=0,
                                         group=variable,
                                         fill=variable, 
                                         color=variable,
                                         alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
  #  xlab("\nSpecific Leaf Area (m2.kg-1)")+
  #  ylab("Count\n")+
    ggtitle("Mangrove\n") +
    theme(legend.title = element_blank())+
    ylim(-25,50)+
    xlim(0,25))

      # sla stdev ----
mangr_slastd_h <- j_mangr_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(mangr_slastd_h_plot <- ggplot(mangr_slastd_h, 
                                        aes(x=sla_std,group=dataset,
                                            fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
mangr_slastd_h_data <- as.data.table(ggplot_build(mangr_slastd_h_plot)$data[1])
mangr_slastd_h_data <- mangr_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
mangr_slastd_h_data_c <- mangr_slastd_h_data[group==2]
mangr_slastd_h_data_b <- mangr_slastd_h_data[group==1]
mangr_slastd_h_merge <- merge(mangr_slastd_h_data_c, 
                              mangr_slastd_h_data_b, 
                                       by=c("xmin","xmax"), 
                                       suffixes = c(".c",".b"), allow.cartesian=TRUE)
mangr_slastd_h_merge <- mangr_slastd_h_merge[,Difference:=count.c-count.b]
setnames(mangr_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
mangr_slastd_h_merge <- melt(mangr_slastd_h_merge, 
                                      id.vars = c("xmin","xmax"), 
                                      measure.vars = c("Cardamom","Difference",
                                                       "Butler"))
(mangr_slastd_h_diff <- ggplot(mangr_slastd_h_merge, 
                                        aes(xmin=xmin, xmax=xmax,
                                            ymax=value,ymin=0,
                                            group=variable,
                                            fill=variable, 
                                            color=variable,
                                            alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    #  xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    #  ylab("Count\n")+
    ggtitle("Mangrove\n") +
    theme(legend.title = element_blank())+
    ylim(-70, 70)+
    xlim(0, 65)) 

# 13) flooded grassland and savanna ----
      # sla mean ----
flo_grass_sav_h <- j_flo_g_sav_sla %>%
  gather(dataset, sla, -x, -y)
(flo_grass_sav_h_plot <- ggplot(flo_grass_sav_h, aes(sla, 
                                                     group = dataset,
                                                     fill = dataset))+
    geom_histogram(bins = 100, alpha = 0.4)+
    theme_ipsum())
flo_grass_sav_h_plot_data <- as.data.table(ggplot_build(flo_grass_sav_h_plot)$data[1])
flo_grass_sav_h_plot_data <- dplyr::select(flo_grass_sav_h_plot_data, 
                                           count, xmin, xmax, ymin, ymax, group)
flo_grass_sav_h_plot_data_c <- flo_grass_sav_h_plot_data[group == 1]
flo_grass_sav_h_plot_data_b <- flo_grass_sav_h_plot_data[group == 2]
flo_grass_sav_h_plot_merge <- merge(flo_grass_sav_h_plot_data_c,
                                    flo_grass_sav_h_plot_data_b, 
                                    by=c("xmin","xmax"), 
                                    suffixes = c(".c",".b"), 
                                    allow.cartesian=TRUE)
flo_grass_sav_h_plot_merge <- flo_grass_sav_h_plot_merge[,Difference:=count.c-count.b]
setnames(flo_grass_sav_h_plot_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
flo_grass_sav_h_plot_merge <- melt(flo_grass_sav_h_plot_merge, 
                           id.vars = c("xmin","xmax"), 
                           measure.vars = c("Cardamom","Difference",
                                            "Butler"))
(flo_grass_sav_h_plot_diff <- ggplot(flo_grass_sav_h_plot_merge, 
                             aes(xmin=xmin, xmax=xmax,
                                 ymax=value,ymin=0,
                                 group=variable,
                                 fill=variable, 
                                 color=variable,
                                 alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Flooded grassland and savanna\n") +
    theme(legend.title = element_blank())+
    ylim(-25,50)+
    xlim(0,25))

      # sla stdev ---- 
flo_grass_sav_slatd_h <- j_flo_g_sav_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(flo_grass_sav_slatd_h_plot <- ggplot(flo_grass_sav_slatd_h, 
                               aes(x=sla_std,group=dataset,
                                   fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
flo_grass_sav_slatd_h_data <- as.data.table(ggplot_build(flo_grass_sav_slatd_h_plot)$data[1])
flo_grass_sav_slatd_h_data <- flo_grass_sav_slatd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
flo_grass_sav_slatd_h_data_c <- flo_grass_sav_slatd_h_data[group==2]
flo_grass_sav_slatd_h_data_b <- flo_grass_sav_slatd_h_data[group==1]
flo_grass_sav_slatd_h_merge <- merge(flo_grass_sav_slatd_h_data_c, 
                                     flo_grass_sav_slatd_h_data_b, 
                              by=c("xmin","xmax"), 
                              suffixes = c(".c",".b"), allow.cartesian=TRUE)
flo_grass_sav_slatd_h_merge <- flo_grass_sav_slatd_h_merge[,Difference:=count.c-count.b]
setnames(flo_grass_sav_slatd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
flo_grass_sav_slatd_h_merge <- melt(flo_grass_sav_slatd_h_merge, 
                             id.vars = c("xmin","xmax"), 
                             measure.vars = c("Cardamom","Difference",
                                              "Butler"))
(flo_grass_sav_slatd_h_diff <- ggplot(flo_grass_sav_slatd_h_merge, 
                               aes(xmin=xmin, xmax=xmax,
                                   ymax=value,ymin=0,
                                   group=variable,
                                   fill=variable, 
                                   color=variable,
                                   alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Flooded grassland and savanna\n") +
    theme(legend.title = element_blank())+
    ylim(-70, 70)+
    xlim(0, 65))

# 14) tropical and subtropical grassland, savanna, shrubland ----
      # sla mean ----
trpsbtrp_grass_sav_shr_h <- j_trpsbtrp_g_sav_shr_sla %>%
  gather(dataset, sla, -x, -y)
trpsbtrp_grass_sav_shr_h_plot <- ggplot(trpsbtrp_grass_sav_shr_h,
                                         aes(sla, group = dataset,
                                             fill = dataset))+
  geom_histogram(bins = 100, alpha = 0.4)
trpsbtrp_grass_sav_shr_h_data <- as.data.table(ggplot_build(trpsbtrp_grass_sav_shr_h_plot)$data[1])
trpsbtrp_grass_sav_shr_h_data <- dplyr::select(trpsbtrp_grass_sav_shr_h_data,
                                               count,xmin,xmax,ymin,ymax,group)
trpsbtrp_grass_sav_shr_h_data_c <- trpsbtrp_grass_sav_shr_h_data[group == 1]
trpsbtrp_grass_sav_shr_h_data_b <- trpsbtrp_grass_sav_shr_h_data[group == 2]
trpsbtrp_grass_sav_shr_h_merge <- merge(trpsbtrp_grass_sav_shr_h_data_c, 
                                        trpsbtrp_grass_sav_shr_h_data_b, 
                                        by=c("xmin","xmax"), 
                                        suffixes = c(".c",".b"), 
                                        allow.cartesian=TRUE)
trpsbtrp_grass_sav_shr_h_merge <- trpsbtrp_grass_sav_shr_h_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_grass_sav_shr_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom","Butler"))
trpsbtrp_grass_sav_shr_h_merge <- melt(trpsbtrp_grass_sav_shr_h_merge, 
                           id.vars = c("xmin","xmax"), 
                           measure.vars = c("Cardamom","Difference",
                                            "Butler"))
(trpsbtrp_grass_sav_shr_h_diff <- ggplot(trpsbtrp_grass_sav_shr_h_merge, 
                             aes(xmin=xmin, xmax=xmax,
                                 ymax=value,ymin=0,
                                 group=variable,
                                 fill=variable, 
                                 color=variable,
                                 alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area (m2.kg-1)")+
  #  ylab("Count\n")+
    ggtitle("Tropical and subtropical grassland, savanna and shrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-80,150)+
    xlim(0,40))

      # sla stdev ----
trpsbtrp_grass_sav_shr_slastd_h <- j_trpsbtrp_g_sav_shr_slastd %>%
  gather(key = "dataset",value="sla_std", -x,-y)
(trpsbtrp_grass_sav_shr_slastd_h_plot <- ggplot(trpsbtrp_grass_sav_shr_slastd_h, 
                                      aes(x=sla_std,group=dataset,
                                          fill=dataset))+
    geom_histogram(bins = 100, alpha=0.4)+
    theme_ipsum())
trpsbtrp_grass_sav_shr_slastd_h_data <- as.data.table(ggplot_build(trpsbtrp_grass_sav_shr_slastd_h_plot)$data[1])
trpsbtrp_grass_sav_shr_slastd_h_data <- trpsbtrp_grass_sav_shr_slastd_h_data %>%
  dplyr::select(count, xmin, xmax, group) 
trpsbtrp_grass_sav_shr_slastd_h_data_c <- trpsbtrp_grass_sav_shr_slastd_h_data[group==2]
trpsbtrp_grass_sav_shr_slastd_h_data_b <- trpsbtrp_grass_sav_shr_slastd_h_data[group==1]
trpsbtrp_grass_sav_shr_slastd_h_merge <- merge(trpsbtrp_grass_sav_shr_slastd_h_data_c, 
                                               trpsbtrp_grass_sav_shr_slastd_h_data_b, 
                                     by=c("xmin","xmax"), 
                                     suffixes = c(".c",".b"), allow.cartesian=TRUE)
trpsbtrp_grass_sav_shr_slastd_h_merge <- trpsbtrp_grass_sav_shr_slastd_h_merge[,Difference:=count.c-count.b]
setnames(trpsbtrp_grass_sav_shr_slastd_h_merge, old = c("count.c","count.b"),
         new = c("Cardamom", "Butler"))
trpsbtrp_grass_sav_shr_slastd_h_merge <- melt(trpsbtrp_grass_sav_shr_slastd_h_merge, 
                                    id.vars = c("xmin","xmax"), 
                                    measure.vars = c("Cardamom","Difference",
                                                     "Butler"))
(trpsbtrp_grass_sav_shr_slastd_h_diff <- ggplot(trpsbtrp_grass_sav_shr_slastd_h_merge, 
                                      aes(xmin=xmin, xmax=xmax,
                                          ymax=value,ymin=0,
                                          group=variable,
                                          fill=variable, 
                                          color=variable,
                                          alpha = 0.7))+
    geom_rect()+
    theme_classic()+
    scale_fill_viridis(discrete = TRUE)+
    scale_color_manual(values=c("black","black","black"))+
    xlab("\nSpecific Leaf Area StDev (m2.kg-1)")+
    ylab("Count\n")+
    ggtitle("Tropical and subtropical grassland, savanna and shrubland\n") +
    theme(legend.title = element_blank())+
    ylim(-200, 200)+
    xlim(0, 70))

## HISTOGRAM DIFF PANELLED FOR 14 BIOMES sla mean ----
  # major biomes 1
(panel_majbiome_hist_diff <- ggarrange(taiga_sla_hist_diff_plot,
                                    tundra_sla_hist_diff_plot,
                                    tmp_b_m_sla_h_diff, 
                                    trpsbtrp_m_broad_h_diff,
                                    ncol = 2,
                                    nrow = 2, common.legend = TRUE,
                                    align = "hv"))
ggsave("./figures/panel_majbiome_diff_hist1.png", panel_majbiome_hist_diff, 
       width = 40, height = 25, units = "cm", dpi = 500)

  # major biomes 2
(panel_majbiome_hist_diff_2 <- ggarrange(temp_con_sla_hist_diff_plot,
                                      trpsbtrp_d_broad_h_diff,
                                      trpsbtrp_conif_h_diff,
                                      med_f_w_scr_h_plot_diff,
                                      ncol = 2, 
                                      nrow = 2,
                                      common.legend = TRUE,
                                      align = "hv"))
ggsave("./figures/panel_majbiome_diff_hist2.png", panel_majbiome_hist_diff_2,
       width = 40, height = 25, units = "cm", dpi = 500)

  # other biomes 1
(panel_obiome1_hist_diff <- ggarrange(des_x_scr_h_plot_merge_diff,
                                     temp_grass_sav_shr_h_diff,
                                     trpsbtrp_grass_sav_shr_h_diff,
                                     ncol = 2,
                                     nrow = 2,
                                     common.legend = TRUE,
                                     align = "v"))
ggsave("./figures/panel_obiome1_diff_hist.png", panel_obiome1_hist_diff,
       width = 45, height = 20, units = "cm", dpi = 500)

  # other biomes 2
(panel_obiome2_hist_diff <- ggarrange(mont_grass_shr_h_diff,
                                      mangr_h_plot_diff,
                                      flo_grass_sav_h_plot_diff,
                                      ncol = 2,
                                      nrow = 2,
                                      common.legend = TRUE,
                                      align = "hv"))
ggsave("./figures/panel_obiome2_diff_hist.png", panel_obiome2_hist_diff,
       width = 45, height = 20, units = "cm", dpi = )

## HISTOGRAM DIFF PANELLED FOR 14 BIOMES sla stdev ----
  # major biomes 1
(panel_majbiome1_hist_slastd_diff <- ggarrange(taiga_slastd_hist_diff_plot,
                                           tundra_slastd_h_diff,
                                           temp_c_slastd_h_diff,
                                           trpsbtrp_m_broad_slastd_h_diff,
                                           temp_broad_mix_slastd_h_diff,
                                           ncol = 2,
                                           nrow = 3,
                                           common.legend = TRUE,
                                           align = "v"))
ggsave("./figures/panel_majbiome1_hist_slastd_diff.png", 
       panel_majbiome1_hist_slastd_diff,
       width = 45, height = 20, units = "cm", dpi = 500) 

  # major biomes 2
(panel_majbiome2_hist_slastd_diff <- ggarrange(trpsbtrp_dry_broad_slastd_h_diff,
                                              trpsbtrp_conif_slastd_h_diff,
                                              med_f_w_scr_slastd_h_diff,
                                              ncol = 2,
                                              nrow = 2,
                                              common.legend = TRUE,
                                              align = "v"))
ggsave("./figures/panel_majbiome2_hist_slastd_diff.png", 
       panel_majbiome2_hist_slastd_diff,
       width = 45, height = 20, units = "cm", dpi = 500)

  # other biomes 1
(panel_obiome1_hist_slastd_diff <- ggarrange(des_x_scr_slastd_h_diff,
                                            temp_grass_sav_shr_slastd_h_diff,
                                            trpsbtrp_grass_sav_shr_slastd_h_diff,
                                            ncol = 2,
                                            nrow = 2,
                                            common.legend = TRUE,
                                            align = "v"))

ggsave("./figures/panel_obiomes1_hist_slastd_diff.png", panel_obiome_hist_slastd_diff,
       width = 45, height = 20, units = "cm", dpi = 500)

  # other biomes 2
(panel_obiome2_hist_slastd_diff <- ggarrange(mont_grass_shr_slastd_h_diff,
                                             mangr_slastd_h_diff,
                                             flo_grass_sav_slatd_h_diff,
                                            ncol = 2,
                                            nrow = 2,
                                            common.legend = TRUE,
                                            align = "v"))

ggsave("./figures/panel_obiomes2_hist_slastd_diff.png", panel_obiome_hist_slastd_diff,
       width = 45, height = 20, units = "cm", dpi = 500)



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
(heatsc_trpsbtrp_con_slastd <- heatscatter(trpsbtrp_con_std_c_n, 
                                           trpsbtrp_con_std_b_n, 
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
# 1) taiga ----
  # sla mean 
taiga_sla_stat <- j_taiga_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = Metrics::bias(butler,new_b_val),
         rmse_av = Metrics::rmse(butler, new_b_val),
         rmse_row = sqrt(Metrics::se(butler, new_b_val)))

  # sla stdev
taiga_slastd_stat <- j_taiga_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))

# 2) tundra ----
  # sla mean 
tundra_sla_stat <- j_tundra_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = bias(butler, new_b_val),
         rmse_av = rmse(butler, new_b_val),
         rmse_row = rmse(butler, new_b_val))

  # sla stdev 
tundra_slastd_stat <- j_tundra_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))

# 3) temp conif forest ----
  # sla mean
temp_con_sla_stat <- j_tmp_c_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = bias(butler, new_b_val),
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)))

  # sla stdev 
temp_con_slastd_stat <- j_tmp_c_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))

# 4) temp broad mix forest ----
  # sla mean 
temp_b_m_sla_stat <- j_tmp_b_m_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = bias(butler, new_b_val),
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)))

  # sla stdev
temp_b_m_slastd_stat <- j_tmp_b_m_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))

# 5) tropical and subtropical dry broadleaf ----
  # sla mean 
trpsbtrp_d_b_sla_stat <- j_trp_sbtrp_d_b_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = bias(butler, new_b_val),
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)))

  # sla stdev 
trpsbtrp_d_b_slastd_stat <- j_trp_sbtrp_d_b_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))

# 6) tropical and subtropical conif forest ----
  # sla mean 
trpsbtrp_con_sla_stat <- j_trp_sbtrp_c_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = bias(butler, new_b_val),
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)))

  # sla stdev
trpsbtrp_con_slastd_stat <- j_trp_sbtrp_c_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))

# 7) tropical subtropical moist broadleaf ----
  # sla mean 
trpsbtrp_m_b_sla_stat <- j_trp_sbtrp_m_br_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = bias(butler, new_b_val),
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)))

  # sla stdev 
trpsbtrp_m_b_slastd_stat <- j_trp_sbtrp_m_br_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))

# 8) mediterranean forests, woodlands, scrub ----
  # sla mean
med_f_w_scr_sla_stat <- j_med_f_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area) %>%
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
         bias = bias(butler, new_b_val),
         rmse_av = rmse(butler, new_b_val),
         rmse_row = sqrt(se(butler, new_b_val)))

  # sla stdev 
med_f_w_scr_slastd_stat <- j_med_f_slastd %>%
  rename("cardamom_std" = Standard_Deviation, 
         "butler_std" = specific.leaf.area) %>%
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
         bias = bias(butler_std, new_bstd_val),
         rmse_av = rmse(butler_std, new_bstd_val),
         rmse_row = sqrt(se(butler_std, new_bstd_val)))


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

#### CHECK FOR LM INCLUDING THE BIOMES ####


