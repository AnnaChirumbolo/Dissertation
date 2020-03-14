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
# calculation of bias / RMSE / R2 - visual and tabular repr. r2 in the negative values...normal?

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

## category Components plot ----

categoryComponentsPlot(cardamom_sla, butler_sla, units = "m^2/kg")

# Cross-tabulate two RasterLayer objects, or mulitiple layers in a RasterStack ----
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



## R SQUARED CALCULATIONS AND CSV OUTPUT ----
  
  ## steps to calculate R2:
  #rss <- sum((preds - actual) ^ 2)  ## residual sum of squares
  #tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  #rsq <- 1 - rss/tss

  # SLA mean
sla_rsq <- joined_sla %>%
  filter(cardamom!=0, butler!=0) %>%
  summarise(rss = sum((cardamom - butler)^2),     ## residual sum of squares
            tss = sum((butler - mean(butler))^2), ## total sum of squares
            rsq = 1- rss/tss)
  # SLA STDEV
slastd_rsq <- joined_sla_std %>%
  filter(cardamom_std!=0, butler_std!=0) %>%
  summarise(rss = sum((cardamom_std - butler_std)^2),
            tss = sum((butler_std - mean(butler_std))^2),
            rsq = 1- rss/tss)

# a negative R2 means that the chosen model (with its constraints) fits the 
# data really poorly.

rsq_results <- rbind(sla_rsq, slastd_rsq)

#write.csv(rsq_results, "R2_results.csv")

### table with r2 and rmse ----

stat_world <- cbind(rsq_results,rmse_all)
stat_world <- stat_world %>%
  dplyr::select(-tss,-rss) %>%
  mutate(sla=c("Mean", "StDev"))

formattable(stat_world) # to create the table

## the R2 when im doing a linear regression... ----

lm_sla <- lm(butler ~ cardamom, data = joined_sla)
plot(lm_sla)
summary(lm_sla)$r.squared
# Multiple R-squared:  0.0003106,	Adjusted R-squared:  0.0002266 

lm_slastd <- lm(butler_std  ~ cardamom_std, data = joined_sla_std)
plot(lm_slastd)
summary(lm_slastd) 
# Multiple R-squared:  0.0002966,	Adjusted R-squared:  0.0002127 

######################splitting of the world####################################

## do splitting by latitudes ----

  # tropics SLA MEAN ----
trpmat <- matrix(data <- c(-180,180,-23.5,23.5), nrow = 2, ncol = 2,
                 byrow = TRUE)
trpext <- extent(trpmat)
trp <- crop(cardamom_sla, trpext)
plot(trp[[1]]) # tropical lats
  # tropics SLA STD ----
trpstd<- crop(cardamom_sla_std, trpext)
plot(trpstd[[1]], asp=NA) # height to fix when (and if) saving it as png
  
  # temperate ----
tmpmat <- matrix(data <- c(-180,180,-66.5,66.5), nrow = 2, ncol = 2, 
                 byrow = TRUE)
tmpext <- extent(tmpmat)
tmp <- crop(cardamom_sla, tmpext)
plot(tmp[[1]])
trp_na[] <- NA
tmp <- mask(tmp,trp_na) 




