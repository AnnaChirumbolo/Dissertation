###############################################################################
#                 BSc Hons Dissertation                                       #
#                 Assessing degree of consistency between                     # 
#                 functional leaf trait estimates from diff methods           #
###############################################################################

## things to do 
  
  # first observations
# maps - add scales
# diff map - change scale colouring
# scatter plot - make heatmap 
# correlation chart - use density maps to find difference between histograms, to be calculated 

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

library(ncdf4)
library(RColorBrewer)
library(ggplot2)
library(raster)
library(tidyverse)
library(SimDesign)
library(Metrics)
library(gplots)
library(reshape2)
library(purrr)
library("grid")
library(diffeR)
library(PerformanceAnalytics)
library(viridis)
library(hrbrthemes)
library(lattice)

### opening netcdf, to data frame ----
cardamom_sla <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc", varname="sla")
cardamom_sla_std <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc", varname="Standard_Deviation")

cardamom_sla_df <- raster::as.data.frame(cardamom_sla, xy = TRUE)
cardamom_sla_std_df <- raster::as.data.frame(cardamom_sla_std, xy=TRUE)

butler_sla <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", varname="sla")
butler_sla_std <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", varname ="sla_std")
butler_sla_df <- raster::as.data.frame(butler_sla, xy = TRUE) 
butler_sla_std_df <- raster::as.data.frame(butler_sla_std, xy=TRUE)

### new df manipulation ----

  # CARDAMOM 1x1
cardamom_df <- cardamom_df %>%
  filter(sla!=0)
cardamom_df <- cardamom_df %>%
  mutate(degree=seq.int(nrow(cardamom_df)))
# see what min and max values are for sla in cardamom
setMinMax (cardamom_sla[[1]]) # values: 3.239088, 62.14139  (min, max)
cardamom_sla_layer[cardamom_sla_layer > 40] <- NA

  # BUTLER 1x1

butler_nc <- butler_nc %>%
  rename("sla"=specific.leaf.area)
butler_nc <- butler_nc %>%
  filter(specific.leaf.area!=0)
butler_nc <- butler_nc %>%
  mutate(degree = seq.int(nrow(butler_nc)))
# see what min and max values are for sla in butler
setMinMax(butler_sla[[1]]) #values: 5.076668, 32.43943  (min, max)

# butler <- read.csv("Butler_sla.csv") # this dataset is at 0.5x0.5 res, 
# don't need it
# butler <- butler %>%
#  rename("sla.m2.kg"=spat_1_pft_sla, "sla_std" = spat_1_pft_sla_sd)
# butler <- write.csv(butler,"Butler_sla.csv")

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

png("./figures/plot_sla_SAMESCALE.png", width = 50, height = 20, units = "cm", res = 200)
par(mfrow=c(1,2), oma = c(0,3,8,0) + 0.1, mar = c(7,0,2,8) + 0.1)
plot(cardamom_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")), zlim=c(0,63),
     main="Cardamom\n")
plot(butler_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")), zlim=c(0,63), 
     main="Butler\n",
     legend.args=list(text='\nSpecific Leaf Area (m2.kg-1)', side=4, font=2, line=2.3))
dev.off() # saved to folder 

### Joining datasets ----

joined_sla <- left_join(cardamom_sla_df, butler_sla_df) 
joined_sla_std <- left_join(cardamom_sla_std_df, butler_sla_std_df)

joined_sla <- joined_sla %>%
  rename("cardamom" = sla, "butler" = specific.leaf.area)
joined_sla_std <- joined_sla_std %>%
  rename("cardamom_std" = Standard_Deviation, "butler_std" = specific.leaf.area)

  # removing the NAs from the joine datasets 
joined_sla_noNA <- joined_sla %>%
  filter(cardamom!=0&butler!=0)
joined_sla_nocoord <- joined_sla %>%
  select(-x,-y) 
joined_sla_nocoordNA <- joined_sla_nocoord %>%
  filter(cardamom!=0 & butler != 0)

joined_sla_std_noNA <- joined_sla_std %>%
  filter(cardamom_std !=0 & butler_std != 0)
joined_sla_std_nocoord <- joined_sla_std %>%
  select(-x,-y)
joined_sla_std_nocoordNA <- joined_sla_nocoord %>%
  filter(cardamom!=0 & butler != 0)

### 1) DATA EXPLORATION 1 - MAKING SCATTERPLOT (CORRELATION) AND HEATMAP ----

  # scatterplot (correlation chart) 
corr_char_sla <- chart.Correlation(joined_sla_nocoord, histogram=TRUE, pch=19)
corr_char_sla_std <- chart.Correlation(joined_sla_std_nocoord, histogram = TRUE, pch=19)
# both too many data points to be sure of what's going on - too clustered 

  # 2D density chart (correlation charts)

# raster function 
(heatmap_trial <- ggplot(joined_sla_nocoord, aes(cardamom, butler))+
    stat_density_2d(aes(fill=..density..), geom = "raster", contour = FALSE)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)))
ggsave("./figures/2d_density.png", last_plot(), width = 30, height = 20, units="cm", 
       dpi = 300)

(heatmap_trial_1 <- ggplot(joined_sla_noNA, aes(cardamom, butler))+
    stat_density_2d(aes(fill=..level..), geom = "polygon")+
    theme_classic())
ggsave("./figures/2d_density_polygon.png", last_plot(), width = 30, height = 20,
       units = "cm", dpi = 300)

# hexbin version of 2d map 

(hexbin_map_corr <- ggplot(joined_sla, aes(x=cardamom, y=butler) ) +
  geom_hex(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_classic()+
  xlab("\nCardamom")+
  ylab("Butler\n"))
ggsave("./figures/hexbin_map-count200.png", hexbin_map_corr, width = 30, 
       height = 20, units = "cm", dpi = 300)

# heatmap (with lattice package and viridis package for colour palette)

  # modifying the dataset (joined_sla) to wide format for the matrix
colnames(joined_sla) <- joined_sla[1,]
rownames(joined_sla_noNA) <- joined_sla_noNA[,1]

joined_sla_matrix <- joined_sla_noNA %>% 
  gather(dataset, sla,-x,-y) %>%
  mutate(x_ = make.unique(as.character(x))) %>%
  mutate(y_ = make.unique(as.character(y))) %>%
  select(-x,-y, -dataset)

rownames(joined_sla_matrix) <- joined_sla_matrix[,3]
joined_sla_matrix <- joined_sla_matrix[,-3]

joined_sla_matrix <- joined_sla_matrix %>%
  spread(x_,sla)

levelplot(joined_sla_nocoordNA, col.regions = terrain.colors(100)) # try cm.colors() or terrain.colors()
View(joined_sla)

### DENSITY PLOT ### ----

new_joined_density <- new_joined %>%
  gather(key="dataset", value="sla",-y,-x)

(density_plot <- ggplot(new_joined_density,aes(x=sla, group=dataset,
                                               fill=dataset)) +
  geom_density(adjust=1.5, alpha=0.4)+
  theme_ipsum()+
  scale_fill_discrete(name = "Specific Leaf Area", 
                      labels = c("Cardamom", "Butler")))
## one way of visualising the distribution of data by themselves 

### scatterplot attempt 2 ----

(degree_scatter <- ggplot()+
   geom_point(cardamom_df, mapping =  aes(x,sla,color="red"))+
   geom_point(butler_nc,mapping=aes(x,specific.leaf.area, color="blue")))


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

joined_wide <- dcast(new_joined,dataset+y~x, value.var="sla")

joined_wide <- as.matrix(joined_wide)
rownames(joined_wide) <- joined_wide[, 1]
joined_wide <- joined_wide[, -1];
str(joined_wide)
heatmap.2(joined_wide,Colv = NA, Rowv = NA)


### COMPARING RESULT ESTIMATES ---- 

difference_data <- map2(cardamom_df, butler_nc, setdiff) %>% 
  map_int(length)

bias <- Metrics::bias(cardamom_sla, butler_sla) ## let's try it?? ASK CC?
bias # gives overall rel bias, and because there's lots of NAs it gives NA as result?
butler_sla <- butler_nc$sla  
butler_sla <- as.data.frame(butler_sla) 
butler_sla
cardamom_sla <- cardamom_df$sla
cardamom_sla <- as.data.frame(cardamom_sla) 

## RMSE function??? ----

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

rmse <- RMSE(butler_nc$sla, cardamom_df$sla)
rmse

## PACKAGE differR ----

categoryComponentsPlot(cardamom_sla, butler_sla, units = "m^2/kg")

# Cross-tabulate two RasterLayer objects, or mulitiple layers in a RasterStack 
# or RasterBrick to create a contingency table.
crosstab <- crosstab(cardamom_sla,butler_nc)
crosstab <- as.data.frame(crosstab)
heatmap(crosstab)

cardamom_sla[30,30]
crosstabm <- crosstabm(cardamom_sla, butler_nc, percent = TRUE)
heatmap(crosstabm)

sla_only <- new_joined %>%
  select(-x,-y) 
 
sla_only<- rcorr(sla_only, type="pearson")

## CORRELATION CHART BETWEEN SLA PARAMETERS ----

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

# Get some colors
sla_only <- as.matrix(sla_only)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = sla_only, col = col, symm = TRUE)


