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

### opening netcdf, to data frame ----
cardamom_sla <- raster("CARDAMOM_2001_2010_LCMA_zeros.nc", varname="sla")
cardamom_sla_std <- raster("CARDAMOM_2001_2010_LCMA_zeros.nc", varname="Standard_Deviation")

cardamom_sla_df <- raster::as.data.frame(cardamom_sla, xy = TRUE)
cardamom_sla_std_df <- raster::as.data.frame(cardamom_sla_std, xy=TRUE)

butler_sla <- raster("Butler_Leaftraits_Processed_1x1_zeros.nc", varname="sla")
butler_sla_std <- raster("Butler_Leaftraits_Processed_1x1_zeros.nc", varname ="sla_std")

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

png("plot_sla_DIFFSCALES.png", width = 50, height = 20, units = "cm", res = 200)
par(mfrow=c(1,2), oma = c(0,3,8,0) + 0.1, mar = c(7,0,2,8) + 0.1)
plot(cardamom_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")))
plot(butler_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")))
dev.off()

### visualising sla from both datasets SAME SCALE ----

png("plot_sla_SAMESCALE.png", width = 50, height = 20, units = "cm", res = 200)
par(mfrow=c(1,2), oma = c(0,3,8,0) + 0.1, mar = c(7,0,2,8) + 0.1)
plot(cardamom_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")), zlim=c(0,63))
plot(butler_sla[[1]], asp=NA, col = rev(brewer.pal(10, "RdBu")), zlim=c(0,63))
dev.off() # saved to folder 

### SCATTERPLOT - joining datasets ----

new_joined <- left_join(cardamom_df, butler_df) 

new_joined <- new_joined %>%
  filter(sla!=0 & specific.leaf.area!=0)

###### LET'S MAKE THE SCATTERPLOT - ATTEMPT 1 ###### ----

(map_scatter <- ggplot()+
  # geom_point(cardamom_df, mapping=aes(sla,y,color="red", alpha=0.7))+
  #geom_point(butler_nc, mapping=aes(specific.leaf.area,y,color="blue",alpha=0.7))+
   geom_point(cardamom_df,mapping = aes(x,y,color="red",alpha=0.7))+
   facet_grid(sla~.,scale="free"))
   
geom_point(butler_nc, mapping = aes(x,specific.leaf.area,color="blue",alpha=0.7))
 #  ylim(-90,+90)+
 #  xlim(-180,+180))# not ideal... 

## need to use degree 1x1 as my factor for correlation - for every degree of lat/lon
## there is a value of sla that is correlated and that is from both datasets...
## TO FIGURE OUT TOMORROW.

plot(cardamom_df$sla, col="red" )
par(new=TRUE)
plot(butler_nc$specific.leaf.area, col="green")

butler$butler_sla <- butler$sla.m2.kg
cardamom_df$cardamom_sla <- cardamom_df$sla
butler <- butler %>%
  select(-sla.m2.kg)
cardamom_df <- cardamom_df %>%
  select(-sla)s

joined_sla <- joined_sla%>%
  rename("b_sla" = sla.m2.kg, "c_sla"=sla) 

joined_sla <- joined_sla %>%
  mutate(c_sla1=c_sla, b_sla1=b_sla) %>%
  gather(key=dataset, value=sla1, -c_sla,-b_sla)

col = rep('black', length(new_joined))
col[new_joined!=butler_sla] <- 'blue'

(sla_scatterplot <- ggplot(new_joined, aes(x=sla,y=butler_sla))+
    geom_point(stat = "identity") +
    theme_classic())


### DENSITY PLOT ### ----
library(viridis)
install.packages("hrbrthemes")
library(hrbrthemes)

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
install.packages("diffeR")
library(diffeR)

cardamom_sla <- raster("CARDAMOM_2001_2010_LCMA_zeros.nc", varname="sla")
butler_nc <- raster("Butler_Leaftraits_Processed_1x1_zeros.nc", varname="sla")

categoryComponentsPlot(cardamom_sla, butler_nc, units = "m^2/kg")

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
sla_only <- sla_only %>%
  rename("cardamom"=sla, "butler"=specific.leaf.area)
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)

png("correlation_chart.png", width=30, height=20, units = "cm", res = 300)
corr_char <- chart.Correlation(sla_only, histogram=TRUE, pch=19)
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


