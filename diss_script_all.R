###############################################################################
#                 BSc Hons Dissertation                                       #
#                 Assessing degree of consistency between                     # 
#                 functional leaf trait estimates from diff methods           #
###############################################################################

## things to do 
  
# MODEL = BUTLER 
# OBSERVATIONS = CARDAMOM = REFERENCE

# scatter plot - make heatmap --| ASK CODING CLUB / MADE A HEATSCATTER (BETTER REPR) / ASK IF I COULD DO IT WITH THE MAP OF THE WORLD?

# for results (to put in table and have the figures showing)
# matrix of map, see how they differ from each other by each grid point 
# calculation of bias / RMSE / R2 - visual and tabular repr.

## Packages ####
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
## libraries ####
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
library(rworldmap)
library(rgeos)
library(maptools)
library(cleangeo)
library(tiff)
library(colorspace)
library(dichromat)

## custom functions ####

##OVERLAP PLOTS (modif from overlapping package) 
my.final.plot.sla <- function (x, OV = NULL){
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
    geom_density(aes(fill = AREA$group), alpha = 0.35) + 
    xlab("\nSLA mean (m2.kg-1)") + 
    ylab("")+
    theme_classic()+
    theme(legend.title = element_blank())+
    scale_color_brewer(palette = "Set1")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))
}

my.overlap.sla <- function (x, nbins = 1024, plot = FALSE, partial.plot = FALSE, 
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
          geom_line() + geom_line(aes(x, dd2$y2)) + 
          geom_line(aes(x,dd2$ovy), color = "red") + 
          geom_line(aes(x,dd2$ally), color = "blue") + 
          ggtitle(comptitle) + 
          xlab("") + 
          ylab("") + 
          theme(plot.title = element_text(hjust = 0.5),
                legend.title = element_blank())
        print(gg)
      }
      DD <- rbind(DD, dd2)
      COMPTITLE <- c(COMPTITLE, comptitle)
    }
  }
  names(xpoints) <- names(OV) <- COMPTITLE
  if (plot) 
    print(my.final.plot.sla(x, OV))
  return(list(DD = DD, OV = OV, xpoints = xpoints))
}

my.final.plot.std <- function (x, OV = NULL){
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
    geom_density(aes(fill = AREA$group), alpha = 0.35) + 
    xlab("\nSLA StDev (m2.kg-1)") + 
    ylab("")+
    theme_classic()+
    theme(legend.title = element_blank())+
    scale_color_brewer(palette = "Set1")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)) 
}
my.overlap.std <- function (x, nbins = 1024, plot = FALSE, partial.plot = FALSE, 
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
          geom_line() + geom_line(aes(x, dd2$y2)) + 
          geom_line(aes(x,dd2$ovy), color = "red") + 
          geom_line(aes(x,dd2$ally), color = "blue") + 
          ggtitle(comptitle) + 
          xlab("") + 
          ylab("") + 
          theme(plot.title = element_text(hjust = 0.5),
                legend.title = element_blank())
        print(gg)
      }
      DD <- rbind(DD, dd2)
      COMPTITLE <- c(COMPTITLE, comptitle)
    }
  }
  names(xpoints) <- names(OV) <- COMPTITLE
  if (plot) 
    print(my.final.plot.std(x, OV))
  return(list(DD = DD, OV = OV, xpoints = xpoints))
}

## extracting p-values from lm summary (credits to Stephen Turner: 
#https://gettinggeneticsdone.blogspot.com/2011/01/rstats-function-for-extracting-f-test-p.html)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    if (modelobject[["df.residual"]]!=0){ # if the residual is not = 0 
      f <- summary(modelobject)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      attributes(p) <- NULL
      return(p)
    }
  else
    modelobject <- NULL
}

# creating function to perform calculations of stats iteratively
lm.f <- function(x){
  lm <- lapply(x, function(dat) lm(dat[,4] ~ dat[,3],data=dat))
}

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
    df$rmse_row <- sqrt((df$butler - df$cardamom)^2)
    df
  })
}

# function to covert raster to data frame
mask.to.df <- function(x){
  new.list <- list()
  for (i in 1:length(x)){
    df <- raster::as.data.frame(x[[i]], xy =TRUE)
    name.df <- paste("df",names(x)[i],sep = ".")
    new.list[[name.df]] <- df
  }
  new.list
}

# joining dataframes 
join.f <- function(x,y){
  new.list <- list()
  join <- mapply(left_join, x, y,SIMPLIFY = FALSE)
  join <- lapply(join, na.omit)
  name.df <- names(x)[i]
  new.list[[name.df]] <- join
}

# masking the biomes by continent 
mask.biome.f <- function(x,y){
  new.list <- list()
  for (i in 1:length(x)){
    for (j in 1:length(y)){
      mask.biome <- raster::mask(x[[i]],y[[j]])
      if (!is.infinite(mask.biome@data@min)&!is.infinite(mask.biome@data@max)){
        name <- paste(names(x)[i],names(y)[j],sep = " ")
        new.list[[name]] <- mask.biome
      }
      else
        NULL
    }
  }
  new.list
}

# function for kableExtra, from: Michael Harper https://stackoverflow.com/questions/28166168/how-to-change-fontface-bold-italics-for-a-cell-in-a-kable-table-in-rmarkdown
format_cells <- function(df, rows ,cols, value = c("italics", "bold", "strikethrough")){
  
  # select the correct markup
  map <- setNames(c("*", "*", "~~"), c("italics", "bold", "strikethrough"))
  markup <- map[value]  
  
  for (r in rows){
    for(c in cols){
      
      # Make sure values are not factors
      df[[c]] <- as.character( df[[c]])
      
      # Update formatting
      df[r, c] <- paste0(markup, df[r, c], markup)
    }
  }
  
  return(df)
}

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
                           varname= "sla_std")

butler_sla <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", 
                     varname="sla")
butler_sla_std <- raster("./DATA/Butler_Leaftraits_Processed_1x1_zeros.nc", 
                         varname="sla_std")

# also opening cardamom lcma estimate for mapping, for write-up
cardamom_lcma <- raster("./DATA/CARDAMOM_2001_2010_LCMA_zeros.nc",
                        varname="lcma")
lcma <-levelplot(cardamom_lcma, margin=F, par.settings=RdBuTheme)
butler.sla <- levelplot(butler_sla, margin=F, par.settings=RdBuTheme)
png("./figures/lcma.sla.png", width = 40, height = 30, units = "cm",
    res = 500)
(panel.lcma.sla <- ggarrange(lcma, butler.sla, ncol = 1, nrow =2))
text.ppt <- grid::grid.text('LCMA mean (gC.m-2)', rot=90,
                            y=unit(0.8, "npc"), 
                            x=unit(0.90, "npc"))
text.temp <- grid::grid.text('SLA mean (m2.kg-1)', rot=90,
                             y=unit(0.25, "npc"), 
                             x=unit(0.90, "npc"))

dev.off()

# making .nc files into data frames 
cardamom_sla_df <- raster::as.data.frame(cardamom_sla, xy = TRUE)
cardamom_sla_std_df <- raster::as.data.frame(cardamom_sla_std, xy=TRUE)
butler_sla_df <- raster::as.data.frame(butler_sla, xy = TRUE) 
butler_sla_std_df <- raster::as.data.frame(butler_sla_std, xy=TRUE)

# basic data manipulation
cardamom_sla_df <- cardamom_sla_df %>%
  rename("cardamom" = specific.leaf.area)
cardamom_sla_std_df <- cardamom_sla_std_df %>%
  rename("cardamom_std" = sla_std)
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
setMinMax(cardamom_sla_std[[1]]) #values     : 6.686308, 117.3275  (min, max)
setMinMax(butler_sla_std[[1]]) #values     : 2.167324, 17.66437  (min, max)
plot(cardamom_sla_std[[1]], asp=NA, col=rev(brewer.pal(10,"RdBu")), 
     zlim=c(0,120), xlab="\nLongitude", 
     ylab="Latitude")
plot(butler_sla_std[[1]], asp=NA, col=rev(brewer.pal(10,"RdBu")),
     zlim=c(0,120), xlab="\nLongitude",
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
gl_stip_locs <- as.data.frame(gl_stip_locs[gl_stip_locs[, "layer"] == 1,])
gl_stip_locs_95 <- (butler_sla[[1]]>=cardamom_25th[[1]])*
  (butler_sla[[1]]<=cardamom_95th[[1]])
gl_stip_locs_95 <- rasterToPoints(gl_stip_locs_95)
gl_stip_locs_95 <- as.data.frame(gl_stip_locs_95[gl_stip_locs_95[, 
                                                                 "layer"] ==1,])
diff_pc <- as.matrix(anti_join(gl_stip_locs_95,gl_stip_locs))

par(mar = c(4, 2.5, 0.5, 2.5),mfcol=c(1,1))
png("./figures/global_analysis/stippling_world.png", width = 40, height = 25, 
    units = "cm", res = 500)
plot(butler_sla, asp = NA, col = rev(brewer.pal(10, "RdYlBu")),
     xlab="Longitude", ylab="Latitude", 
     legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                        side=4, font=1, line=2.3),
     main="Butler Sla Mean (stippling)\n")
    # sub = "Base map with Butler SLA mean (m2.kg-1). Stippling (black): mean SLA
    # Butler values falling within 25pc-75pc confidence interval of Cardamom mean 
    # SLA values.Stippling (green): SLA mean Butler values falling within the 
    # 25pc-95pc confidence intervals of Cardamom mean SLA values.")
points(gl_stip_locs, pch = 18, cex=0.5)
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
png("./figures/global_analysis/stippling_std_world.png", width = 40, height = 25, 
    units = "cm", res = 500)
plot(butler_sla_std, asp = NA, col = rev(brewer.pal(10, "RdYlBu")),
     xlab="\nLongitude", ylab="Latitude", 
     legend.args = list(text="\n\nSLA StDev (m2.kg-1)", 
                        side=4, font=1, line=2.3),
     main="Butler Sla StDev (stippling)\n")
points(stip_locs_std, pch = 18, cex=0.5)
dev.off()

#### HEATSCATTERS ####
  # creating numeric vectors of sla mean and sla stdev to input in heatscatter ----
cardamom_sla_num <- cardamom_sla_df$cardamom
butler_sla_num <- butler_sla_df$butler
cardamom_sla_std_num <- cardamom_sla_std_df$cardamom_std
butler_sla_std_num <- butler_sla_std_df$butler_std

  # LSD package ----
png("./figures/global_analysis/heatscatter_sla_panelled.png", width = 50, height = 25, 
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


#### PERCENTAGE OVERLAP DENSITY PLOTS ####
# sla mean ----
joined_sla <- joined_sla %>%
  filter(cardamom!=0, butler!=0)
sla_n <- list(cardamom = joined_sla$cardamom,
       butler = joined_sla$butler) 
png("./figures/global_analysis/global_density_overlap.png", width = 40, height = 20,
    units = "cm", res = 400)
sla_overl<- my.overlap.sla(sla_n, plot = TRUE) 
dev.off()
    
# sla stdev ----
joined_sla_std <- joined_sla_std %>%
  filter(cardamom_std!=0, butler_std!=0)
slastd_n <- list(cardamom_std = joined_sla_std$cardamom_std,
                 butler_std = joined_sla_std$butler_std)
png("./figures/global_analysis/global_std_density_overlap.png", width = 40,
    height = 20, units = "cm", res = 400)
sla_std_overl <- my.overlap.std(slastd_n, plot = T)
dev.off()

### GLOBAL DATA STATS: r2, rmse, bias ----
# SLA MEAN
global.r2.lm  <- lm(butler ~ cardamom,data = joined_sla)
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
         p_val = lmp(global.r2.lm),
         rmse_av = rmse(butler, cardamom),
         rmse_row = sqrt(se(butler, cardamom)),
         bias_av = bias(butler, cardamom),
         bias_row = butler-cardamom)

# SLA STD
global_slastd_lm <- lm(butler_std ~ cardamom_std, data = joined_sla_std )
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
         sla_r2 = sum_sqrd_dist_std_b / sum_sqrd_diff_b,
         p_val = lmp(global_slastd_lm),
         rmse_av = rmse(butler_std, cardamom_std),
         rmse_row = sqrt(se(butler_std, cardamom_std)),
         bias_av = bias(butler_std, cardamom_std),
         bias_row = butler_std-cardamom_std)




####                          splitting of the world                          ####


#### SPLITTING BY LATITUDE ####

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
trpSLA <- left_join(tropicsSLA_df,trp_b_df, by = c("x","y")) 
trpSTD <- left_join(tropicsSTD_df, trp_b_std_df, by=c("x","y"))
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
sbtrp_c <- merge(sbtrpN_df, sbtrpS_df, by=c("x", "y", 
                                            "specific.leaf.area"), 
                 all=TRUE) 

  # sla std 
sbtrp_std_N <- crop(cardamom_sla_std, sbtrpextN)
sbtrpSTD_c_N <- raster::as.data.frame(sbtrp_std_N, xy=TRUE)
sbtrp_std_S <- crop(cardamom_sla_std, sbtrpextS)
sbtrpSTD_c_S <- raster::as.data.frame(sbtrp_std_S, xy=TRUE)
sbtrpSTD_c <- merge(sbtrpSTD_c_N, sbtrpSTD_c_S, 
                    by=c("x","y","sla_std"), all = TRUE)

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
sbtrp_joined_SLA <- left_join(sbtrp_c, sbtrp_b, by=c("x","y"))
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
tmp_c_df <- merge(tmp_c_N_df, tmp_c_S_df, by = c("x","y",
                                                 "specific.leaf.area"), 
                  all = TRUE)

  # sla std 
tmpSTD_c_N <- crop(cardamom_sla_std, tmpextN)
tmpSTD_c_S <- crop(cardamom_sla_std, tmpextS)
tmpSTD_c_N_df <- raster::as.data.frame(tmpSTD_c_N, xy =TRUE)
tmpSTD_c_S_df <- raster::as.data.frame(tmpSTD_c_S, xy =TRUE)
tmpSTD_c_df <- merge(tmpSTD_c_N_df, tmpSTD_c_S_df, 
                     by = c("x","y","sla_std"), 
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
tmp_sla <- left_join(tmp_c_df,tmp_b_df, by=c("x","y"))
tmp_slastd <- left_join(tmpSTD_c_df, tmpSTD_b_df, by=c("x","y"))
  
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
plN <- left_join(plN_c_df, plN_b_df, by=c("x","y"))
plN_std <- left_join(plN_slastd_c_df, plN_slastd_b_df,
                     by=c("x","y"))


#### PERCENTAGE OVERLAP DENSITY PLOTS ####
  # sla mean ----
lat.list <- list(Tropics=trpSLA,Subtropics=sbtrp_joined_SLA,
                 Temperate=tmp_sla, Pole=plN)
lat.list.n <- list()
lat.density.plot <- list()
for (i in 1:length(lat.list)){
  lat.list[[i]] <- lat.list[[i]]%>%
    filter(specific.leaf.area.x!=0,specific.leaf.area.y!=0) %>%
    rename("cardamom" = specific.leaf.area.x, 
           "butler" = specific.leaf.area.y) 
  list <- list(cardamom=lat.list[[i]][["cardamom"]],
               butler=lat.list[[i]][["butler"]])
  name <- paste(names(lat.list)[i],"n",sep = "_")
  lat.list.n[[name]] <- list
  for (j in 1:length(lat.list.n)){
    plot_name <- paste("./figures/lat_analysis/lat",
                       names(lat.list)[i],"density_overlap.png",
                       sep = "_")
    png(plot_name, width = 40, height = 20, units = "cm", res = 400)
    lat.density.plot[[j]] <- my.overlap.sla(lat.list.n[[j]], plot = TRUE)
    dev.off()
  }
}

  # sla stdev ----
lat.std.list <- list(Tropics.std=trpSTD,Subtropics.std=sbtrp_joined_STD,
                   Temperate.std=tmp_slastd, Pole.std=plN_std)
lat.std.list.n <- list()
lat.density.std.plot <- list()
for (i in 1:length(lat.std.list)){
  lat.std.list[[i]] <- lat.std.list[[i]]%>%
    filter(sla_std!=0,specific.leaf.area!=0) %>%
    rename("cardamom_std" = sla_std, 
           "butler_std" = specific.leaf.area) 
  list <- list(cardamom_std=lat.std.list[[i]][["cardamom_std"]],
               butler_std=lat.std.list[[i]][["butler_std"]])
  name <- paste(names(lat.std.list)[i],"n",sep = "_")
  lat.std.list.n[[name]] <- list
  for (j in 1:length(lat.std.list.n)){
    plot_name <- paste("./figures/lat_analysis/lat_std",names(lat.std.list)[i],
                       "density_overlap.png",sep = "_")
    png(plot_name, width = 40, height = 20, units = "cm", res = 400)
    lat.density.std.plot[[j]] <- my.overlap.std(lat.std.list.n[[j]], 
                                                plot = TRUE)
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
lat.list.raster <- list(tropics=trp_b,subtropics_N=sbtrpN_b,
                        subtropics_S=sbtrpS_b,temperate_N=tmp_b_N,
                        temperate_S=tmp_b_S,pole_N=plN_b)

for (i in 1:length(lat.list.raster)){
  crop.25 <- crop(cardamom_25th,extent(lat.list.raster[[i]]))
  crop.75 <- crop(cardamom_75th, extent(lat.list.raster[[i]]))
  crop.95 <- crop(cardamom_95th, extent(lat.list.raster[[i]]))
  stp.25 <- raster::mask(crop.25, lat.list.raster[[i]])
  stp.75 <- raster::mask(crop.75, lat.list.raster[[i]])
  stp.95 <- raster::mask(crop.95, lat.list.raster[[i]])
  stp.name25 <- paste("stp25", names(lat.list.raster)[i],sep = ".")
  stp.name75 <- paste("stp75", names(lat.list.raster)[i],sep = ".")
  stp.name95 <- paste("stp95",names(lat.list.raster)[i],sep=".")
  lat.stp.25[[stp.name25]] <- stp.25
  lat.stp.75[[stp.name75]] <- stp.75
  lat.stp.95[[stp.name95]] <- stp.95
  stp.locs.75.calc <- (lat.list.raster[[i]]>=lat.stp.25[[i]])*
    (lat.list.raster[[i]]<= lat.stp.75[[i]])
  stp.locs.95.calc <- (lat.list.raster[[i]]>=lat.stp.25[[i]])*
    (lat.list.raster[[i]]<= lat.stp.95[[i]])
  stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
  stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
  stp.locs.75.p <- as.data.frame(stp.locs.75.p[stp.locs.75.p[,
                                                             "layer"] == 1,])
  stp.locs.95.p <- as.data.frame(stp.locs.95.p[stp.locs.95.p[, 
                                                             "layer"] == 1,])
  stp.locs.75n <- paste("stp.locs",names(lat.list.raster)[i],"75",
                        sep = ".")
  stp.locs.95n <- paste("stp.locs",names(lat.list.raster)[i],"95",
                        sep = ".")
  stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
  stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
  diff.stp.locs.calc <- as.matrix(anti_join(stp.locs.95[[i]],
                                            stp.locs.75[[i]]))
  diff.stp.locs.n <- paste("diff",names(lat.list.raster)[i],sep = ".")
  diff.stp.locs[[diff.stp.locs.n]] <- diff.stp.locs.calc
  png(paste("./figures/lat_analysis/lat.stippling",
            names(lat.list.raster)[i],".png",
            sep = "."),
      width = 50,height = 25,units = "cm",res = 500)
  plot(lat.list.raster[[i]],col=rev(brewer.pal(10,"RdYlBu")),
       legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                          side=4, font=1, line=2.3),
       main="Butler Sla Mean (stippling)\n")
  points(stp.locs.75[[i]],pch=18,cex=0.8)
  if (length(diff.stp.locs[[i]])!=0){
  points(diff.stp.locs[[i]],pch=23,cex=0.9,col="darkgreen",bg="green")
  }
  else
    NULL
  dev.off()
}

  # sla stdev 
lat.stp.25 <- list()
lat.stp.75 <-list()
lat.stp.95 <- list()
stp.locs.75 <- list()
stp.locs.95<- list()
diff.stp.locs <- list()
lat.std.list.raster<- list(tropics.std=trp_b_std,subtropics.std_N=sbtrpSTD_b_N,
                     subtropics.std_S=sbtrpSTD_b_S,temperate.std_N=tmpSTD_b_N,
                     temperate.std_S=tmpSTD_b_S,pole.std_N=plN_slastd_b) 
for (i in 1:length(lat.std.list.raster)){
  crop.25 <- crop(cardamom_25th,extent(lat.std.list.raster[[i]]))
  crop.75 <- crop(cardamom_75th, extent(lat.std.list.raster[[i]]))
  crop.95 <- crop(cardamom_95th, extent(lat.std.list.raster[[i]]))
  stp.25 <- raster::mask(crop.25, lat.std.list.raster[[i]])
  stp.75 <- raster::mask(crop.75, lat.std.list.raster[[i]])
  stp.95 <- raster::mask(crop.95, lat.std.list.raster[[i]])
  stp.name25 <- paste("stp25", names(lat.std.list.raster)[i],sep = ".")
  stp.name75 <- paste("stp75", names(lat.std.list.raster)[i],sep = ".")
  stp.name95 <- paste("stp95",names(lat.std.list.raster)[i],sep=".")
  lat.stp.25[[stp.name25]] <- stp.25
  lat.stp.75[[stp.name75]] <- stp.75
  lat.stp.95[[stp.name95]] <- stp.95
  stp.locs.75.calc <- (lat.std.list.raster[[i]]>=lat.stp.25[[i]])*
    (lat.std.list.raster[[i]]<= lat.stp.75[[i]])
  stp.locs.95.calc <- (lat.std.list.raster[[i]]>=lat.stp.25[[i]])*
    (lat.std.list.raster[[i]]<= lat.stp.95[[i]])
  stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
  stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
  stp.locs.75.p <- as.data.frame(stp.locs.75.p[stp.locs.75.p[,
                                                             "layer"] == 1,])
  stp.locs.95.p <- as.data.frame(stp.locs.95.p[stp.locs.95.p[, 
                                                             "layer"] == 1,])
  stp.locs.75n <- paste("stp.locs",names(lat.std.list.raster)[i],"75",
                        sep = ".")
  stp.locs.95n <- paste("stp.locs",names(lat.std.list.raster)[i],"95",
                        sep = ".")
  stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
  stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
  diff.stp.locs.calc <- as.matrix(anti_join(stp.locs.95[[i]],
                                            stp.locs.75[[i]]))
  diff.stp.locs.n <- paste("diff",names(lat.std.list.raster)[i],sep = ".")
  diff.stp.locs[[diff.stp.locs.n]] <- diff.stp.locs.calc
  png(paste("./figures/lat_analysis/lat.std.stippling",
            names(lat.std.list.raster)[i],".png",
            sep = "."),
      width = 50,height = 25,units = "cm",res = 500)
  plot(lat.std.list.raster[[i]],col=rev(brewer.pal(10,"RdYlBu")),
       legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                          side=4, font=1, line=2.3),
       main="Butler Sla Mean (stippling)\n")
  points(stp.locs.75[[i]],pch=18,cex=0.8)
  points(diff.stp.locs[[i]],pch=23,cex=0.9,col="darkgreen",bg="green")
  dev.off()
}


#### LATITUDE STATS ####
  # sla mean
for (i in 1:length(lat.list)){
  colnames(lat.list[[i]]) <- sub("specific.leaf.area.x","cardamom",
                                 colnames(lat.list[[i]]))
  colnames(lat.list[[i]]) <- sub("specific.leaf.area.y","butler",
                                    colnames(lat.list[[i]]))
}
lat.list <-   lapply(lat.list, na.omit)
lat.lm <- lm.f(lat.list)
lat.p <- lapply(lat.lm, lmp)
lat.stats <- stats.f(lat.list)

# sla stdev
for (i in 1:length(lat.std.list)){
  colnames(lat.std.list[[i]]) <- sub("sla_std","cardamom",
                                       colnames(lat.std.list[[i]]))
  colnames(lat.std.list[[i]]) <- sub("specific.leaf.area","butler",
                                       colnames(lat.std.list[[i]]))
}

lat.std.list <- lapply(lat.std.list, na.omit)
lat.std.lm <- lm.f(lat.std.list)
lat.std.p <- lapply(lat.std.lm, lmp)
lat.std.stats <- stats.f(lat.std.list)


#### HEATSCATTER LATS SLA MEAN ####
png("./figures/lat_analysis/heatsc_latitude_sla.png", width = 40, height = 30, 
    units = "cm", res = 500)
par(mfcol = c(2,2))
# TROPICS 
trp_c_sla_n <- lat.list$Tropics[["cardamom"]]
trp_b_sla_n <- lat.list$Tropics[["butler"]]
(heatsc_sla_trp <- heatscatter(trp_c_sla_n, trp_b_sla_n, pch = 19, 
                                cexplot = 0.5, colpal="spectral", 
                                #disco() for all color options / could set to colorblind
                                add.contour=TRUE, main = "Tropics SLA Mean",
                                xlab="", 
                                ylab="Butler"))

# SUBTROPICS
sbtrp_c_sla_n <- lat.list$Subtropics[["cardamom"]]
sbtrp_b_sla_n <- lat.list$Subtropics[["butler"]]

(heatsc_sla_sbtrp <- heatscatter(sbtrp_c_sla_n, sbtrp_b_sla_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Subtropics",
                                  xlab="Cardamom", 
                                  ylab="Butler"))

# TEMPERATE
tmp_sla_c_n <- lat.list$Temperate[["cardamom"]]
tmp_sla_b_n <- lat.list$Temperate[["butler"]]
(heatsc_sla_tmp <- heatscatter(tmp_sla_c_n, tmp_sla_b_n, pch = 19, 
                               cexplot = 0.5, colpal="spectral", 
                               #disco() for all color options / could set to colorblind
                               add.contour=TRUE, main = "Temperate",
                               xlab="", 
                               ylab=""))

# N POLE
pl_sla_c_n <- lat.list$Pole[["cardamom"]]
pl_sla_b_n <- lat.list$Pole[["butler"]]
(heatsc_sla_pl <- heatscatter(pl_sla_c_n, pl_sla_b_n, pch = 19, 
                               cexplot = 0.5, colpal="spectral", 
                               add.contour=TRUE, main = "N pole",
                               xlab="Cardamom", 
                               ylab=""))

# dev.off----
dev.off()


#### HEATSCATTER LATS SLA STDEV ####
png("./figures/lat_analysis/heatsc_latitude_std.png", width = 40, height = 30, 
    units = "cm", res = 500)
par(mfcol = c(2,2))
# TROPICS 
trp_c_std_n <- lat.std.list$Tropics.std[["cardamom_std"]]
trp_b_std_n <- lat.std.list$Tropics.std[["butler_std"]]
(heatsc_slastd_trp <- heatscatter(trp_c_std_n, trp_b_std_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Tropics SLA StDev",
                                  xlab="", 
                                  ylab="Butler"))
# SUBTROPICS
sbtrp_c_slastd_n <- lat.std.list$Subtropics.std[["cardamom_std"]]
sbtrp_b_slastd_n <- lat.std.list$Subtropics.std[["butler_std"]]

(heatsc_slastd_sbtrp <- heatscatter(sbtrp_c_slastd_n, sbtrp_b_slastd_n, pch = 19, 
                                    cexplot = 0.5, colpal="spectral", 
                                    #disco() for all color options / could set to colorblind
                                    add.contour=TRUE, main = "Subtropics",
                                    xlab="Cardamom", 
                                    ylab="Butler"))
# TEMPERATE
tmp_slastd_c_n <- lat.std.list$Temperate.std[["cardamom_std"]]
tmp_slastd_b_n <- lat.std.list$Temperate.std[["butler_std"]]
(heatsc_slastd_tmp <- heatscatter(tmp_slastd_c_n, tmp_slastd_b_n, pch = 19, 
                                  cexplot = 0.5, colpal="spectral", 
                                  #disco() for all color options / could set to colorblind
                                  add.contour=TRUE, main = "Temperate",
                                  xlab="", 
                                  ylab=""))

# N POLE
pl_slastd_c_n <- lat.std.list$Pole.std[["cardamom_std"]]
pl_slastd_b_n <- lat.std.list$Pole.std[["butler_std"]]
(heatsc_slastd_pl <- heatscatter(pl_slastd_c_n, pl_slastd_b_n, pch = 19, 
                                 cexplot = 0.5, colpal="spectral", 
                                 #disco() for all color options / could set to colorblind
                                 add.contour=TRUE, main = "N Pole",
                                 xlab="Cardamom", 
                                 ylab=""))

dev.off()

#### SPLITTING BY BIOME #####

# OPEN DATASET: ECOREGIONS17 ####
ecoregions17 <- st_read("./DATA/Ecoregions2017/Ecoregions2017.shp")
st_crs(ecoregions17)
st_bbox(ecoregions17)
#      xmin       ymin       xmax       ymax 
#-179.99999  -89.89197  180.00000   83.62313 

#unique(ecoregions17[[4]])
#ecoregions17_geom <- st_geometry(ecoregions17)
#ecoregions17_geom[[1]]

# splitting ecosystem17 dataset by biome, list of dataframes
filtered.biome <- split(ecoregions17, f=ecoregions17$BIOME_NAME)
filtered.biome$`N/A` <- NULL # removing ecoregion "rocks and ice", equivalent of N/A biome in dataset 
                              # check https://ecoregions2017.appspot.com/ for confirmation
                              # confirmed

# creating mask of cardamom and butler based on biome ----
  # sla mean
m.cardamom.biome <- list()
for (i in 1:length(filtered.biome)){
    mask.biome <- raster::mask(cardamom_sla, filtered.biome[[i]])
    if (!is.infinite(mask.biome@data@min)&!is.infinite(mask.biome@data@max)){
      name.biome <- names(filtered.biome)[i]
      m.cardamom.biome[[name.biome]] <- mask.biome
    }
    else
      NULL
}

m.butler.biome <- list()
for (i in 1:length(filtered.biome)){
    mask.biome <- raster::mask(butler_sla, filtered.biome[[i]])
    if (!is.infinite(mask.biome@data@min)&!is.infinite(mask.biome@data@max)){
      name.biome <- names(filtered.biome)[i]
      m.butler.biome[[name.biome]] <- mask.biome
    }
    else
      NULL
}

  # sla stdev 
m.cardamom.std.biome <- list()
for (i in 1:length(filtered.biome)){
    mask.biome <- raster::mask(cardamom_sla_std, filtered.biome[[i]])
    if (!is.infinite(mask.biome@data@min)&!is.infinite(mask.biome@data@max)){
      name.biome <- names(filtered.biome)[i]
      m.cardamom.std.biome[[name.biome]] <- mask.biome
    }
    else
      NULL
}

m.butler.std.biome <- list()
for (i in 1:length(filtered.biome)){
    mask.biome <- raster::mask(butler_sla_std, filtered.biome[[i]])
    if (!is.infinite(mask.biome@data@min)&!is.infinite(mask.biome@data@max)){
      name.biome <- names(filtered.biome)[i]
      m.butler.std.biome[[name.biome]] <- mask.biome
    }
    else
      NULL
}

#### VISUAL AND STAT ANALYSIS BY BIOME ####

  # sla mean
biome.cardamom.df <- mask.to.df(m.cardamom.biome) 
biome.cardamom.df <- lapply(biome.cardamom.df, function(x){
  rename(x, "cardamom"=specific.leaf.area)
})
biome.butler.df <- mask.to.df(m.butler.biome)
biome.butler.df <- lapply(biome.butler.df, function(x){
  rename(x, "butler"=specific.leaf.area)
})

  # sla stdev
biome.cardamom.df.std <- mask.to.df(m.cardamom.std.biome)
biome.cardamom.df.std <- lapply(biome.cardamom.df.std, function(x){
  rename(x, "cardamom_std"=sla_std)
})
biome.butler.df.std <- mask.to.df(m.butler.std.biome)
biome.butler.df.std <- lapply(biome.butler.df.std, function(x){
  rename(x, "butler_std"=specific.leaf.area)
})

  # sla mean
j.biome.sla <- join.f(biome.cardamom.df,biome.butler.df)
names(j.biome.sla) <- sub(" ","",names(j.biome.sla))
  # sla stdev 
j.biome.slastd <- join.f(biome.cardamom.df.std,biome.butler.df.std)
names(j.biome.slastd) <- sub(" ","",names(j.biome.slastd))

#### PERCENTAGE OVERLAP DENSITY PLOTS ####
  # sla mean ----
biome.list.n <- list()

for (i in 1:length(j.biome.sla)){
  j.biome.sla[[i]] <- j.biome.sla[[i]]%>%
    filter(cardamom!=0,butler!=0) 
  list <- list(cardamom=j.biome.sla[[i]][["cardamom"]],
               butler=j.biome.sla[[i]][["butler"]])
  name <- names(j.biome.sla)[i]
  biome.list.n[[name]] <- list
  for (j in 1:length(biome.list.n)){
  png(paste0("./figures/", i, ".png"), 
      width = 40,height = 20, units = "cm", res = 500)
  my.overlap.sla(biome.list.n[[j]], plot = TRUE)
  dev.off()
  }
} # something wrong with naming of the plots... still havent figured out
# also for stdev - something to do witht the list? hm

# sla stdev ----
biome.std.list.n <- list()

for (i in 1:length(j.biome.slastd)){
  j.biome.slastd[[i]] <- j.biome.slastd[[i]]%>%
    filter(cardamom_std!=0,butler_std!=0) 
  list <- list(cardamom=j.biome.slastd[[i]][["cardamom_std"]],
               butler=j.biome.slastd[[i]][["butler_std"]])
  name <- paste(names(j.biome.slastd)[i],"n",sep = "_")
  biome.std.list.n[[name]] <- list
  for (j in 1:length(biome.std.list.n)){
    plot_name <- paste("./figures/biome_std",i,
                       "png",sep=".")
    png(plot_name, width = 40, height = 20, units = "cm", res = 400)
    my.overlap.std(biome.std.list.n[[j]], plot = TRUE)
    dev.off()
  }
}

#### STIPPLING BY BIOME ####
## creation of function to do the automatisation of stippling by biome
biome.stp.25 <- list()
biome.stp.75 <-list()
biome.stp.95 <- list()
stp.locs.75 <- list()
stp.locs.95<- list()
diff.stp.locs <- list()
for (i in 1:length(m.butler.biome)){
    stp.25 <- raster::mask(cardamom_25th, m.butler.biome[[i]])
    stp.75 <- raster::mask(cardamom_75th, m.butler.biome[[i]])
    stp.95 <- raster::mask(cardamom_95th, m.butler.biome[[i]])
    stp.name25 <- paste("stp25", names(m.butler.biome)[i],sep = ".")
    stp.name75 <- paste("stp75", names(m.butler.biome)[i],sep = ".")
    stp.name95 <- paste("stp95",names(m.butler.biome)[i],sep=".")
    biome.stp.25[[stp.name25]] <- stp.25
    biome.stp.75[[stp.name75]] <- stp.75
   biome.stp.95[[stp.name95]] <- stp.95
   stp.locs.75.calc <- (m.butler.biome[[i]]>=biome.stp.25[[i]])*2
    (m.butler.biome[[i]]<= biome.stp.75[[i]])
  stp.locs.95.calc <- (m.butler.biome[[i]]>=biome.stp.25[[i]])*
    (m.butler.biome[[i]]<= biome.stp.95[[i]])
  stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
  stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
  stp.locs.75.p <- as.data.frame(stp.locs.75.p[stp.locs.75.p[, "layer"] == 1,])
  stp.locs.95.p <- as.data.frame(stp.locs.95.p[stp.locs.95.p[, "layer"] == 1,])
  stp.locs.75n <- paste("stp.locs",names(m.butler.biome)[i],"75",sep = ".")
  stp.locs.95n <- paste("stp.locs",names(m.butler.biome)[i],"95",sep = ".")
  stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
  stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
  diff.stp.locs.calc <- as.matrix(anti_join(stp.locs.95[[i]],
                                            stp.locs.75[[i]]))
  diff.stp.locs.n <- paste("diff",names(m.butler.biome)[i],sep = ".")
  diff.stp.locs[[diff.stp.locs.n]] <- diff.stp.locs.calc
  png(paste("./figures/stippling",i,"png",sep = "."),
      width = 50,height = 25,units = "cm",res = 500)
  plot(m.butler.biome[[i]],asp=NA,col=rev(brewer.pal(10,"RdYlBu")),
       legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                          side=4, font=1, line=2.3),
       main="Butler Sla Mean (stippling)\n")
  points(stp.locs.75[[i]],pch=18,cex=0.8)
}

# stippling biomes for sla stdev
biome.stp.25 <- list()
biome.stp.75 <-list()
biome.stp.95 <- list()
stp.locs.75 <- list()
stp.locs.95<- list()
diff.stp.locs <- list()
for (i in 1:length(m.butler.std.biome)){
    stp.25 <- raster::mask(cardamom_25th, m.butler.std.biome[[i]])
    stp.75 <- raster::mask(cardamom_75th, m.butler.std.biome[[i]])
    stp.95 <- raster::mask(cardamom_95th, m.butler.std.biome[[i]])
    stp.name25 <- paste("stp25", names(m.butler.std.biome)[i],sep = ".")
    stp.name75 <- paste("stp75", names(m.butler.std.biome)[i],sep = ".")
    stp.name95 <- paste("stp95",names(m.butler.std.biome)[i],sep=".")
    biome.stp.25[[stp.name25]] <- stp.25
    biome.stp.75[[stp.name75]] <- stp.75
    biome.stp.95[[stp.name95]] <- stp.95
    stp.locs.75.calc <- (m.butler.std.biome[[i]]>=biome.stp.25[[i]])*
      (m.butler.std.biome[[i]]<= biome.stp.75[[i]])
    stp.locs.95.calc <- (m.butler.std.biome[[i]]>=biome.stp.25[[i]])*
      (m.butler.std.biome[[i]]<= biome.stp.95[[i]])
    stp.locs.75.p <- rasterToPoints(stp.locs.75.calc)
    stp.locs.95.p <- rasterToPoints(stp.locs.95.calc)
    stp.locs.75.p <- stp.locs.75.p[stp.locs.75.p[, "layer"] == 1,]
    stp.locs.95.p <- stp.locs.95.p[stp.locs.95.p[, "layer"] == 1,]
    stp.locs.75n <- paste("stp.locs",names(m.butler.std.biome)[i],"75",
                          sep = ".")
    stp.locs.95n <- paste("stp.locs",names(m.butler.std.biome)[i],"95",
                          sep = ".")
    stp.locs.75[[stp.locs.75n]]<- stp.locs.75.p
    stp.locs.95[[stp.locs.95n]]<- stp.locs.95.p
    png(paste("./figures/biome_analysis/stippling.std",i,"png",
              sep = "."),
        width = 50,height = 25,units = "cm",res = 500)
    plot(m.butler.std.biome[[i]],asp=NA,col=rev(brewer.pal(10,"RdYlBu")),
         legend.args = list(text="\n\nSLA Mean (m2.kg-1)", 
                            side=4, font=1, line=2.3),
         main="Butler Sla Mean (stippling)\n")
    points(stp.locs.75[[i]],pch=18,cex=0.8)
    dev.off()
}

# HEATSCATTER SLA MEAN MAJ BIOMES ----
j.biome.sla.n <- lapply(j.biome.sla, function(x){
  x[c("cardamom","butler")]
})
j.biome.sla.n <- unlist(j.biome.sla.n, recursive = F)

png("./figures/biome_analysis/panel_heatsc_biome_sla.png", width = 50, height = 30,
    units = "cm", res = 500)
par(mfcol=c(2,4))
# 1) taiga ----
(heatsc_taiga_sla <- heatscatter(j.biome.sla.n$`df.BorealForests/Taiga.cardamom`, 
                                 j.biome.sla.n$`df.BorealForests/Taiga.butler`, 
                                 pch = 19, cexplot = 0.5, colpal="spectral", 
                                 add.contour=TRUE, main = "Taiga",
                                 xlab=" ", 
                                 ylab="Butler"))
  
# 2) tundra ----
(heatsc_taiga_sla <- heatscatter(j.biome.sla.n$df.Tundra.cardamom, 
                                 j.biome.sla.n$df.Tundra.butler, 
                                 pch = 19, cexplot = 0.5, colpal="spectral", 
                                 add.contour=TRUE, main = "Tundra",
                                 xlab="Cardamom", 
                                 ylab="Butler"))

# 3) temp conif forest ----
(heatsc_temp_con <- heatscatter(j.biome.sla.n$`df.TemperateConifer Forests.cardamom`, 
                                j.biome.sla.n$`df.TemperateConifer Forests.butler`, 
                                 pch = 19, cexplot = 0.5, colpal="spectral", 
                                 add.contour=TRUE, main = "Temperate coniferous",
                                 xlab=" ", 
                                 ylab=" "))

# 4) temp broad mix forest ----
(heatsc_temp_con <- heatscatter(j.biome.sla.n$`df.TemperateBroadleaf & Mixed Forests.cardamom`, 
                                j.biome.sla.n$`df.TemperateBroadleaf & Mixed Forests.butler`, 
                                pch = 19, cexplot = 0.5, colpal="spectral", 
                                add.contour=TRUE, main = "Temperate broad/mixed",
                                xlab="Cardamom", 
                                ylab=" "))

# 5) tropical and subtropical dry broadleaf ----
(heatsc_trpsbtrp_d_broad <- heatscatter(j.biome.sla.n$`df.Tropical& Subtropical Dry Broadleaf Forests.cardamom`, 
                                        j.biome.sla.n$`df.Tropical& Subtropical Dry Broadleaf Forests.butler`, 
                                pch = 19, cexplot = 0.5, colpal="spectral", 
                                add.contour=TRUE, 
                                main = "Tropical/subtropical dry broadleaf",
                                xlab=" ", 
                                ylab=" "))

# 6) tropical and subtropical conif forest ----
(heatsc_trpsbtrp_con <- heatscatter(j.biome.sla.n$`df.Tropical& Subtropical Coniferous Forests.cardamom`, 
                                    j.biome.sla.n$`df.Tropical& Subtropical Coniferous Forests.butler`, 
                                pch = 19, cexplot = 0.5, colpal="spectral", 
                                add.contour=TRUE, main = "Tropical/subtropical coniferous",
                                xlab="Cardamom", 
                                ylab=" "))

# 7) tropical subtropical moist broadleaf ----
(heatsc_trpsbtrp_m_broad <- heatscatter(j.biome.sla.n$`df.Tropical& Subtropical Moist Broadleaf Forests.cardamom`, 
                                        j.biome.sla.n$`df.Tropical& Subtropical Moist Broadleaf Forests.butler`, 
                                    pch = 19, cexplot = 0.5, colpal="spectral", 
                                    add.contour=TRUE, 
                                    main = "Tropical/subtropical moist broadleaf",
                                    xlab=" ", 
                                    ylab=" "))

# 8) mediterranean forests, woodlands, scrub ----
(heatsc_med_f_w_scr <- heatscatter(j.biome.sla.n$`df.MediterraneanForests, Woodlands & Scrub.cardamom`, 
                                   j.biome.sla.n$`df.MediterraneanForests, Woodlands & Scrub.butler`, 
                                        pch = 19, cexplot = 0.5, colpal="spectral", 
                                        add.contour=TRUE, 
                                   main = "Mediterranean\nwoodland and scrubland",
                                        xlab="Cardamom", 
                                        ylab=" "))

# dev.off ----
dev.off()



# HEATSCATTER SLA STDEV MAJ BIOMES ----
j.biome.std.n <- lapply(j.biome.slastd, function(x){
  x[c("cardamom_std","butler_std")]
})
j.biome.std.n <- unlist(j.biome.std.n, recursive = F)

png("./figures/biome_analysis/panel_heatsc_slastd_biomes.png", width = 50,
    height = 30, units = "cm", res = 500)
par(mfcol = c(2,4))
# 1) Taiga ----
(heatsc_taiga_slastd <- heatscatter(j.biome.std.n$`df.BorealForests/Taiga.cardamom_std`, 
                                    j.biome.std.n$`df.BorealForests/Taiga.butler_std`, 
                                    pch = 19, cexplot = 0.5, colpal="spectral", 
                                    add.contour=TRUE, main = "Taiga",
                                    xlab="", 
                                    ylab="Butler"))
# 2) Tundra ----
(heatsc_tundra_slastd <- heatscatter(j.biome.std.n$df.Tundra.cardamom_std, 
                                     j.biome.std.n$df.Tundra.butler_std, 
                                    pch = 19, cexplot = 0.5, colpal="spectral", 
                                    add.contour=TRUE, main = "Tundra",
                                    xlab="Cardamom", 
                                    ylab="Butler"))

# 3) temp conif forest ----
(heatsc_temp_con_slastd <- heatscatter(j.biome.std.n$`df.TemperateConifer Forests.cardamom_std`, 
                                       j.biome.std.n$`df.TemperateConifer Forests.butler_std`, 
                                     pch = 19, cexplot = 0.5, colpal="spectral", 
                                     add.contour=TRUE, main = "Temperate coniferous",
                                     xlab="", 
                                     ylab=""))

# 4) temp broad mix forest ----
(heatsc_temp_b_m_slastd <- heatscatter(j.biome.std.n$`df.TemperateBroadleaf & Mixed Forests.cardamom_std`, 
                                       j.biome.std.n$`df.TemperateBroadleaf & Mixed Forests.butler_std`, 
                                       pch = 19, cexplot = 0.5, 
                                       colpal="spectral", 
                                       add.contour=TRUE, 
                                       main = "Temperate broadleaf mixed",
                                       xlab="Cardamom", 
                                       ylab=""))

# 5) tropical and subtropical dry broadleaf ----
(heatsc_trpsbtrp_d_b_slastd <- heatscatter(j.biome.std.n$`df.Tropical& Subtropical Dry Broadleaf Forests.cardamom_std`, 
                                           j.biome.std.n$`df.Tropical& Subtropical Dry Broadleaf Forests.butler_std`, 
                                       pch = 19, cexplot = 0.5, 
                                       colpal="spectral", 
                                       add.contour=TRUE, 
                                       main = "Tropical/subtropical dry broadleaf",
                                       xlab="", 
                                       ylab=""))

# 6) tropical and subtropical conif forest ----
(heatsc_trpsbtrp_con_slastd <- heatscatter(j.biome.std.n$`df.Tropical& Subtropical Coniferous Forests.cardamom_std`, 
                                           j.biome.std.n$`df.Tropical& Subtropical Coniferous Forests.butler_std`, 
                                           pch = 19, cexplot = 0.5, 
                                           colpal="spectral", 
                                           add.contour=TRUE, 
                                           main = "Tropical/subtropical coniferous",
                                           xlab="Cardamom", 
                                           ylab=""))

# 7) tropical subtropical moist broadleaf ----
(heatsc_trpsbtrp_con_slastd <- heatscatter(j.biome.std.n$`df.Tropical& Subtropical Moist Broadleaf Forests.cardamom_std`, 
                                           j.biome.std.n$`df.Tropical& Subtropical Moist Broadleaf Forests.butler_std`, 
                                           pch = 19, cexplot = 0.5, 
                                           colpal="spectral", 
                                           add.contour=TRUE, 
                                           main = "Tropical/subtropical moist broadleaf",
                                           xlab="", 
                                           ylab=""))

# 8) mediterranean forests, woodlands, scrub ----
(heatsc_med_f_w_scr_slastd <- heatscatter(j.biome.std.n$`df.MediterraneanForests, Woodlands & Scrub.cardamom_std`, 
                                          j.biome.std.n$`df.MediterraneanForests, Woodlands & Scrub.butler_std`, 
                                           pch = 19, cexplot = 0.5, 
                                           colpal="spectral", 
                                           add.contour=TRUE, 
                                           main = "Mediterranean woodland and scrubland",
                                           xlab="Cardamom", 
                                           ylab=""))
# dev.off ----
dev.off()


# STATS BIOMES: r2, rmse, bias ----
  # for sla stdev
for (i in 1:length(j.biome.slastd)){
  colnames(j.biome.slastd[[i]]) <- sub("cardamom_std","cardamom",
                                       colnames(j.biome.slastd[[i]]))
  colnames(j.biome.slastd[[i]]) <- sub("butler_std","butler",
                                       colnames(j.biome.slastd[[i]]))
}

# stats
  # sla mean 
biome.sla.lm <- lm.f(j.biome.sla)
biome.sla.p <- lapply(biome.sla.lm, lmp)
biome.sla.stats <- stats.f(j.biome.sla)

# sla stdev
biome.std.lm <- lm.f(j.biome.slastd)
biome.std.p <- lapply(biome.std.lm, lmp)
biome.slastd.stats <- stats.f(j.biome.slastd)


#### stat analysis by biome in diff continents  ####
world <- getMap()
plot(world)
world <- clgeo_Clean(world)  ## Needed to fix up some non-closed polygons 
cont <- sapply(levels(world$continent),
         FUN = function(i) {
           ## Merge polygons within a continent
           poly <- gUnionCascaded(subset(world, continent==i))
           ## Give each polygon a unique ID
           poly <- spChFIDs(poly, i)
           ## Make SPDF from SpatialPolygons object
           SpatialPolygonsDataFrame(poly,
                                    data.frame(continent=i, row.names=i))
         },
         USE.NAMES=TRUE)

  # sla mean
mask.biome.bycont.c <- mask.biome.f(m.cardamom.biome, cont) # 52 elements
mask.biome.bycont.b <- mask.biome.f(m.butler.biome, cont) # 49 elements
  # butler seems to be missing three biomes when divided by continent, 
  # sla stdev 
mask.biome.std.bycont.c <- mask.biome.f(m.cardamom.std.biome, cont)
mask.biome.std.by.cont.b <- mask.biome.f(m.butler.std.biome, cont)

# turning masked biomes by continent to dataframes
  # sla mean
biome.bycont.c.df <- mask.to.df(mask.biome.bycont.c) # cardamom rasterLayer to df
biome.bycont.b.df <- mask.to.df(mask.biome.bycont.b) # butler rasterLayer to df
  # sla stdev
biome.std.bycont.c.df <- mask.to.df(mask.biome.std.bycont.c)
biome.std.bycont.b.df <- mask.to.df(mask.biome.std.by.cont.b)

  # understand what biomes dont match between datasets when split by continent
diff <- setdiff(names(biome.bycont.c.df),names(biome.bycont.b.df))
View(diff) 
diff.std <- setdiff(names(biome.std.bycont.c.df),names(biome.std.bycont.b.df))
View(diff.std)
# df.Mangroves Africa, df.Mangroves Australia,df.Mangroves South America
# for both sla mean and stdev

  # remove them from cardamom 
biome.bycont.c.df$`df.Mangroves Africa` <- NULL
biome.bycont.c.df$`df.Mangroves Australia` <- NULL
biome.bycont.c.df$`df.Mangroves South America`<-NULL

biome.std.bycont.c.df$`df.Mangroves Africa` <- NULL
biome.std.bycont.c.df$`df.Mangroves Australia` <- NULL
biome.std.bycont.c.df$`df.Mangroves South America` <- NULL

  # join dataframes under one list by biome*continent
biome.bycont.c.df <- lapply(biome.bycont.c.df, function(x){
  rename(x,"cardamom"=specific.leaf.area)
})
biome.bycont.b.df <- lapply(biome.bycont.b.df, function(x){
  rename(x, "butler"=specific.leaf.area)
})
j.biome.bycont <-join.f(biome.bycont.c.df,biome.bycont.b.df)

biome.std.bycont.c.df <- lapply(biome.std.bycont.c.df, function(x){
  rename(x, "cardamom"=sla_std)
})
biome.std.bycont.b.df <- lapply(biome.std.bycont.b.df, function(x){
  rename(x, "butler"=specific.leaf.area)
})
j.biome.std.bycont <- join.f(biome.std.bycont.c.df,biome.std.bycont.b.df)

#### carry out stats for biomes*continent ----
# sla mean
biome.bycont.sla.lm <- lm.f(j.biome.bycont)
biome.bycont.sla.p <- lapply(biome.bycont.sla.lm, lmp)
biome.bycont.stats.sla <- stats.f(j.biome.bycont) 
# sla stdev
biome.bycont.std.lm <- lm.f(j.biome.std.bycont)
biome.bycont.std.p <- lapply(biome.bycont.std.lm, lmp)
biome.bycont.stats.slastd <- stats.f(j.biome.std.bycont) 

#### Preparation of table outputs with statistical results ####

  # one table for global, latitude and biome (not by continent) divisions:
sla.stats.biglist <- list(Global=list(global_sla_stat),Lat=lat.stats,
                          Biome=biome.sla.stats)
slastd.stats.biglist <- list(Global.std = list(global_slastd_stat),
                             Lat.std = lat.std.stats,
                             Biome.std = biome.slastd.stats)

stats.table.sla <- lapply(sla.stats.biglist, function(x){
  lapply(x, function(y){
    unique(y[c("sla_r2","rmse_av","bias_av")])
  })
})

stats.table.std <- lapply(slastd.stats.biglist, function(x){
  lapply(x, function(y){
    unique(y[c("sla_r2","rmse_av","bias_av")])
  })
})

stats.table.sla <- lapply(stats.table.sla, function(x){
  lapply(rapply(x, enquote,how = "unlist"),eval)
})
stats.table.std <- lapply(stats.table.std, function(x){
  lapply(rapply(x,enquote,how = "unlist"),eval)
})

stats.table.sla <- lapply(rapply(stats.table.sla, enquote,
                                   how = "unlist"),eval)
stats.table.std <- lapply(rapply(stats.table.std, 
                                      enquote,how = "unlist"),eval)

stats.table.sla <- as.data.frame(do.call(rbind, stats.table.sla))
stats.table.sla<-setDT(stats.table.sla, keep.rownames = TRUE)[]
stats.table.std <- as.data.frame(do.call(rbind, stats.table.std))
stats.table.std<-setDT(stats.table.std, keep.rownames = TRUE)[]

# p values
p_val.list <- list(Global.sla =unique(global_sla_stat$p_val),
                   Global.std = unique(global_slastd_stat$p_val),
                   Lat.sla = lat.p, 
                   Lat.std = lat.std.p, Biome.sla = biome.sla.p,
                   Biome.std = biome.std.p,
                   Biome.bycont.sla = biome.bycont.sla.p,
                   Biome.bycont.std = biome.bycont.std.p)

p_vals <- lapply(rapply(p_val.list, enquote, how = "unlist"),eval)
p_vals <- as.data.frame(do.call(rbind, p_vals))
p_vals<-setDT(p_vals, keep.rownames = TRUE)[]
p_vals <- p_vals %>%
  mutate(rn=str_replace(rn,"Global","Global..")) %>%
  mutate(rn=str_replace(rn,"df.","")) %>%
  mutate(rn=str_replace(rn,"Biome.bycont","Biomebycont")) %>%
  separate(rn, into = c("Area","Index"),extra = "merge") %>%
  separate(Index,into = c("Stats","Index"),sep = "\\.") %>%
  dplyr::select(Area,Index,Stats,V1) %>%
  group_by(Stats) %>% 
  mutate(grouped_id = row_number()) %>%
  spread(key=Stats,value=V1) %>%
  dplyr::select(-grouped_id) %>%
  rename(sla.p=sla,std.p=std) 
p_vals.bycont <- p_vals %>%
  filter(Area=="Biomebycont")
p_vals <- p_vals %>%
  filter(Area!="Biomebycont") %>%
  arrange_at(1:length(.)) %>%
  arrange(match(Area, c("Global", "Lat", "Biome"))) %>%
  mutate(Index=replace_na(Index,"Global"))

stats.table.sla <- stats.table.sla %>%
  mutate(rn=str_replace(rn,"Global","Global.Global"))%>%
  mutate(rn=str_replace(rn,"df.","")) %>%
  mutate(rn=str_replace(rn,"biome.cont","biomebycont")) %>%
  separate(rn,into = c("Area","Index"),extra = "merge") %>%
  separate(Index,into = c("Index","Stats"),sep = "\\.") %>%
  spread(key = Stats, value = V1) %>%
  rename("Sla Mean (bias)"=bias_av, "Sla Mean (rmse)"=rmse_av,
         "Sla Mean (r2)"=sla_r2) 

stats.table.std <- stats.table.std %>%
  mutate(rn=str_replace(rn,"Global","Global.Global")) %>%
  mutate(rn = str_replace_all(rn,"std.","")) %>%
  mutate(rn=str_replace(rn,"df.","")) %>%
  mutate(rn=str_replace(rn,"biome.cont","biomebycont")) %>%
  separate(rn,into = c("Area","Index"),extra = "merge") %>%
  separate(Index,into = c("Index","Stats"),sep = "\\.") %>%
  spread(key = Stats, value = V1) %>%
  rename("Sla StDev (bias)"=bias_av,"Sla StDev (rmse)"=rmse_av,
         "Sla StDev (r2)"=sla_r2)

stats.table.out <- left_join(stats.table.sla,stats.table.std) %>%
  dplyr::select(Area,Index,`Sla Mean (r2)`,`Sla StDev (r2)`,
                `Sla Mean (rmse)`,`Sla StDev (rmse)`, 
                `Sla StDev (rmse)`, `Sla Mean (bias)`,
                `Sla StDev (bias)`) %>%
  arrange_at(1:length(.)) %>%
  arrange(match(Area, c("Global", "Lat", "Biome"))) %>%
  left_join(.,p_vals) %>%
  format(digits=2,nsmall = 2) 

(stats.out <- stats.table.out[2:8] %>% #no p-vals and first column
   # format_cells(2:5,3,"bold") %>%
  kable(digits = 4, "latex",booktabs = T) %>%
  kable_styling(latex_options = c("striped", "scale_down"),
    full_width = F,position = "center", font_size = 12) %>%
     add_header_above(c(" ","R2" = 2, "RMSE" = 2, "Bias" = 2), 
                     bold = T) %>%
  kableExtra::group_rows("Latitudinal gradient", 2,5) %>%
  kableExtra::group_rows("Biome",6,19) %>% # NEED TO TRY AND FIX THIS ONCE I PUT IT IN THE RMARKDOWN!!
  as_image(stats.table.out, 
           file= "./figures/table_stats-global-lat-biome.png", 
           width = 4))

   # one table for biomes divided by continent:
stats.table.sla.bycont <- lapply(biome.bycont.stats.sla, function(x){
    unique(x[c("sla_r2","rmse_av","bias_av")])
})

stats.table.std.bycont <- lapply(biome.bycont.stats.slastd, function(x){
    unique(x[c("sla_r2","rmse_av","bias_av")])
})

stats.table.sla.bycont <- lapply(stats.table.sla.bycont, function(x){
  lapply(rapply(x, enquote,how = "unlist"),eval)
})
stats.table.std.bycont <- lapply(stats.table.std.bycont, function(x){
  lapply(rapply(x,enquote,how = "unlist"),eval)
})

stats.table.sla.bycont <- lapply(rapply(stats.table.sla.bycont, enquote,
                                 how = "unlist"),eval)
stats.table.std.bycont <- lapply(rapply(stats.table.std.bycont, 
                                 enquote,how = "unlist"),eval)

stats.table.sla.bycont <- as.data.frame(do.call(rbind, stats.table.sla.bycont))
stats.table.sla.bycont<-setDT(stats.table.sla.bycont, keep.rownames = TRUE)[]
stats.table.std.bycont <- as.data.frame(do.call(rbind, stats.table.std.bycont))
stats.table.std.bycont<-setDT(stats.table.std.bycont, keep.rownames = TRUE)[]

stats.table.sla.bycont <- stats.table.sla.bycont %>%
  mutate(rn=str_replace(rn,"df.","")) %>%
  mutate(rn=str_replace(rn,"biome.cont","biomebycont")) %>%
  separate(rn,into = c("Biomes by continent","Stats"),sep = "\\.") %>%
  spread(key = Stats, value = V1) %>%
  rename("Sla Mean (bias)"=bias_av, "Sla Mean (rmse)"=rmse_av,
         "Sla Mean (r2)"=sla_r2) 

stats.table.std.bycont <- stats.table.std.bycont %>%
  mutate(rn = str_replace_all(rn,"std.","")) %>%
  mutate(rn=str_replace(rn,"df.","")) %>%
  mutate(rn=str_replace(rn,"biome.cont","biomebycont")) %>%
  separate(rn,into = c("Biomes by continent","Stats"),sep = "\\.") %>%
  spread(key = Stats, value = V1) %>%
  rename("Sla StDev (bias)"=bias_av,"Sla StDev (rmse)"=rmse_av,
         "Sla StDev (r2)"=sla_r2)

stats.table.bycont.out <- left_join(stats.table.sla.bycont,
                                    stats.table.std.bycont) %>%
  dplyr::select(`Biomes by continent`,`Sla Mean (r2)`,`Sla StDev (r2)`,
                `Sla Mean (rmse)`,`Sla StDev (rmse)`, `Sla StDev (rmse)`, 
                `Sla Mean (bias)`,`Sla StDev (bias)`) %>%
  filter_all(all_vars(!is.na(.))) %>% # removed three biomes which had NAs 
  # (not enough data points to do correlation = only one data point)
  arrange_at(1:length(.))

# table output for biomes split by continent
(stats.bycont.out <- stats.table.bycont.out %>%
    kable(digits = 7, "latex", booktabs = T) %>%
    kable_styling(latex_options = c("striped", "scale_down"),
                  full_width = F,position = "center", font_size = 12) %>%
    column_spec(1, italic = T, border_right = T,include_thead	=F) %>%
    row_spec(0,bold = T)%>%
    add_header_above(c(" ", "R2" = 2, "RMSE" = 2, "Bias" = 2), 
                     bold = T) %>%
    kableExtra::group_rows("Boreal", 1,2) %>%
    kableExtra::group_rows("Desert",3,7) %>%
    kableExtra::group_rows("Flooded",8,10) %>%
    kableExtra::group_rows("Mangroves",11,11) %>%
    kableExtra::group_rows("Mediterranean",12,16) %>%
    kableExtra::group_rows("Montane",17,20) %>%
    kableExtra::group_rows("Temperate",21,30) %>%
    kableExtra::group_rows("Tropical & Subtropical",31,44) %>%
    kableExtra::group_rows("Tundra",45,46) %>%
    as_image(stats.table.out, file= "./figures/table_stats-biome-by-cont.png", 
             width = 5, dpi = 500))

 
#### turning data back to raster for cooler visualisation ####

  # raster visualisation of RMSE values 

  # sla mean
r.rmse <- lapply(sla.stats.biglist, function(x){
  lapply(x, function(y){
    unique(y[c("x","y","rmse_row")])
  })
})

r.rmse <- unlist(r.rmse,recursive=FALSE) # to list of dataframes

spg.rmse <- r.rmse # calling it differently, to make a copy of it as im turning it to sp
#spg <- lapply(spg, coordinates)
spg.rmse <- lapply(spg.rmse, rasterFromXYZ)

  # raster visualisation of BIAS values 

r.bias <- lapply(sla.stats.biglist, function(x){
  lapply(x, function(y){
    unique(y[c("x","y","bias_row")])
  })
})

r.bias <- unlist(r.bias, recursive = F)

spg.bias <- r.bias
spg.bias <- lapply(spg.bias, rasterFromXYZ)

## plotting 
levelplot(spg.bias$`Biome.df.BorealForests/Taiga`,
          par.settings=RdBuTheme, margin =F,asp=NA)
levelplot(spg.rmse$`Biome.df.Boreal Forests/Taiga`,
          par.settings=RdBuTheme, margin=F, asp=NA)

r.stats.sla <- lapply(sla.stats.biglist, function(x){
  lapply(x, function(y){
    unique(y[c("x","y","rmse_row","bias_row")])
  })
})

r.stats.sla <- unlist(r.stats.sla, recursive = F)

spg.stats.sla <- r.stats.sla
spg.stats.sla <- lapply(spg.stats.sla, rasterFromXYZ)

png("./figures/sla.stats.global.png", width = 40, height = 30, units = "cm",
    res = 500)
levelplot(spg.stats.sla$Global, 
          par.settings=RdBuTheme(region=rev(brewer.pal(9,'RdBu'))),
          main="Results for RMSE and Bias for each 1°x1° grid value,
          between Cardamom and Butler SLA mean values", 
          colorkey = list(title=expression(m^2/kg), vjust=2))
dev.off()

  # sla stdev

r.stats.std <- lapply(slastd.stats.biglist, function(x){
  lapply(x, function(y){
    unique(y[c("x","y","rmse_row","bias_row")])
  })
})

r.stats.std <- unlist(r.stats.std, recursive = F)

spg.stats.std <- r.stats.std
spg.stats.std <- lapply(spg.stats.std, rasterFromXYZ)

# plotting standard dev for both rmse and bias - global
png("./figures/std.stats.global.png", width = 40, height = 30, units = "cm",
    res = 500)
levelplot(spg.stats.std$Global, 
          par.settings=RdBuTheme(region=rev(brewer.pal(9,'RdBu'))),
          main="Results for RMSE and Bias for each 1°x1° grid value,
          between Cardamom and Butler SLA StDev values", 
          colorkey = list(title=expression(m^2/kg), vjust=2))
dev.off()

# rasterising world map to overlay other data 
r <- raster(ncols=360, nrows=180)
r.world <- rasterize(world,r,progress="text")
r.world[NA] <- NULL

#levelplot(spg.bias$Global,par.settings=RdBuTheme,margins=T)
par(mar = c(0.7, 2.5, 0.5, 2.5),mfcol=c(2,1))
plot(r.world,col="black",asp=NA,legend=F)
plot(spg.bias$Global,add=T,asp=NA,col=rev(brewer.pal(10,"RdBu")))

plot(r.world,col="black",asp=NA,legend=F)
plot(spg.rmse$Global,add=T,asp=NA,col=rev(brewer.pal(10,"RdBu")))



