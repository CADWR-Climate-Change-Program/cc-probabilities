library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(mvtnorm)
library(png)
library(gridExtra)

setwd("C:\\Users\\warnold\\Local\\repos\\cc-probabilities\\src")

# settings 
filter_GCMs <- FALSE # filter by a GCM list (e.g. CCTAG-20)
temp_increment <- 0.1
precip_increment <- 0.1
base_period_center = 2006

# plotting probability
prob_interval_count <- 100
prob_levels <- c(0.68,.95)

# directories
dir_GCM <- file.path("../data/cmip5//")
families_GCM_list <- file.path("../data/gcm-families.csv")
filter_GCM_list <- file.path("../data/cctagGCM.csv")

# load files
gcm_summary_worksheets <- list.files(dir_GCM)
gcm_families <- read.csv(families_GCM_list)
gcm_families <- gcm_families %>% mutate(GCM = toupper(GCM)) # uppercase GCM names

# set increments
Precip <- as.numeric(seq(-25, 25, by=precip_increment))
Temp <- seq(0, 5, by=temp_increment) 


#############################################################################################################################
##################################### LOAD AND PROCESS DATA #################################################################
#############################################################################################################################

########### LOAD GCM PROCESSED FOR DWR VA - ONLY INCLUDES MEAN T and P CHANGES ##############################################
# Load historical
file_index = (base_period_center-1995)*2
P_hist45 <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[file_index-1]), sheet= "pr_fut", startRow= 1, cols=c(1:13)) # units are mm
T_hist45 <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[file_index-1]), sheet= "tas_fut", startRow= 1, cols=c(1:13))
P_hist85 <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[file_index]), sheet= "pr_fut", startRow= 1, cols=c(1:13)) # units are mm
T_hist85 <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[file_index]), sheet= "tas_fut", startRow= 1, cols=c(1:13))
P_hist45$Total <- rowSums(P_hist45[,-1])
T_hist45$Total <- rowMeans(T_hist45[,-1])
P_hist85$Total <- rowSums(P_hist85[,-1])
T_hist85$Total <- rowMeans(T_hist85[,-1])
# read in changes in P and T for each year (30 year climate windows)
GCM_data <- data.frame()
for (i in 1:length(gcm_summary_worksheets)) {
  name <- strsplit(gcm_summary_worksheets[i], "_")[[1]][5]
  RCP <- sub(".xlsx", "", strsplit(gcm_summary_worksheets[i], "_")[[1]][6])
  PRdata <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(1:13))
  TASdata <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "tas_fut", startRow= 1, cols=c(1:13))
  PRdata$Total <- rowSums(PRdata[, -1])
  TASdata$Total <- rowMeans(TASdata[, -1])
  T_histdata <- if(RCP=="rcp45"){T_hist45}else{T_hist85}
  P_histdata <- if(RCP=="rcp45"){P_hist45}else{P_hist85}
  PRdata$D_pr <- (PRdata$Total-P_histdata$Total)/P_histdata$Total*100
  TASdata$D_tas <- TASdata$Total-T_histdata$Total
  deltas <- full_join(PRdata[, c("GCM", "D_pr")], TASdata[, c("GCM", "D_tas")], by="GCM")
  deltas$RCP <- RCP
  deltas$Year <- name
  if (i%%2==0) {
    i_t = i/2
    deltas$period <- 1995+i_t
    
  } else {
    i_t = ((i-1)/2)+1
    deltas$period <- 1995+i_t
    
  }
  GCM_data <- rbind(GCM_data, deltas)
}
GCM_data <- GCM_data[complete.cases(GCM_data),]
GCM_data <- GCM_data %>% mutate(GCM = toupper(GCM))

write.csv(GCM_data,'../processed/cmip5_DT_DP.csv',row.names=FALSE)


if (filter_GCMs==TRUE) {
  filter_GCMs <- read.csv(filter_GCM_list, header= F)
  GCM_data <- subset(GCM_data, GCM %in% filter_GCMs$V1)
}

GCMs_dataframe_family <- GCM_data %>% 
  left_join(y=gcm_families, by="GCM") %>%
  group_by(GCM.Family, RCP, Year) %>%
  summarise(DT=mean(D_tas), DP=mean(D_pr))

write.csv(GCMs_dataframe_family,'../processed/cmip5_family_DT_DP.csv',row.names=FALSE)

#############################################################################################################################
##################################### PROBABILITIES #########################################################################
#############################################################################################################################

## DT DP ##
gcm_mean <- GCMs_dataframe_family[,c(3:5)] %>%
  group_by(Year) %>%
  summarise_all(.,list(mean=mean))
sigs <- list()
Years <- list()

for(i in 1: length(unique(GCMs_dataframe_family$Year)))
{
  Years[[i]] = subset(GCMs_dataframe_family, Year==unique(GCMs_dataframe_family$Year[i]))
  sigs[[i]]=cov(cbind("DTsig"= Years[[i]]["DT"],"DPsig"= Years[[i]]["DP"]))
}

stress_test_grid <- expand.grid(Temp,Precip)
norm_probs_normalize.dtdp <- list()
norm_probs_normalize_mean_high.dtdp <- list()
norm_probs_normalize_mean_low.dtdp <- list()

for(j in 1:length(unique(GCMs_dataframe_family$Year)))

{
  norm_probs <- dmvnorm(stress_test_grid,mean=as.numeric(gcm_mean[j,2:3]),sigma=sigs[[j]])
  norm_probs_normalize.dtdp[[j]] <- cbind(stress_test_grid,norm_probs/sum(norm_probs))
  colnames(norm_probs_normalize.dtdp[[j]]) <- c("T_lev","P_lev","Biv_Norm_Prob")
  norm_probs_normalize.dtdp[[j]]$period<- GCMs_dataframe_family$Year[[j]]
  
}

saveRDS(norm_probs_normalize.dtdp, file= "../processed/biv_norm_values_dt-dp_cmip5.RDS")
write.csv(norm_probs_normalize.dtdp,'../processed/biv_norm_vals_dt-dp_cmip5.csv')
write.csv(gcm_mean,'../processed/gcm_mean_cmip5_family.csv',row.names=FALSE)
write.csv(sigs,'../processed/gcm_sigs_cmip5_family.csv',row.names=TRUE)


#############################################################################################################################
##################################### PLOTTING  PROBABILITIES ###############################################################
#############################################################################################################################
mylevel <- seq(0,1, length.out = prob_interval_count)^5
mycolors <- c(colorRampPalette(c("white", "#c5e2f0","dark blue"))(prob_interval_count)[-prob_interval_count])


my.filled.contour <- function (x = seq(0, 1, length.out = nrow(z)), 
                                 y = seq(0, 1, length.out = ncol(z)), z, 
                                 xlim = range(x, finite = TRUE),
                                 ylim = range(y, finite = TRUE),
                                 zlim = range(z, finite = TRUE),
                                 levels = pretty(zlim, nlevels), nlevels = 20, 
                                 color.palette = cm.colors,
                                 col = color.palette(length(levels) - 1), 
                                 plot.title, plot.axes, key.title, key.axes, 
                                 asp = NA, xaxs = "i", yaxs = "i", las = 1,
                                 axes = TRUE, frame.plot = axes, ...){
    if (missing(z)) { if (!missing(x)) { if (is.list(x)) {
      z <- x$z
      y <- x$y
      x <- x$x }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      } } else stop("no 'z' matrix specified")}
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0,1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    if (!missing(key.title)) 
      key.title
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    .filled.contour(x, y, z, levels, col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }


climate_period <- 48
plotRCP45 <- subset(GCM_data, RCP=="rcp45" & period==1995+climate_period)
plotRCP85 <- subset(GCM_data, RCP=="rcp85" & period==1995+climate_period)

# dt dp
norm_probs_normalize.dtdp.period <- norm_probs_normalize.dtdp[[climate_period]]
norm_probs_normalize.dtdp.period <- norm_probs_normalize.dtdp.period[order(norm_probs_normalize.dtdp.period$Biv_Norm_Prob),]
norm_probs_normalize.dtdp.period$area <- cumsum(norm_probs_normalize.dtdp.period$Biv_Norm_Prob)

plotGCMs_dataframe_family <- subset(GCMs_dataframe_family, Year==norm_probs_normalize.dtdp[[climate_period]]$period[1])
plotRCP45_family <- subset(plotGCMs_dataframe_family, RCP=="rcp45")
plotRCP85_family <- subset(plotGCMs_dataframe_family, RCP=="rcp85")

EvalData_stress_matrix.dtdp <- norm_probs_normalize.dtdp.period[c("T_lev", "P_lev", "area")] %>% spread(key=T_lev,value="area") 
rownames(EvalData_stress_matrix.dtdp) <- EvalData_stress_matrix.dtdp$P_lev
EvalData_stress_matrix.dtdp<- EvalData_stress_matrix.dtdp[,-1]
prob_plot_title <- paste0("Projected Range of Likely Climate Changes (CMIP5) by ", 
  climate_period+1995,
  "\n(relative to baseline 30-yr period ", base_period_center-14,"-",base_period_center+15,")")

png(paste0("../figures/biv-norm-prob-dt-dp-cmip5-",climate_period+1995,".png"), width = 4500,height = 3800, res=600)
prob_plot <- my.filled.contour(x=Precip, y=Temp, z=as.matrix(EvalData_stress_matrix.dtdp),
   levels=mylevel, col=mycolors, plot.axes={
    #  points(y=plotRCP45$D_tas, x=plotRCP45$D_pr, pch=20, col='lightgreen');
     points(y=plotRCP45_family$DT, x=plotRCP45_family$DP, pch=10, col='green',cex=1);
    #  points(y=plotRCP85$D_tas, x=plotRCP85$D_pr, pch=20, col='tan');
     points(y=plotRCP85_family$DT, x=plotRCP85_family$DP, pch=10, col='red',cex=1);
     contour(x=Precip, y=Temp, z=as.matrix(1-EvalData_stress_matrix.dtdp), 
             add=TRUE, lwd=2.5, labcex=1,
             levels=prob_levels, col='dark blue',
             plot.axes= contour(Precip, Temp, z=1-EvalData_stress_matrix.dtdp));
     {axis(1);axis(2)}},
   main=prob_plot_title, 
   xlab="Change in Precipitation (%)", 
   ylab="Change in Temperature (C)",xlim=c(-25, 25))

print(prob_plot)
dev.off()

