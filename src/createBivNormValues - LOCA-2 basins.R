library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(mvtnorm)
library(png)
library(gridExtra)
library(reshape2)
library(jocre)
library(ellipse)


# settings 
setwd("C:\\Users\\warnold\\Local\\repos\\callite-decision-scaling\\gcm\\src")
dir_GCM <- file.path("../data/loca2-projection-summary-worksheets-indiv//")

# choose region
my_region <- "LakeMillerton"

# plotting probability
base_center<-2006
prob_interval_count <- 100
prob_levels <- c(0.68,.95)
prob_plot_title <- paste0("Projected Range of Likely Climate Changes by ", climate_period+base_center)


# set increments
temp_increment <- 0.25
precip_increment <- 5
Precip <- as.numeric(seq(-25, 25, by=precip_increment))
Temp <- seq(0, 5, by=temp_increment) 

# load files
gcm_summary_worksheets <- list.files(dir_GCM)
gcm_families <- read.csv(families_GCM_list)
gcm_families <- gcm_families %>% mutate(GCM = toupper(GCM)) # uppercase GCM names


#############################################################################################################################
##################################### LOAD AND PROCESS DATA #################################################################
#############################################################################################################################

# Load historical
P_hist <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[1]), sheet= "pr_hist", startRow= 1, cols=c(6:17)) # units are mm
T_hist <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[1]), sheet= "tas_hist", startRow= 1, cols=c(6:17))
P_hist$Total <- rowSums(P_hist)
T_hist$Total <- rowMeans(T_hist)
# read in changes in P and T for each year (30 year climate windows)
i<-1
GCM_data <- data.frame()
for (i in 1:length(gcm_summary_worksheets)) {
  period <- sub(".xlsx", "", strsplit(gcm_summary_worksheets[i], "_")[[1]][[4]])
  
  GCM <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(1))
  Realization <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(4))
  SSP <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(5))

  RegionID <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(2))
  RegionName <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(3))

  PRdata <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(6:17))
  TASdata <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "tas_fut", startRow= 1, cols=c(6:17))
  PRdata$Total <- rowSums(PRdata)
  TASdata$Total <- rowMeans(TASdata)
  PRdata$D_pr <- ((PRdata$Total-P_hist$Total)/P_hist$Total)*100
  TASdata$D_tas <- TASdata$Total-T_hist$Total
  deltas <- as.data.frame(cbind(PRdata[, "D_pr"], TASdata[, "D_tas"]))
  colnames(deltas) <- c("D_pr","D_tas")
  deltas$Model <- GCM$GCM
  deltas$Variant <- Realization$Realization
  deltas$SSP <- SSP$SSP
  deltas$Year <- period
  deltas$period <- base_center+i
  deltas$RegionID <- RegionID$Region.ID
  deltas$RegionName <- RegionName$Region
  GCM_data <- rbind(GCM_data, deltas)
}
GCM_data$Model = toupper(GCM_data$Model)

GCM_data <- filter(GCM_data, RegionName==my_region)

#############################################################################################################################
##################################### PROBABILITIES #########################################################################
#############################################################################################################################

gcm_mean <- GCM_data[,c('period','D_tas','D_pr')] %>%
  group_by(period) %>%
  summarise_all(.,funs(mean))
colnames(gcm_mean)<-c('period','DT','DP')
sigs <- list()
hotelling <- list()
for(i in 1: length(unique(GCM_data$period))){
  period_data <- subset(GCM_data, period==base_center+i)
  sigs[[i]]=cov(cbind("DTsig"=period_data$"D_tas","DPsig"= period_data$"D_pr"))
}
stress_test_grid <- expand.grid(Temp,Precip)
norm_probs_normalize.dtdp <- list()

for(i in 1:length(unique(GCM_data$period))){
  norm_probs <- dmvnorm(stress_test_grid,mean=as.numeric(gcm_mean[i,2:3]),sigma=sigs[[i]])
  norm_probs_normalize.dtdp[[i]] <- cbind(stress_test_grid,norm_probs/sum(norm_probs))
  colnames(norm_probs_normalize.dtdp[[i]]) <- c("T_lev","P_lev","Biv_Norm_Prob")
  norm_probs_normalize.dtdp[[i]]$period<- base_center+i

}
# saveRDS(norm_probs_normalize.dtdp, file= "../processed/biv_norm_values_dt-dp.RDS")
# write.csv(norm_probs_normalize.dtdp,'../processed/biv_norm_vals_dt-dp.csv')

# load connection in environment using ../runbook/00-Run-Setup.Rmd
# df <- do.call(rbind.data.frame, norm_probs_normalize.dtdp)
# dbWriteTable(con,'ref_gcm_bivNormVal_dt_half_degree_dp_minus_20_32models',df,append=FALSE, row.names=FALSE)

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

climate_period <- 2

plot245 <- subset(GCM_data, SSP=="ssp245" & period==base_center+climate_period)
plot370 <- subset(GCM_data, SSP=="ssp370" & period==base_center+climate_period)
plot585 <- subset(GCM_data, SSP=="ssp585" & period==base_center+climate_period)

# dt dp
norm_probs_normalize.dtdp.period <- norm_probs_normalize.dtdp[[climate_period]]
norm_probs_normalize.dtdp.period <- norm_probs_normalize.dtdp.period[order(norm_probs_normalize.dtdp.period$Biv_Norm_Prob),]
norm_probs_normalize.dtdp.period$area <- cumsum(norm_probs_normalize.dtdp.period$Biv_Norm_Prob)

EvalData_stress_matrix.dtdp <- norm_probs_normalize.dtdp.period[c("T_lev", "P_lev", "area")] %>% spread(key=T_lev,value="area") 
rownames(EvalData_stress_matrix.dtdp) <- EvalData_stress_matrix.dtdp$P_lev
EvalData_stress_matrix.dtdp<- EvalData_stress_matrix.dtdp[,-1]
prob_plot_title <- paste0("Projected Range of Likely Climate Changes (CMIP6 LOCA-2) by ", 
  climate_period+base_center,
  "\n(relative to baseline 30-yr period ", base_center-14,"-",base_center+15,")",
  "\n",my_region)

png(paste0("../figures/biv-norm-prob-dt-dp.png"), width = 4500,height = 3800, res=600)
prob_plot <- my.filled.contour(
  x=Precip, y=Temp, z=as.matrix(EvalData_stress_matrix.dtdp),
  levels=mylevel, col=mycolors, plot.axes={
     points(y=plot245$D_tas, x=plot245$D_pr, pch=20, col='green');
     points(y=plot370$D_tas, x=plot370$D_pr, pch=20, col='orange');
     points(y=plot585$D_tas, x=plot585$D_pr, pch=20, col='red');
     contour(x=Precip, y=Temp, z=as.matrix(1-EvalData_stress_matrix.dtdp), 
             add=TRUE, lwd=2.5, labcex=1,
             levels=prob_levels, col='dark blue',
             plot.axes= contour(Precip, Temp, z=1-EvalData_stress_matrix.dtdp));
     {axis(1);axis(2)}},
   main=prob_plot_title, 
   xlab="Change in Precipitation (%)", 
   ylab="Change in Temperature (C)",xlim=c(-25, 25)
  )
print(prob_plot)
dev.off()

