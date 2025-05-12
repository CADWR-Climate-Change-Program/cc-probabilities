library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(mvtnorm)
library(png)
library(fitdistrplus)
library(svglite)
library(magick)
library(magrittr)

# settings 
setwd("C:\\Users\\warnold\\Local\\repos\\cc-probabilities\\src")

# %%
domain <- 'cv-flow-weighted'
gcm_filter <- c('HadGEM3-GC31-LL') 
filter_GCMs <- TRUE
filter_nmem_GCMs <- TRUE
temp_increment <- 0.05
precip_increment <- 0.5

# FOR PREVIOUS CALLITE MODELING OF -20 to +30 PR and 0 to +4 T
# temp_increment <- 0.5
# precip_increment <- 10

# plotting probability
base_center<-2006
prob_interval_count <- 100
prob_levels <- c(0.68,.95)

# load summary 30y monthly gcm files
dir_GCM <- file.path("../data/loca2/loca2-projection-summary-worksheets-flow//")
gcm_summary_worksheets <- list.files(dir_GCM)

# set increments
Precip <- as.numeric(seq(-25, 25, by=precip_increment))
Temp <- seq(0, 5, by=temp_increment) 

# FOR PREVIOUS CALLITE MODELING OF -20 to +30 PR and 0 to +4 T
# Precip <- as.numeric(seq(-20, 30, by=precip_increment))
# Temp <- seq(0, 4, by=temp_increment) 

#############################################################################################################################
##################################### LOAD AND PROCESS DATA #################################################################
#############################################################################################################################

lm_fits = read.csv(paste0('../processed/loca2_lmfits_1981_',domain,'.csv'))
lm_pr_hist <- (lm_fits$slope * (base_center-14)) + lm_fits$intercept

# Load historical 
P_hist <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[1]), sheet= "pr_hist", startRow= 1, cols=c(4:15)) # units are mm
T_hist <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[1]), sheet= "tas_hist", startRow= 1, cols=c(4:15))
P_hist$Total <- rowSums(P_hist)
T_hist$Total <- rowMeans(T_hist)
GCM <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[1]), sheet= "pr_hist", startRow= 1, cols=c(1))
Realization <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[1]), sheet= "pr_hist", startRow= 1, cols=c(2))
SSP <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[1]), sheet= "pr_hist", startRow= 1, cols=c(3))

# read in changes in P and T for each year (30 year climate windows)
GCM_data <- as.data.frame(cbind(P_hist$Total, T_hist$Total))
colnames(GCM_data) <- c("Pr_total","Tas_mean")
GCM_data$D_pr <- 0
GCM_data$D_pr_lm <- 0
GCM_data$D_tas <- 0
GCM_data$Model <- GCM$GCM
GCM_data$Variant <- Realization$Realization
GCM_data$SSP <- SSP$SSP
GCM_data$Year <- "Hist(1992-2021)"
GCM_data$period <- base_center

for (i in 1:length(gcm_summary_worksheets)) {
  period <- sub(".xlsx", "", strsplit(gcm_summary_worksheets[i], "_")[[1]][[4]])
  
  GCM <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(1))
  Realization <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(2))
  SSP <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(3))

  PRdata <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "pr_fut", startRow= 1, cols=c(4:15))
  TASdata <- read.xlsx(paste0(dir_GCM, gcm_summary_worksheets[i]), sheet= "tas_fut", startRow= 1, cols=c(4:15))
  PRdata$Total <- rowSums(PRdata)
  TASdata$Total <- rowMeans(TASdata)
  PRdata$D_pr <- ((PRdata$Total-P_hist$Total)/P_hist$Total)*100
  TASdata$D_tas <- TASdata$Total-T_hist$Total
  deltas <- as.data.frame(cbind(PRdata[, "D_pr"], TASdata[, "D_tas"]))
  colnames(deltas) <- c("D_pr","D_tas")

  lm_pr_fut <- ((lm_fits$slope * (base_center+i-14)) + lm_fits$intercept)
  lm_fits$D_pr_lm <- ((lm_pr_fut - lm_pr_hist) / lm_pr_hist)*100

  deltas$Pr_total <- PRdata$Total 
  deltas$Tas_mean <- TASdata$Total 
  deltas$Model <- GCM$GCM
  deltas$Variant <- Realization$Realization
  deltas$SSP <- SSP$SSP
  deltas$Year <- period
  deltas$period <- base_center+i

  deltas <- inner_join(deltas, lm_fits[,c('m','v','s','D_pr_lm')], by=c("Model"="m","Variant"="v","SSP"="s"))

  GCM_data <- rbind(GCM_data, deltas)
}

if (filter_GCMs ==TRUE) {
  GCM_data <- filter(GCM_data,!Model  %in% gcm_filter)
}

GCM_data$Model = toupper(GCM_data$Model)

write.csv(GCM_data,paste0('../processed/loca2_DT_DP_',domain,'.csv'),row.names=FALSE)

#############################################################################################################################
##################################### MODEL VARIANT COLLAPSE ################################################################
#############################################################################################################################

GCM_models_nmem <- subset(GCM_data[,c('Model','Variant','SSP','period')],period==2006) %>%
  group_by(Model,SSP) %>%
  dplyr::count() %>%
  subset(.,n>1) 

GCM_models_nmem$model_ssp <- paste0(GCM_models_nmem$Model,"_", GCM_models_nmem$SSP)
GCM_data$model_ssp <- paste0(GCM_data$Model,"_", GCM_data$SSP)

if (filter_nmem_GCMs ==TRUE) {
  GCM_data_filtered <- filter(GCM_data, model_ssp  %in% c(GCM_models_nmem$model_ssp))
  
  GCM_hist_means <- subset(GCM_data_filtered, period==base_center) %>% 
    group_by(Model,SSP) %>% #keep ```Variant```` in groupby if need non-variant avg
    summarise_all(.,funs(mean)) 
  GCM_hist_means$Year <- "Hist(1992-2021)"
  GCM_data <- GCM_hist_means

  for (i in 1:length(gcm_summary_worksheets)) {
    period <- sub(".xlsx", "", strsplit(gcm_summary_worksheets[i], "_")[[1]][[4]])
    GCM_fut_mean <- subset(GCM_data_filtered, period==base_center+i) %>% 
      group_by(Model,SSP) %>%  #keep ```Variant```` in groupby if need non-variant avg
      summarise_all(.,funs(mean)) 
    GCM_fut_mean$D_pr <- (GCM_fut_mean$Pr_total - GCM_hist_means$Pr_total)/GCM_hist_means$Pr_total*100
    GCM_fut_mean$D_tas <- GCM_fut_mean$Tas_mean - GCM_hist_means$Tas_mean
    GCM_fut_mean$Year <- period
    GCM_data <- rbind(GCM_data, GCM_fut_mean)
  }
}

write.csv(GCM_data,paste0('../processed/loca2_DT_DP_varavg_',domain,'.csv'),row.names=FALSE)


#############################################################################################################################
##################################### PROBABILITIES #########################################################################
#############################################################################################################################
pr_type <- 'D_pr_lm' #D_pr_lm: linear model precip | D_pr: 30y precip

gcm_mean <- GCM_data[,c('period','D_tas',pr_type)] %>%
  group_by(period) %>%
  summarise_all(.,funs(mean))
colnames(gcm_mean)<-c('period','DT','DP')
sigs <- list()
t_norm_fits <- list()
p_norm_fits <- list()

for(i in 1:(length(unique(GCM_data$period))-1)){
  period_data <- subset(GCM_data, period==base_center+i)
  sigs[[i]]=cov(cbind("DTsig"=period_data$"D_tas","DPsig"= period_data[pr_type]))
  t_norm_fits[[i]] <- fitdist(period_data$"D_tas", "norm")
  p_norm_fits[[i]] <- fitdist(as.numeric(unlist(period_data[pr_type])), "norm")
}

stress_test_grid <- expand.grid(Temp,Precip)
norm_probs_normalize.dtdp <- list()

for(i in 1:(length(unique(GCM_data$period))-1)){
  norm_probs <- dmvnorm(stress_test_grid,mean=as.numeric(gcm_mean[i,2:3]),sigma=sigs[[i]])
  norm_probs_normalize.dtdp[[i]] <- cbind(stress_test_grid,norm_probs/sum(norm_probs))
  colnames(norm_probs_normalize.dtdp[[i]]) <- c("T_lev","P_lev","Biv_Norm_Prob")
  norm_probs_normalize.dtdp[[i]]$period<- base_center+i

}

saveRDS(norm_probs_normalize.dtdp, file=paste0('../processed/biv_norm_values_dt-dp_loca2_',domain,'.RDS'))
write.csv(norm_probs_normalize.dtdp,paste0('../processed/biv_norm_vals_dt-dp_loca2_',domain,'.csv'))
write.csv(gcm_mean,paste0('../processed/gcm_mean_loca2_varavg_lm_',domain,'.csv'),row.names=FALSE)
write.csv(sigs,paste0('../processed/gcm_sigs_loca2_varavg_lm_',domain,'.csv'),row.names=TRUE)

# FOR PREVIOUS CALLITE MODELING OF -20 to +30 PR and 0 to +4 T
norm_probs_normalize.dtdp.table <- norm_probs_normalize.dtdp[[1]]
for(i in 2:(length(unique(GCM_data$period))-1)){
  norm_probs_normalize.dtdp.table <- rbind(norm_probs_normalize.dtdp.table,norm_probs_normalize.dtdp[[i]])
}
write.csv(norm_probs_normalize.dtdp.table,paste0('../processed/biv_norm_vals_dt-dp_for_DB_loca2_',domain,'.csv'))

#############################################################################################################################
##################################### PLOTTING  PROBABILITIES ###############################################################
#############################################################################################################################
mylevel <- seq(0,1, length.out = prob_interval_count)^5
mycolors <- c(colorRampPalette(c("white", "#d0d0d0","#373737"))(prob_interval_count)[-prob_interval_count])


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

# climate_period <- 14 #2020
# climate_period <- 34 #2040
# climate_period <- 37 #2043
# climate_period <- 44 #2050
# climate_period <- 54 #2060
climate_period <- 64 #2070
# climate_period <- 74 #2080

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
prob_plot_title <- paste0("Projected range of likely climate changes by ", climate_period+base_center,
  "\nrelative to the baseline 30-yr period ", base_center-14,"-",base_center+15,"")
# prob_plot_title <- ""

svg(paste0("../figures/loca2-biv-norm-prob-dt-dp-lm-varavg-",base_center+climate_period,"_",domain,".svg"), width = 5.5,height = 7)
# png(paste0("../figures/_loca2-biv-norm-prob-dt-dp-lm-varavg-",base_center+climate_period,"_",domain,".png"), width = 350,height = 550)
prob_plot <- my.filled.contour(
  x=Precip, y=Temp, z=as.matrix(EvalData_stress_matrix.dtdp),
  levels=mylevel, col=mycolors, plot.axes={
     points(y=as.numeric(unlist(plot245['D_tas'])), x=as.numeric(unlist(plot245[pr_type])), cex=1.5, pch=20, col='#86289c');
     points(y=as.numeric(unlist(plot370['D_tas'])), x=as.numeric(unlist(plot370[pr_type])), cex=1.5, pch=20, col='#efa823');
     points(y=as.numeric(unlist(plot585['D_tas'])), x=as.numeric(unlist(plot585[pr_type])), cex=1.5, pch=20, col='#e01a1a');
     contour(x=Precip, y=Temp, z=as.matrix(1-EvalData_stress_matrix.dtdp), 
             add=TRUE, lwd=2.5, labcex=1,vfont=c("sans serif", "bold"),
             levels=prob_levels, col='#373737',
             plot.axes= contour(Precip, Temp, z=1-EvalData_stress_matrix.dtdp));
     {axis(1);axis(2)}},
   main=prob_plot_title, 
   xlab="Change in Precipitation (%)", 
   ylab="Change in Temperature (C)",
   xlim=c(-25, 25),ylim=c(0,5)
  )
print(prob_plot)
dev.off()


################################################### normal dist for t & p ########################################################

# png(paste0("../figures//t_norm_plot.png"), width = 4500,height = 4500, res=600)
# t_norm_plot <- denscomp(t_norm_fits[[climate_period]], demp=TRUE, xlab='Change in Temperature (C)')
# print(t_norm_plot)
# dev.off()

# png(paste0("../figures//p_norm_plot.png"), width = 4500,height = 4500, res=600)
# p_norm_plot <- denscomp(p_norm_fits[[climate_period]], demp=TRUE, xlab='Change in Precipitation (%)')
# print(p_norm_plot)
# dev.off()

# quantile(t_norm_fits[[climate_period]], p=c(0.95, 0.8, 0.65, 0.5,0.05))
# quantile(p_norm_fits[[climate_period]], p=c(0.95, 0.8, 0.65, 0.5,0.05))
# GCM_models_nmem


################################################### create animation ###########################################################

for(climate_period in 1:64) {
# climate_period <- 64

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
  prob_plot_title <- paste0("Projected range of likely climate changes by ", climate_period+base_center,
    "\nrelative to the baseline 30-yr period ", base_center-14,"-",base_center+15,"")

  png(paste0("../figures/animation/loca2-biv-norm-prob-dt-dp-lm-varavg-",base_center+climate_period,"_",domain,".png"), width = 500,height = 700)
  prob_plot <- my.filled.contour(
    x=Precip, y=Temp, z=as.matrix(EvalData_stress_matrix.dtdp),
    levels=mylevel, col=mycolors, plot.axes={
      points(y=as.numeric(unlist(plot245['D_tas'])), x=as.numeric(unlist(plot245[pr_type])), cex=1.5, pch=20, col='#86289c');
      points(y=as.numeric(unlist(plot370['D_tas'])), x=as.numeric(unlist(plot370[pr_type])), cex=1.5, pch=20, col='#efa823');
      points(y=as.numeric(unlist(plot585['D_tas'])), x=as.numeric(unlist(plot585[pr_type])), cex=1.5, pch=20, col='#e01a1a');
      contour(x=Precip, y=Temp, z=as.matrix(1-EvalData_stress_matrix.dtdp), 
              add=TRUE, lwd=2.5, labcex=1,vfont=c("sans serif", "bold"),
              levels=prob_levels, col='#373737',
              plot.axes= contour(Precip, Temp, z=1-EvalData_stress_matrix.dtdp));
      {axis(1);axis(2)}},
    main=prob_plot_title, 
    xlab="Change in Precipitation (%)", 
    ylab="Change in Temperature (C)",xlim=c(-15, 15),ylim=c(0,5)
    )
  print(prob_plot)
  dev.off()
}

list.files(path='../figures/animation/', pattern = '*.png', full.names = TRUE) %>% 
        image_read() %>% # reads each path file
        image_join() %>% # joins image
        image_animate(fps=4) %>% # animates, can opt for number of loops
        image_write("../figures/loca2-biv-norm-prob-dt-dp-lm-varavg-animation.gif") # write to current dir


######################################### old format ##############################################################################
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
    # .filled.contour(x, y, z, levels, col)
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

# prob_plot_title <- ""
prob_plot_title <- paste0("Projected range of likely climate changes (CMIP6/LOCA-2 all variants)\nby ", climate_period+base_center,
  " relative to the baseline 30-yr period ", base_center-14,"-",base_center+15,"")
png(paste0("../figures/loca2-biv-norm-prob-dt-dp-lm-novaravg-",base_center+climate_period,"_",domain,".png"), width = 4500,height = 3800, res=600)
# png(paste0("../figures/_loca2-biv-norm-prob-dt-dp-lm-varavg-",base_center+climate_period,"_",domain,".png"), width = 350,height = 550)
prob_plot <- my.filled.contour(
  x=Precip, y=Temp, z=as.matrix(EvalData_stress_matrix.dtdp),
  levels=mylevel, col=mycolors, plot.axes={
     points(y=as.numeric(unlist(plot245['D_tas'])), x=as.numeric(unlist(plot245[pr_type])), cex=1.5, pch=20, col='green');
     points(y=as.numeric(unlist(plot370['D_tas'])), x=as.numeric(unlist(plot370[pr_type])), cex=1.5, pch=20, col='orange');
     points(y=as.numeric(unlist(plot585['D_tas'])), x=as.numeric(unlist(plot585[pr_type])), cex=1.5, pch=20, col='red');
    #  contour(x=Precip, y=Temp, z=as.matrix(1-EvalData_stress_matrix.dtdp), 
    #          add=TRUE, lwd=2.5, labcex=1,vfont=c("sans serif", "bold"),
    #          levels=prob_levels, col='#373737',
    #          plot.axes= contour(Precip, Temp, z=1-EvalData_stress_matrix.dtdp));
     {axis(1);axis(2)}},
   main=prob_plot_title, 
   xlab="Change in Precipitation (%)", 
   ylab="Change in Temperature (C)",
   cex.lab=1.4,
   xlim=c(-25, 25),ylim=c(0,5)
  )
print(prob_plot)
dev.off()
