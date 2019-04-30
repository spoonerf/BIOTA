library(raster)
library(tidyverse)
library(here)

pred<-read.csv(paste0(here(),"/Fiona_help_Hansen_dataset/PREDICTS_NatPlusCrop_forestBiome_site_level_new.csv"))
pred<-data.frame(pred$Longitude, pred$Latitude)
names(pred)=c("lon","lat")

hansen<-raster(paste0(here(), "/Fiona_help_Hansen_dataset/Hansen_reclass.tiff"))


hansen_dist_fun<-function(x){
  
  x<-data.frame(x[1], x[2])
  colnames(x)<-c("lon", "lat")
  extract_out<-raster::extract(hansen, x)
  extract_out[is.na(extract_out) | extract_out == 255] <- 0
  dist_out<-NULL
  
  if(extract_out == 1){
    
    dist_out<-0
    
  } else {
    
    xmin<-x$lon -0.1
    xmax<-x$lon +0.1
    ymin<-x$lat-0.1
    ymax<-x$lat +0.1
    test<-crop(hansen, c(xmin,xmax,ymin,ymax))
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
    test_pts<-test_pts[test_pts[,3] ==1,]
  
    test_pts <- matrix(test_pts, ncol = 3)
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
  }
  
  if(is.null(dist_out)){
    
    xmin<-x$lon -0.5
    xmax<-x$lon +0.5
    ymin<-x$lat -0.5
    ymax<-x$lat +0.5
    test<-crop(hansen, c(xmin,xmax,ymin,ymax))
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
    test_pts<-test_pts[test_pts[,3] ==1,]
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
    
    }
  
  if(is.null(dist_out)){
    
    xmin<-x$lon -1
    xmax<-x$lon +1
    ymin<-x$lat -1
    ymax<-x$lat +1
    test<-crop(hansen, c(xmin,xmax,ymin,ymax))
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
    test_pts<-test_pts[test_pts[,3] ==1,]
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
    }
  
  if(is.null(dist_out)){
    #have capped the maximum cropping distance to 4 degrees so that we don't end up cropping beyond the extent of the original raster e.g. values >180 degrees
    xmin<-x$lon -4
    xmax<-x$lon +4
    ymin<-x$lat -4
    ymax<-x$lat +4
    test<-crop(hansen, c(xmin,xmax,ymin,ymax))
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
    test_pts<-test_pts[test_pts[,3] ==1,]
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
  }
  
  if(is.null(dist_out)){
    #have capped the maximum cropping distance to 4 degrees so that we don't end up cropping beyond the extent of the original raster e.g. values >180 degrees
    xmin<-x$lon -10
    xmax<-x$lon +10
    ymin<-x$lat -10
    ymax<-x$lat +10
    #making sure the extents of the crop don't go over the edge of the hansen data
    xmin[xmin < -180]<- -180
    xmax[xmax > 180]<- 180
    ymin[xmin < -60]<- -60
    ymax[xmin < 80]<- 80
    test<-crop(hansen, c(xmin,xmax,ymin,ymax))
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
    test_pts<-test_pts[test_pts[,3] ==1,]
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
  }
  
  if(is.null(dist_out)){
    dist_out<-paste0("Further than 10 degrees from nearest forest!")  
  }
  
  dist_out_df<-cbind(x, dist_out)
  print(dist_out_df)
  return(dist_out_df)
  
}
#taking only unique locations to minimise processing time

pred_u<-as.matrix(unique(pred))

fun_out<-apply(pred_u[3680:nrow(pred_u),], MARGIN = 1 ,FUN = hansen_dist_fun)
all_dist<-do.call("rbind", fun_out)

pred<-read.csv(paste0(here(),"/Fiona_help_Hansen_dataset/PREDICTS_NatPlusCrop_forestBiome_site_level_new.csv"))
all_dist_out<-merge(pred, all_dist, by=c("lon", "lat"))

write.csv(all_dist_out,"hans_min_dist_sp_1.csv", row.names=FALSE)


