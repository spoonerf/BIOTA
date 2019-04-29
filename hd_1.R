library(raster)
library(tidyverse)
library(here)

pred<-read.csv(paste0(here(),"/Fiona_help_Hansen_dataset/PREDICTS_NatPlusCrop_forestBiome_site_level_new.csv"))
pred<-data.frame(pred$Longitude, pred$Latitude)
names(pred)=c("lon","lat")

hansen<-raster(paste0(here(), "/Fiona_help_Hansen_dataset/Hansen_reclass.tiff"))
#min_dist<-vector()



hansen_dist_fun<-function(x){
  
  x<-data.frame(x[1], x[2])
  colnames(x)<-c("lon", "lat")
  extract_out<-raster::extract(hansen, x)
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
    dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
  }
  
  if(is.null(dist_out)){
    
    xmin<-x$lon -2
    xmax<-x$lon +2
    ymin<-x$lat -2
    ymax<-x$lat +2
    test<-crop(hansen, c(xmin,xmax,ymin,ymax))
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
    test_pts<-test_pts[test_pts[,3] ==1,]
    dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
  }
  
  if(is.null(dist_out)){
    
    xmin<-x$lon -4
    xmax<-x$lon +4
    ymin<-x$lat -4
    ymax<-x$lat +4
    test<-crop(hansen, c(xmin,xmax,ymin,ymax))
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
    test_pts<-test_pts[test_pts[,3] ==1,]
    dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
  }
  
  
  dist_out_df<-cbind(x, dist_out)
  print(dist_out_df)
  return(dist_out_df)
  
}

pred_u<-as.matrix(unique(pred))

test_fun<-apply(pred_u, MARGIN = 1 ,FUN = hansen_dist_fun)



write.csv(predog,"hans_min_dist_sp_1.csv", row.names=FALSE)
