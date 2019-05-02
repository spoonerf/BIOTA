library(raster)
library(tidyverse)

pred<-read.csv("Fiona_help_Hansen_dataset/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")
pred<-data.frame(pred$Longitude, pred$Latitude)
names(pred)=c("lon","lat")

hansen<-raster("D:/Fiona/BIOTA/Fiona_help_Hansen_dataset/Hansen_reclass_90.tiff")

buff_crop<-function(x,buff){
  
  xmin<-x$lon -buff
  xmax<-x$lon +buff
  ymin<-x$lat-buff
  ymax<-x$lat +buff
  xmin[xmin < -180]<- -180
  xmax[xmax > 180]<- 180
  ymin[xmin < -60]<- -60
  ymax[xmin < 80]<- 80
  test<-crop(hansen, c(xmin,xmax,ymin,ymax))
  #subset to the cropped area
  test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  test_pts<-test_pts[test_pts[,3] ==1,]
  test_pts <- matrix(test_pts, ncol = 3)
  return(test_pts)
}


hansen_dist_fun<-function(x){
  
  x<-data.frame(x[1], x[2])
  colnames(x)<-c("lon", "lat")
  extract_out<-raster::extract(hansen, x)
  extract_out[is.na(extract_out) | extract_out == 255] <- 0
  dist_out<-NULL
  
  if(extract_out == 1){
    
    dist_out<-0
    
  } else {
    
    test_pts <- buff_crop(x, 0.1)

    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
  }
  
  if(is.null(dist_out)){
    
    test_pts<-buff_crop(x, 0.5)
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
  }
  
  if(is.null(dist_out)){
    
    test_pts<-buff_crop(x, 1)
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
  }
  
  if(is.null(dist_out)){
    test_pts<-buff_crop(x, 4)
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    if(nrow(test_pts>= 1)){
      dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
    }
  }
  
  if(is.null(dist_out)){
    test_pts<-buff_crop(x, 10)

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

fun_out<-apply(pred_u, MARGIN = 1 ,FUN = hansen_dist_fun)
all_dist<-do.call("rbind", fun_out)

pred<-read.csv("PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")
all_dist_out<-merge(pred, all_dist, by=c("lon", "lat"))

write.csv(all_dist_out,"hans_min_dist_sp_90.csv", row.names=FALSE)
