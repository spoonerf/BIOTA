wd<-getwd()
.libPaths(c(wd,.libPaths()))

library(raster)
library(tidyverse)
library(SearchTrees)

install.packages('SearchTrees', repos="http://cran.r-project.org")

dir.create("/home/ucfafsp/Scratch/temp_files")

rasterOptions(tmpdir = "/home/ucfafsp/Scratch/temp_files")

pred<-read.csv("PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")
pred<-data.frame(pred$Longitude, pred$Latitude)
names(pred)=c("lon","lat")

hansen<-raster("Hansen_reclass_90.tiff")

buff_crop<-function(x,buff){
  
  xmin<-x$lon -buff
  xmax<-x$lon +buff
  ymin<-x$lat-buff
  ymax<-x$lat +buff
  xmin[xmin < -180]<- -180
  xmax[xmax > 180]<- 180
  ymin[ymin < -60]<- -60
  ymax[ymax > 80]<- 80
  
  #e<-extent(xmin,xmax,ymin,ymax)
  
  test<-crop(hansen, extent(c(xmin,xmax,ymin,ymax)))
  #  test<-ff(getValues(test))
  #subset to the cropped area
  test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  
  test_pts<-test_pts[test_pts[,3] ==1,]
  test_pts <- matrix(test_pts, ncol = 3)
  unlink(dirname(rasterTmpFile()), recursive=TRUE) 
  return(test_pts)
}


hansen_dist_fun<-function(x){
  
  x<-data.frame(x[1], x[2])
  colnames(x)<-c("lon", "lat")
  extract_out<-raster::extract(hansen, x)
  extract_out[is.na(extract_out) | extract_out == 255] <- 0
  dist_out<-NULL
  
  print(x)
  
  if(extract_out == 1){
    
    dist_out<-0
    
    } else {
    
    test_pts <- buff_crop(x, 0.1)
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      #dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
      tree<- createTree(coordinates(test_pts))
      inds<-knnLookup(tree, newdat = coordinates(x), columns = 1:2, k = 1)
      dist_out<-pointDistance(x, test_pts[inds,1:2], lonlat = TRUE)/1000
    }
  }
  
  if(is.null(dist_out)){
    
    test_pts<-buff_crop(x, 0.5)
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      #dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
      tree<- createTree(coordinates(test_pts))
      inds<-knnLookup(tree, newdat = coordinates(x), columns = 1:2, k = 1)
      dist_out<-pointDistance(x, test_pts[inds,1:2], lonlat = TRUE)/1000      
    }
  }
  
  if(is.null(dist_out)){
    
    test_pts<-buff_crop(x, 1)
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      #dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
      tree<- createTree(coordinates(test_pts))
      inds<-knnLookup(tree, newdat = coordinates(x), columns = 1:2, k = 1)
      dist_out<-pointDistance(x, test_pts[inds,1:2], lonlat = TRUE)/1000
  
    }
  }
  
  if(is.null(dist_out)){
    test_pts<-buff_crop(x, 2)
    
    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      #dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
      tree<- createTree(coordinates(test_pts))
      inds<-knnLookup(tree, newdat = coordinates(x), columns = 1:2, k = 1)
      dist_out<-pointDistance(x, test_pts[inds,1:2], lonlat = TRUE)/1000
    }
  }
  
  if(is.null(dist_out)){
    test_pts<-buff_crop(x, 5)	    

    if(nrow(test_pts)==0){
      dist_out<-NULL
    }
    
    if(nrow(test_pts>= 1)){
      #dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
      tree<- createTree(coordinates(test_pts))
      inds<-knnLookup(tree, newdat = coordinates(x), columns = 1:2, k = 1)
      dist_out<-pointDistance(x, test_pts[inds,1:2], lonlat = TRUE)/1000
        }
  }
  
  if(is.null(dist_out)){
    dist_out<-paste0("Further than 5 degrees from nearest forest!")  
  }
  
  dist_out_df<-cbind(x, dist_out)
  print(dist_out_df)
  unlink(dirname(rasterTmpFile()), recursive=TRUE)  
  return(dist_out_df)
 
}

#taking only unique locations to minimise processing time

pred_u<-as.matrix(unique(pred))

fun_out<-apply(pred_u, MARGIN = 1 ,FUN = hansen_dist_fun)
all_dist<-do.call("rbind", fun_out)

#pred<-read.csv("PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")
all_dist_out<-merge(pred, all_dist, by=c("lon", "lat"))

#saveRDS(fun_out, "rds_out_hans_min_dist_sp_90_tree.RDS")
write.csv(all_dist_out, "hans_min_dist_sp_90_tree.csv", row.names =FALSE)

