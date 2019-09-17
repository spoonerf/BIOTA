wd<-getwd()    
.libPaths(c(wd,.libPaths()))

install.packages('SearchTrees', repos="http://cran.r-project.org")

library(raster)
library(tidyverse)
library(SearchTrees)

dir.create("/home/ucfafsp/Scratch/temp_files")  # Rasters are often stored in a temp folder so if handling a lot of them you need to be 
                                                  #careful your temp folder doesn't fill up

rasterOptions(tmpdir = "/home/ucfafsp/Scratch/temp_files")

pred<-read.csv("PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")
pred<-data.frame(pred$Longitude, pred$Latitude)
names(pred)=c("lon","lat")

hansen<-raster("Hansen_reclass_90.tiff")

buff_crop<-function(x,buff){
  
  xmin<-x$lon -buff   #Adding the buffer to the point
  xmax<-x$lon +buff
  ymin<-x$lat-buff
  ymax<-x$lat +buff
  xmin[xmin < -180]<- -180    #Ensuring the extent of the crop does not go beyond the bounds of the original raster
  xmax[xmax > 180]<- 180
  ymin[ymin < -60]<- -60
  ymax[ymax > 80]<- 80
  
  
  test<-crop(hansen, extent(c(xmin,xmax,ymin,ymax)))
  test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  
  test_pts<-test_pts[test_pts[,3] ==1,]
  test_pts <- matrix(test_pts, ncol = 3)
  unlink(dirname(rasterTmpFile()), recursive=TRUE)    #Deletes temp folders and files which can accumulate rapidly in raster analysis 
  return(test_pts)
}

bs<-c(0.1,0.5,1,2,5)   # range of buffers to iterate through in degrees


hansen_dist_fun<-function(x){
  
  x<-data.frame(x[1], x[2])
  colnames(x)<-c("lon", "lat")
  extract_out<-raster::extract(hansen, x)
  extract_out[is.na(extract_out) | extract_out == 255] <- 0
  
  dist_out<-NULL
  
  if(extract_out == 1){   
    
    dist_out<-0   #If the cell the point is on is forest then the distance to forest is zero
    
   } else {
    
     i <- 1
    while(is.null(dist_out)){
      test_pts <- buff_crop(x, bs[i])  #If it isn't then we expand the buffer using a while loop. It iterates throught the buffer sizes until "dist_out" is no longer null   
      
      if(nrow(test_pts)==0){
        dist_out<-NULL
        }
      
      if(nrow(test_pts) >= 1){
        #dist_out<-min(pointDistance(x, test_pts[,1:2], lonlat = TRUE)/1000)
        tree<- createTree(coordinates(test_pts))
        inds<-knnLookup(tree, newdat = coordinates(x), columns = 1:2, k = 1)
        dist_out<-pointDistance(x, test_pts[inds,1:2], lonlat = TRUE)/1000   ##Distance out is calculated in km 
      }
      #print(bs[i])
      bs_out<-bs[i]
      i = i + 1   #if the distance out is "NULL" then iterate up to a larger buffer size
       }
      
    }
  
  
  if(is.null(dist_out)){
    dist_out<-paste0("Further than 5 degrees from nearest forest!")  
  }
  
  dist_out_df<-cbind(x, bs_out, dist_out)
  print(dist_out_df)
  unlink(dirname(rasterTmpFile()), recursive=TRUE)  
  return(dist_out_df)
 
}

pred_u<-as.matrix(unique(pred)) #taking only unique locations to minimise processing time

fun_out<-apply(pred_u[639:nrow(pred_u),], MARGIN = 1 ,FUN = hansen_dist_fun)  #Applying the function through the matrix of lon/lats
all_dist<-do.call("rbind", fun_out)  #rbinding the results together

#pred<-read.csv("PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")
all_dist_out<-merge(pred, all_dist, by=c("lon", "lat"))

#saveRDS(fun_out, "rds_out_hans_min_dist_sp_90_tree.RDS")
write.csv(all_dist_out, "hans_min_dist_sp_90_tree.csv", row.names =FALSE)

