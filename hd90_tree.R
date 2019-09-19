####################################################
####################################################
#####                                          #####   
##### A function for calculating the distance  #####
##### between a set of locations e.g. PREDICTS #####  
##### sites and the nearest non-NA cell in a   #####  
##### raster, e.g. a forest raster.            #####
#####                                          #####
#####                                          #####  
#####                                          #####
#####                                          #####
####################################################
####################################################
# wd <- getwd()
# .libPaths(c(wd, .libPaths()))

#install.packages('SearchTrees', repos = "http://cran.r-project.org")

library(raster)
library(tidyverse)
library(SearchTrees)
library(parallel)

#dir.create("/home/ucfafsp/Scratch/temp_files")  # Use if running on Myriad - rasters are often stored in a temp folder so if handling a lot of them you need to be
#careful your temp folder doesn't fill up

#rasterOptions(tmpdir = "/home/ucfafsp/Scratch/temp_files")  #Setting the raster temp file directory

bs <-
  c(0.1, 0.5, 1, 2)   # Range of buffer sizes for the function to iterate, measured in degrees. I
                      # If a point is further than the largest buffer then it will throw a message
                      # like "Further than X degrees to the nearest forest!

pred <-
  read.csv("PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")  # Reading in the PREDICTS data but you just need lon/lat coords
pred <- data.frame(pred$Longitude, pred$Latitude)
names(pred) = c("lon", "lat")

hansen <- raster("Hansen_reclass_90.tiff")   #  A raster to get distances from, should ideally be binary 1s and 0s where we are looking for distances to 1s 

buff_crop <- function(x, buff) {
  xmin <- x$lon - buff   #Adding the buffer to the point
  xmax <- x$lon + buff
  ymin <- x$lat - buff
  ymax <- x$lat + buff
  xmin[xmin < -180] <-
    -180    #Ensuring the extent of the crop does not go beyond the bounds of the original raster
  xmax[xmax > 180] <- 180
  ymin[ymin < -60] <- -60
  ymax[ymax > 80] <- 80
  
  
  test <- crop(hansen, raster::extent(c(xmin, xmax, ymin, ymax)))
  test_pts <- rasterToPoints(test, function(x) {
    !is.na(x)
  })
  
  test_pts <- test_pts[test_pts[, 3] == 1, ]   ##Here we are looking for raster cells with the value of one but this could be changed to be more flexible
  test_pts <- matrix(test_pts, ncol = 3)
  unlink(dirname(rasterTmpFile()), recursive = TRUE)    #Deletes temp folders and files which can accumulate rapidly in raster analysis
  return(test_pts)
}

hansen_dist_fun <- function(x) {
  x <- data.frame(x[1], x[2])
  colnames(x) <- c("lon", "lat")
  extract_out <- raster::extract(hansen, x)
  extract_out[is.na(extract_out) | extract_out == 255] <- 0   # This may need changing if your raster has different values. Here we are looking for 1s only
  
  dist_out <- NULL
  
  if (extract_out == 1) {
    dist_out <-
      0   #If the cell the point is on is forest then the distance to forest is zero
    bs_out <- 0
    
  } else {
    i <- 1
    while (is.null(dist_out)) {
      
      if(i > length(bs)){
        dist_out <- paste0("Location is further than ", max(bs), " degrees from nearest forest!")
        bs_out <- NA
      } else {
        
      
      test_pts <-
        buff_crop(x, bs[i])  #If the cell isn't forest then we crop the raster with a buffer using a while loop.
                              #It iterates throught the buffer sizes until  we find a cell of forest and therefore "dist_out" is no longer null
      if (nrow(test_pts) == 0) {
        dist_out <- NULL
      }
      
      if (nrow(test_pts) >= 1) {
        tree <- createTree(coordinates(test_pts))
        inds <-
          knnLookup(tree,
                    newdat = coordinates(x),
                    columns = 1:2,
                    k = 1)
        dist_out <-
          pointDistance(x, test_pts[inds, 1:2], lonlat = TRUE) / 1000   ##Distance out is calculated in km
      }
      
      bs_out <- bs[i]
      i = i + 1   #if the distance out is "NULL" then iterate up to a larger buffer size
      }
    }
      
  }
  

  dist_out_df <- cbind(x, bs_out, dist_out)
  print(dist_out_df)
  unlink(dirname(rasterTmpFile()), recursive = TRUE)  #Deletes raster temp files so they don't clog up the memory
  return(dist_out_df)
  
}

pred_u <-
  as.matrix(unique(pred)) #taking only unique locations to minimise processing time

fun_out <-
  apply(pred_u, MARGIN = 1 , FUN = hansen_dist_fun)  #Applying the function through the matrix of lon/lats


# cl <- makeCluster(getOption("cl.cores", detectCores()-1))
# clusterExport(cl=cl, varlist=c("hansen", "bs", "buff_crop", "hansen_dist_fun", "raster"), envir=environment())
# 
# fun_out <-
#   parApply(cl, pred_u[azores_rows,], MARGIN = 1 , FUN = hansen_dist_fun)  #Applying the function through the matrix of lon/lats
# 
all_dist <- do.call("rbind", fun_out)  #rbinding the results together

#write.csv(dist_out, "hans_min_dist_sp_90_tree.csv", row.names = FALSE)


