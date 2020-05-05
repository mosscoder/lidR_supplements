##### LAS Processing helper functions #####

# Check if something is a raster 

is.Raster <- function(x)
{
  return((class(x)[1]=="RasterLayer" || class(x)[1]=="RasterBrick" || class(x)[1]=="RasterStack"))
}

# Plot classified point cloud or raster as a static hillshade map

hillPlot <- function(x, algo = tin(), gamma = 0.35, main = NULL) {
  if(is.Raster(x)) {
    ter <- x
  } else {
    ter <- grid_terrain(x, res = 0.15, algorithm = algo)
  }
  slope <- raster::terrain(ter, opt = 'slope')
  aspect <- raster::terrain(ter, opt = 'aspect')
  hill <- hillShade(slope, aspect)
  plot(hill, col = gray.colors(100, start = 0, gamma = gamma), main = main)
}

# Generate Zhang parameters for four-window progression

make4posZhang <- function(b, c, thMin = 0.1, dhMax = 5) {
  ws <- c(b)
  for(i in 1:3) {
    ws[i + 1] <- 2*i*b + 1
  }
  
  th2 <- (ws[2] - ws[1])*c + thMin
  th3 <- (ws[3] - ws[2])*c + thMin
  th4 <- (ws[4] - ws[3])*c + thMin
  
  th <- c(thMin, th2, th3, ifelse(th4 >= dhMax, dhMax, th4))
  
  zList <- list(ws, th)
  names(zList) <- c('ws', 'th')
  zList
}

# Conduct pmf iteratively

pmfIter <- function(las, iter, zhang) {
  if(iter == 1) {
    lasground(las, algorithm = pmf(zhang$ws, zhang$th))
    } else { 
  
  las@data$ind <- 1:nrow(las@data)
  outTemplate <- las
  outTemplate@data$Classification <- 1L
  las@data$Classification <- 2L
  
  for(i in 1:iter) {
    las <- lasground(lasfilter(las, Classification == 2L), algorithm = pmf(zhang$ws, zhang$th))
  }
  
  groundOnly <- lasfilter(las, Classification == 2L)
  
  outTemplate@data$Classification[groundOnly$ind] <- 2L
  
  outTemplate
  }
}

# Remove lone canopy points within specified radius
filterLoneCanopy <-
  function(las, k = 100, r = 0.5, plotLoners = TRUE) {
    dat <- data.frame(
      ind = seq_along(las@data$X),
      las@data[, c('X', 'Y', 'Classification')],
      loners = 0L
    )
    
    cDat <- dat %>% filter(Classification == 1L)
    gDat <- dat %>% filter(Classification == 2L)
    
    rm(dat)
    gc()
    
    cRows <- nrow(cDat)
    
    if (k > cRows) {
      k <- cRows
    }
    
    canNeighbs <- RANN::nn2(
      data = cDat[, c('X', 'Y')],
      query = cDat[, c('X', 'Y')],
      k = k,
      searchtype = 'radius',
      radius = r,
      treetype = 'kd'
    )$nn.idx[,-1]
    
    groundNeighbs <- RANN::nn2(
      data = gDat[, c('X', 'Y')],
      query = cDat[, c('X', 'Y')],
      k = k,
      searchtype = 'radius',
      radius = r,
      treetype = 'kd'
    )$nn.idx[,-1]
    
    cDat$loners <- ifelse(rowSums(canNeighbs) == 0 &
                                 rowSums(groundNeighbs) > 0, 1, 0)
    rm(canNeighbs, groundNeighbs)
    gc()
    
    if (plotLoners == TRUE) {
      cLas <- lasfilter(las, Classification == 1L)
      cLas@data$loners <- 0L
      cLas@data$loners <- cDat$loners
      plot(cLas, color = 'loners')
    }
    
    las@data$loners <- 0L
    las@data$loners[cDat$ind] <- cDat$loners
    las <- lasfilter(las, loners != 1L)
    
    return(las)
    rm(list = ls())
    gc()
    
  }

expandCanopy <- function(las, k = 100, r = 0.5, plotExpansion = TRUE) {
  dat <- data.frame(
    ind = seq_along(las@data$X),
    las@data[, c('X', 'Y', 'Classification')],
    near = 2L
  )
  
  cDat <- dat %>% filter(Classification == 1L)
  gDat <- dat %>% filter(Classification == 2L)
  
  cRows <- nrow(cDat)
  
  rm(dat)
  gc()
  
  if (k > cRows) {
    k <- cRows
  }
  
  newNeighbs <- RANN::nn2(
    data = gDat[, c('X', 'Y')],
    query = cDat[, c('X', 'Y')],
    k = k,
    searchtype = 'radius',
    radius = r,
    treetype = 'kd'
    
  )
  
  whichNear <- which(newNeighbs$nn.dists > 0 & newNeighbs$nn.dists <= r)
  
  gDat$near[newNeighbs$nn.idx[whichNear]] <- 1L
   
  rm(newNeighbs)
  gc()
  
  if(plotExpansion == TRUE){
    las@data$near <- 0L
    las@data$near[gDat$ind] <- gDat$near
    las@data$near <- ifelse(las@data$near == 2L, 0L, 1L)
    plot(las, color = 'near')
  }
  
  las@data$Classification[gDat$ind] <- gDat$near
  
  return(las)
  
  rm(list = ls())
  gc()
}

flipGroundOutliers <- function(las, k = 100, cutoff = 0.5, plotDumpers = TRUE, buffer = 0) {
  
  dat <- data.frame(
    ind = seq_along(las@data$X),
    las@data[, c('X', 'Y', 'Z', 'Classification')],
    dumpers = 0L
  )
  
  gDat <- dat %>% filter(Classification == 2L)
  
  xmn <- las$X %>% min()
  xmx <- las$X %>% max()
  ymn <- las$Y %>% min()
  ymx <- las$Y %>% max()
  
  gRows <- nrow(gDat)
  
  if (k > gRows) { k <- gRows }
  
  groundNeighbs <- RANN::nn2(
    data = gDat[, c('X', 'Y')],
    query = gDat[, c('X', 'Y')],
    k = k,
    searchtype = 'priority',
    treetype = 'kd'
  )$nn.idx[,-1]
  
  zMat <- matrix(
    nrow = nrow(groundNeighbs),
    ncol = ncol(groundNeighbs),
    byrow = FALSE,
    data = gDat$Z[groundNeighbs[seq_along(groundNeighbs)]]
  )
  
  groundDumpers <-
    ifelse(gDat$Z - rowMeans(zMat, na.rm = T) > cutoff, 1L, 2L)
  
  rm(groundNeighbs, zMat)
  gc
  
  gDat$gDumpers <- groundDumpers
  
  gDat <- gDat %>% filter(
                    X > xmn + buffer &
                      X < xmx - buffer &
                      Y > ymn + buffer &
                      Y < ymx + buffer)
  
  if (plotDumpers == TRUE) {
    las@data$gDumpers <- 0L
    las@data$gDumpers[gDat$ind] <- ifelse(gDat$gDumpers == 2L, 0L, 1L)
    gLas <- lasfilter(las, Classification == 2L)
    plot(gLas, color = 'gDumpers')
  }
  
  las@data$Classification[gDat$ind] <- gDat$gDumpers
  
  return(las)
  
  rm(list = ls())
  gc()
  
}

bicubicSplineInterp <- function(cluster, res, buffer) {
  las <- readLAS(cluster)
  if (is.empty(las))
    return(NULL)
  
  bbox  <- raster::extent(cluster)
  bufferedExtent <- raster::extent(las)
  
  dir <- paste0(tempdir(), '/', bbox[1], '_', bbox[3])
  dir.create(dir)
  
  xyz <-lasfilter(las,
                  Classification == 2L)@data[ , c('X', 'Y', 'Z')] 
  xyz <- st_as_sf(as.data.frame(xyz), 
                  coords = c('X', 'Y'), 
                  crs = crs(las))
  point_dir <- paste0(dir, "/irr_pts.shp")
  st_write(xyz, point_dir, driver = "ESRI Shapefile", overwrite = TRUE)
  
  ras_dir <- paste0(dir, '/interp.sgrd')
  
  system(
    paste0(
      'saga_cmd grid_spline 4', # http://www.saga-gis.org/saga_tool_doc/2.1.3/grid_spline_4.html
      ' -SHAPES ', point_dir, # Location of point cloud in vector form
      ' -FIELD Z', # Layer containing Z
      ' -METHOD 1', # Default 'refinement', seems faster than alternative
      ' -LEVEL_MAX 14', # Highest level ~looks best, doesn't have huge performance disadvantage
      ' -TARGET_USER_XMIN ', bufferedExtent[1],
      ' -TARGET_USER_XMAX ', bufferedExtent[2],
      ' -TARGET_USER_YMIN ', bufferedExtent[3],
      ' -TARGET_USER_YMAX ', bufferedExtent[4],
      ' -TARGET_USER_SIZE ', res, 
      ' -TARGET_OUT_GRID ', ras_dir
    )
  )
  
  interp <- raster(paste0(dir, '/interp.sdat'))
  template <- raster(bufferedExtent, resolution = 0.2, crs = crs(las))
  interp <- resample(interp, template)
  interp <- raster::crop(interp, bbox)
  unlink(dir, recursive=TRUE)
  
  return(interp)
  
}
