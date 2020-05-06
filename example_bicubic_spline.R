library(lidR)
library(MBA)
library(future)

source('https://raw.githubusercontent.com/mosscoder/lidR_supplements/master/las_helper_functions.R')

set_lidr_threads(1L)
plan('multisession', workers = 4L)

t <- tempfile()

fUrl <- 'https://github.com/mosscoder/lidR_supplements/blob/master/classified_example.las.zip?raw=true'

download.file(url = fUrl,
              destfile = t)

ctg <- catalog(unzip(t)[1])

res <- 0.2

opt_chunk_buffer(ctg) <- 300 #absurd buffer needed to avoid edge artifcats
opt_chunk_size(ctg) <- 300
opt_stop_early(ctg) <- TRUE
opt_select(ctg) <- 'xyz'
opt_filter(ctg) <- '-keep_class 2'

opt <- list(raster_alignment = res,
            automerge = TRUE)

dtm <- catalog_apply(
  ctg,
  FUN = bicubicSplineInterp,
  res = res,
  .options = opt
)

hillPlot(readLAS(ctg@data$filename), main = 'TIN interpolation') 
hillPlot(dtm, main = 'Bicubic spline interpolation') 
