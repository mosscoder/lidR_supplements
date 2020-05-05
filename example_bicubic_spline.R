library(lidR)
source('https://raw.githubusercontent.com/mosscoder/lidR_supplements/master/las_helper_functions.R')
setwd('~/Downloads')

set_lidr_threads(1L)
plan('multisession', workers = 4L)

t <- tempfile()

fUrl <- 'https://github.com/mosscoder/lidR_supplements/blob/master/classified_example.las.zip?raw=true'

download.file(url = fUrl,
              destfile = t)

ctg <- catalog(unzip(t)[1])

  buffer <- 30
  res <- 0.2
  
  opt_chunk_buffer(ctg) <- buffer
  opt_chunk_size(ctg) <- 200
  opt_stop_early(ctg) <- TRUE
  opt_select(ctg) <- 'xyzc'
  opt_output_files(ctg) <- './bicubic/{XLEFT}_{YBOTTOM}_{ID}'
  
  opt <- list(raster_alignment = res,
              automerge = TRUE)
  
  dtm <- catalog_apply(ctg, 
                FUN = bicubicSplineInterp, 
                buffer = buffer,
                res = res, 
                .options = opt)

hillPlot(readLAS(ctg@data$filename), main = 'TIN interpolation') 
hillPlot(dtm, main = 'Bicubic spline interpolation') 
