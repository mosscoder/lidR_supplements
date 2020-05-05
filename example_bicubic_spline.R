library(lidR)
source('~/mpgPostdoc/projects/bareEarth/analysis/support/las_helper_functions.R')

setwd('~/Downloads')

t <- tempfile()

download.file(url = 'https://storage.googleapis.com/mpg-files/classified_example.las.zip',
              destfile = t)

ctg <- catalog(unzip(t)[1])

buffer <- 30
res <- 0.2

opt_chunk_buffer(ctg) <- buffer
opt_chunk_size(ctg) <- 200
opt_stop_early(ctg) <- TRUE
opt_select(ctg) <- 'xyzc'
opt_output_files(ctg) <- './bicubic/{XLEFT}_{YBOTTOM}_{ID}'

set_lidr_threads(1L)
data.table::setDTthreads(1L)
plan('multisession', workers = 4L)

opt <- list(raster_alignment = res,
            automerge = TRUE)

dtm <- catalog_apply(ctg, 
              FUN = bicubicSplineInterp, 
              buffer = buffer,
              res = res, 
              .options = opt)

hillPlot(dtm)
