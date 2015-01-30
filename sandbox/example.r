# Example R script


library(methods)

cat('Hello, world!\n')

cat('Load oalib\n')

dyn.load(paste("oalibR", .Platform$dynlib.ext, sep=""))

cat('source oalib.R\n')
source("oalib.R")

cat('cacheMetaData\n')
cacheMetaData(1)

cat('Done\n')

