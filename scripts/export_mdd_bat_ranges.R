#!/usr/bin/env Rscript
# Convert the checksum-verified downscaled MDD order file, without publishing geometry.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
library(terra)
x <- readRDS(args[1])
if (inherits(x,"PackedSpatVector")) x <- unwrap(x)
stopifnot(inherits(x,"SpatVector"), nrow(x)>500L, "sciname" %in% names(x))
stopifnot(all(!is.na(x$sciname)), all(nzchar(x$sciname)), nzchar(crs(x)))
cat("MDD_RANGE_SCHEMA",paste(names(x),collapse="|"),"rows",nrow(x),"unique names",length(unique(x$sciname)),"\n")
writeVector(x[,"sciname"], args[2], filetype="GPKG", overwrite=TRUE)
