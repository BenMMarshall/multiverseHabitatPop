#!/usr/bin/env Rscript

# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

# targets::tar_make()
# targets::tar_make_clustermq(workers = 2) # nolint
targets::tar_make_clustermq(names = contains("OUT"), workers = 2) # nolint
# targets::tar_make_future(workers = 2) # nolint

# Examine -----------------------------------------------------------------

# targets::tar_make(starts_with("sampDuraFreq"))
# targets::tar_make(ends_with("esults"))
# targets::tar_make_clustermq(names = ends_with("esults"), workers = 8, log_worker = TRUE)
# targets::tar_make_clustermq(workers = 2, log_worker = TRUE)
# targets::tar_visnetwork(names = starts_with("combined"))
# mani <- targets::tar_manifest()
# pattern <- paste0("ssfOUT.*(",
#                   paste(sample(1:4, 2, replace = FALSE), collapse = "|"),
#                   ")_BADGER")
# mani$name[grep("ssfOUT.*4_BADGER", mani$name)]
# mani$name[grep(pattern, mani$name)]

# "poisOUT_K_30_12"
# "twoStepOUT_K_30_12"