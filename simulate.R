
library(fitdistrplus)
# library(ExomeDepth) # uncomment this
library(GenomicRanges)
library(tidyverse)

source("resampleBamCounts.R")
source("computeAcc.R")

simulate <- function(weighting = FALSE, weighted_regions = NULL, tp = 10^-4, wrtest = NULL) {
  sim <- resampleBamCounts(use_pkg_data = TRUE)
  ExomeCount.dafr <- sim$counts
  truth <- sim$truth
  if (weighting) {
    split <- ExomeCount.dafr %>% 
      mutate(
        region = floor(start /10000000)
      )
    split <- split(split, split$region)
    calls <- NULL
    for (i in 1:length(split)) {
      ExomeCount.dafr <- split[[i]]
      my.test <- ExomeCount.dafr$Exome4
      my.ref.samples <- c('Exome1', 'Exome2', 'Exome3')
      my.reference.set <- as.matrix(ExomeCount.dafr[, my.ref.samples])
      my.choice <- select.reference.set (test.counts = my.test,
                                         reference.counts = my.reference.set,
                                         bin.length = ExomeCount.dafr$width/1000,
                                         n.bins.reduced = 10000)
      my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])
      my.reference.selected <- apply(X = my.matrix,
                                     MAR = 1,
                                     FUN = sum)
      all.exons <- new('ExomeDepth',
                       test = my.test,
                       reference = my.reference.selected,
                       formula = 'cbind(test, reference) ~ 1')
      if (i %in% weighted_regions) {
        tp = tp
      } else {
        tp = 10^-4
      }
      all.exons <- CallCNVs(x = all.exons,
                            transition.probability = tp,
                            chromosome = ExomeCount.dafr$seqnames,
                            start = ExomeCount.dafr$start,
                            end = ExomeCount.dafr$end,
                            name = ExomeCount.dafr$names)
      calls <- calls %>% bind_rows(
        tibble(all.exons@CNV.calls)
      )
    }
  } else {
    my.test <- ExomeCount.dafr$Exome4
    my.ref.samples <- c('Exome1', 'Exome2', 'Exome3')
    my.reference.set <- as.matrix(ExomeCount.dafr[, my.ref.samples])
    my.choice <- select.reference.set (test.counts = my.test,
                                       reference.counts = my.reference.set,
                                       bin.length = ExomeCount.dafr$width/1000,
                                       n.bins.reduced = 10000)
    my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])
    my.reference.selected <- apply(X = my.matrix,
                                   MAR = 1,
                                   FUN = sum)
    all.exons <- new('ExomeDepth',
                     test = my.test,
                     reference = my.reference.selected,
                     formula = 'cbind(test, reference) ~ 1')
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = tp,
                          chromosome = ExomeCount.dafr$seqnames,
                          start = ExomeCount.dafr$start,
                          end = ExomeCount.dafr$end,
                          name = ExomeCount.dafr$names)
    calls <- tibble(all.exons@CNV.calls)
  }
  if (!is.null(wrtest)) {
    computeAcc(calls, truth, weighted_regions = wrtest)
  } else {
    computeAcc(calls, truth)
  }
}