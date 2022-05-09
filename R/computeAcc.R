
library(fitdistrplus)
library(ExomeDepth)
library(GenomicRanges)
library(tidyverse)

computeAcc <- function(calls, truth, weighted_regions = NULL) {
  if (!is.null(weighted_regions)) {
    sens <- calls %>% 
      select(start.p, end.p, type) %>% 
      inner_join(truth %>% 
                   filter(floor(start / 10000000) %in% weighted_regions))
  } else {
    sens <- calls %>% 
      select(start.p, end.p, type) %>% 
      inner_join(truth)
  }
  total_sens = nrow(sens) / nrow(truth)
  total_spec = nrow(sens) / nrow(calls)
  dels <- sens %>% 
    filter(type == "deletion")
  delsens <- nrow(dels) / (nrow(truth %>% filter(type == "deletion")))
  delspec = nrow(dels) / (nrow(calls %>% filter(type == "deletion")))
  dups <- sens %>% 
    filter(type == "duplication")
  dupsens <- nrow(dups) / (nrow(truth %>% filter(type == "duplication")))
  dupspec = nrow(dups) / (nrow(calls %>% filter(type == "duplication")))
  return(
    list(
      spec = total_spec,
      sens = total_sens,
      spec_dup = dupspec,
      spec_del = delspec,
      sens_dup = dupsens,
      sens_del = delsens
    )
  )
}