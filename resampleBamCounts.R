
library(fitdistrplus)
# library(ExomeDepth) # uncomment this
library(GenomicRanges)
library(tidyverse)

resampleBamCounts <- function(
  use_pkg_data = FALSE,
  n_regions = round(rnorm(1, 26000, 5000)),
  width_exp_rate = 0.0037,
  count_exp_rate = 0.0083,
  pct_zero_avg = .17,
  pct_zero_sd = .02,
  names = NULL, # TODO
  chromosome = NULL, # TODO
  avg_bin_size = 3,
  sd_bin_size = 5.5,
  n_deletions = 10,
  n_duplications = 10,
  avg_deletion_prop = 0.35,
  sd_deletion_prop = 0.13,
  avg_duplication_prop = 1.7,
  sd_duplication_prop = 0.30
) {
  if (!use_pkg_data) {
    frame <- tibble(.rows = n_regions) %>% 
      mutate(
        seqnames = rep("chr1", n_regions),
        start = rep(0, n_regions),
        end = rep(0, n_regions),
        width = round(rexp(n_regions, rate = width_exp_rate)),
        strand = "*",
        Exome1 = round(rexp(n_regions, rate = count_exp_rate)),
        Exome2 = round(rexp(n_regions, rate = count_exp_rate)),
        Exome3 = round(rexp(n_regions, rate = count_exp_rate)),
        Exome4 = round(rexp(n_regions, rate = count_exp_rate)),
        names = rep("gene", n_regions),
        rn = row_number(),
        bin_end = NA_integer_
      ) %>% 
      mutate(
        Exome4 = ifelse(rn %in% sample(1:n_regions, round(n_regions*.9)),
                        Exome1, Exome4),
        start = cumsum(c(12012, round(rexp(n_regions - 1, rate = 0.0001)))),
        end = start + width
      )
    zeros <- sample(frame$rn, size = round(rnorm(1, pct_zero_avg, pct_zero_sd) * n_regions))
    frame <- frame %>% 
      mutate(
        Exome1 = ifelse(rn %in% zeros, 0, Exome1),
        Exome2 = ifelse(rn %in% zeros, 0, Exome2),
        Exome3 = ifelse(rn %in% zeros, 0, Exome3),
        Exome4 = ifelse(rn %in% zeros, 0, Exome4),
      )
    duplication_starts <- sample(frame$rn, size = n_duplications)
    duplication_bins <- round(rnorm(n = n_duplications, mean = avg_bin_size, 
                                    sd = sd_bin_size))
    duplication_bins <- ifelse(duplication_bins < 0, 0, duplication_bins)
    dups <- tibble(rn = duplication_starts, bins = as.integer(duplication_bins))
    dupframe <- frame %>% 
      left_join(dups)
    
    dupframe <- dupframe %>% 
      mutate(
        bin_end = rn + bins
      )
    
    bins <- dupframe %>% 
      filter(!is.na(bin_end))
    
    CNVs <- bins %>% 
      select(start, end) %>% 
      mutate(
        CNV = "duplication"
      )
    
    duplication_rns <- NULL
    for (i in 1:nrow(bins)) {
      s <- bins[i,]$rn
      e <- bins[i,]$bin_end
      binsize <- bins[i,]$bins
      update <- dupframe %>% 
        filter(between(rn, s, e)) %>% 
        mutate(
          mult = rnorm(binsize + 1, avg_duplication_prop, sd_duplication_prop)
        ) %>% 
        mutate(
          Exome4 = round(((Exome1 + Exome2 + Exome3) / 3 + 5) * mult -
                           5)
        )
      duplication_rns <- duplication_rns %>% c((min(update$rn) - 100):(max(update$rn) + 100))
      dupframe <- dupframe %>% 
        filter(!between(rn, s, e)) %>% 
        bind_rows(update) %>% 
        arrange(rn)
    }
    
    deletion_starts <- sample(setdiff(frame$rn, duplication_rns), 
                              size = n_deletions)
    deletions_bins <- round(rnorm(n = n_deletions, mean = avg_bin_size, 
                                  sd = sd_bin_size))
    deletions_bins <- ifelse(deletions_bins < 0, 0, deletions_bins)
    dels <- tibble(rn = deletion_starts, bins = as.integer(deletions_bins))
    delframe <- frame %>% 
      left_join(dels)
    
    delframe <- delframe %>% 
      mutate(
        bin_end = rn + bins
      )
    
    bins <- delframe %>% 
      filter(!is.na(bin_end))
    
    CNVs <- CNVs %>% 
      bind_rows(
        bins %>% 
          select(start, end) %>% 
          mutate(
            CNV = "deletion"
          )
      )
    
    deletions_rns <- NULL
    for (i in 1:nrow(bins)) {
      s <- bins[i,]$rn
      e <- bins[i,]$bin_end
      binsize <- bins[i,]$bins
      update <- delframe %>% 
        filter(between(rn, s, e)) %>% 
        mutate(
          mult = rnorm(binsize + 1, avg_deletion_prop, sd_deletion_prop)
        ) %>% 
        mutate(
          Exome4 = round(((Exome1 + Exome2 + Exome3) / 3 + 5) * mult),
          Exome4 = ifelse(Exome4 < 0, 0, Exome4)
        )
      deletions_rns <- deletions_rns %>% c((min(update$rn) - 100):(max(update$rn) + 100))
      delframe <- delframe %>% 
        filter(!between(rn, s, e)) %>% 
        bind_rows(update) %>% 
        arrange(rn)
    }
    accuracy <- 
      return(
        list(
          counts = frame %>% 
            full_join(delframe) %>% 
            full_join(dupframe) %>% 
            group_by(rn) %>% 
            filter(
              (n() == 1) | !is.na(mult)
            ) %>% 
            ungroup() %>% 
            arrange(rn) %>% 
            select(-(bin_end:mult)),
          truth = CNVs
        )
      )
  } else {
    data(ExomeCount, package = "ExomeDepth")
    ExomeCount <- data.frame(ExomeCount)
    n <- nrow(ExomeCount)
    new <- ExomeCount %>% 
      mutate(
        Exome1 = round(Exome1 * rnorm(n, 1, 0.05)),
        Exome2 = round(Exome2 * rnorm(n, 1, 0.05)),
        Exome3 = round(Exome3 * rnorm(n, 1, 0.05)),
        Exome4 = round(Exome4 * rnorm(n, 1, 0.05)),
        # across(Exome1:Exome4, function(x) ifelse(x < 0, 0, x))
      )
    return(
      list(
        counts = new,
        truth = read_csv("ground-truth.csv")
      )
    )
  }
}
