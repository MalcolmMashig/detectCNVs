
source("simulate.R")

spec <- c()
sens <- c()
delspec <- c()
delsens <- c()
dupspec <- c()
dupsens <- c()
for (i in 1:30) {
  print(i)
  acclist <- simulate(weighting = FALSE, 
                      wrtest = c(0, 1, 2, 3, 4, 6, 11, 14, 20, 22))
  spec <- spec %>% c(acclist$spec)
  sens <- sens %>% c(acclist$sens)
  delspec <- delspec %>% c(acclist$spec_del)
  dupspec <- dupspec %>% c(acclist$spec_dup)
  delsens <- delsens %>% c(acclist$sens_del)
  dupsens <- dupsens %>% c(acclist$sens_dup)
}

for (tpi in c(0.00015, 0.0002, 0.00025)) {
  spec2 <- c()
  sens2 <- c()
  delspec2 <- c()
  delsens2 <- c()
  dupspec2 <- c()
  dupsens2 <- c()
  for (i in 1:30) {
    print(i)
    acclist <- simulate(weighting = FALSE, 
                        wrtest = c(0, 1, 2, 3, 4, 6, 11, 14, 20, 22),
                        tp = tpi)
    spec2 <- spec2 %>% c(acclist$spec)
    sens2 <- sens2 %>% c(acclist$sens)
    delspec2 <- delspec2 %>% c(acclist$spec_del)
    dupspec2 <- dupspec2 %>% c(acclist$spec_dup)
    delsens2 <- delsens2 %>% c(acclist$sens_del)
    dupsens <- dupsens %>% c(acclist$sens_dup)
  }
  if (tpi == 0.00015) {
    list0015 <- list(
      spec2,
      sens2,
      delspec2,
      delsens2,
      dupspec2,
      dupsens2
    )
  }
  if (tpi == 0.0002) {
    list002 <- list(
      spec2,
      sens2,
      delspec2,
      delsens2,
      dupspec2,
      dupsens2
    )
  }
  if (tpi == 0.00025) {
    list0025 <- list(
      spec2,
      sens2,
      delspec2,
      delsens2,
      dupspec2,
      dupsens2
    )
  }
}

results <- data.frame(
  Algorithm = c("Baseline", "Baseline", 
                "Account for Active Regions", "Account for Active Regions"),
  Accuracy = c(mean(sens), mean(spec), 
               mean(list0015[[2]]), mean(list0015[[1]])),
  Measure = c("Sensitivity", "Specificity", "Sensitivity", "Specificity")
)
# write_csv(results, "heuristicEval-results.csv")

ggplot(results, aes(Measure, Accuracy, fill = Algorithm)) + 
  geom_bar(stat = 'identity', position = position_dodge()) + 
  ylim(0, 1) +
  theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
