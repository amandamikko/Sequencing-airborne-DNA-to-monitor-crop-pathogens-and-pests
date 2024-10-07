library(compositions)
library(zCompositions)

args <- list("../kraken_reports/all_weeks_count_values_removed_zero_inflated.csv",
             "../kraken_reports/all_weeks_count_values_removed_zero_inflated_gbm.csv")

read.counts.l <- read.csv(unlist(args[1]),row.names=1, header= TRUE)

read.counts.l.gbm <- cmultRepl(X = t(read.counts.l), label = 0, method = "GBM")
write_csv(read.counts.l.gbm, unlist(args[2]))
