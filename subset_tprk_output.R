#subset tprk allreads

library(tidyverse)

reads <- read_csv("allreads_filtered.csv")
absolute <- reads[, c(1, 2, (seq(4, ncol(reads), 2)))]
rel <- reads[, c(1, 2, (seq(3, ncol(reads), 2)))]

write.csv(absolute, "allreads_filt_abs.csv", row.names = FALSE, quote=FALSE)
write.csv(rel, "allreads_filt_rel.csv", row.names = FALSE, quote=FALSE)