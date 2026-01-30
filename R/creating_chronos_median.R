library(readr)
library(dplyr)
chronos <- read_csv("CRISPR_(DepMap_Public_25Q3+Score,_Chronos)_subsetted.csv")

chronos_median <- chronos %>%
  summarise(across(where(is.numeric), median))

chronos_median_final <- as.data.frame(t(chronos_median))
chronos_median_final$essential <- c('0')
colnames(chronos_median_final)[-1] <- "Essential_label"

chronos_median_final$Essential_label <- ifelse(
  chronos_median_final$Chronos_median < -0.5,
  1,
  0
)
table(chronos_median_final$Essential_label)
head(chronos_median_final[order(chronos_median_final$Chronos_median), ])
write_csv2(chronos_median_final, "chronos_median")
