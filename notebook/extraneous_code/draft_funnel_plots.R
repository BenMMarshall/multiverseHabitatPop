
targets::tar_load("areaBasedResults")

library(dplyr)
library(ggplot2)

areaBasedResults %>% 
  ggplot() +
  geom_point(aes(x = sampleSize, y = companaP))
