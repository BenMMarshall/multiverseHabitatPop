
targets::tar_load("areaBasedAvailUse_15_1")

areaBasedAvailUse_15_1

# ARGUEMENTS
# avialUseData = areaBasedAvailUse,
# sampleGroups = optionsList_samples,
# optionsList = optionsList_areaMethods

library(adehabitatHS)
library(dplyr)

use <- areaBasedAvailUse_15_1 %>% 
  filter(type == "III",
         method == "MCP",
         contour == 95,
         aPoints == 1,
         spSamp == "rd") %>% 
  select(used_c0, used_c2)

avail <- areaBasedAvailUse_15_1 %>% 
  filter(type == "III",
         method == "MCP",
         contour == 95,
         aPoints == 1,
         spSamp == "rd") %>% 
  select(avail_c0, avail_c2)

names(use) <- c("c0", "c2")
names(avail) <- c("c0", "c2")

# eisera gives individual selection, but require multiple individuals to run
eiseraOUT <- eisera(used = use, available = avail, scannf = FALSE)
eiseraOUT$wij


companaOUT <- compana(used = use, avail = avail,
                      test = "randomisation")

companaOUT$test
