
# Code orignally from:
#' @article{doi:10.1111/1365-2656.13087,
#'   author = {Muff, Stefanie and Signer, Johannes and Fieberg, John},
#'   title = {Accounting for individual-specific variation in habitat-selection studies: Efficient estimation of mixed-effects models using Bayesian or frequentist computation},
#'   journal = {Journal of Animal Ecology},
#'   volume = {89},
#'   number = {1},
#'   pages = {80-92},
#'   doi = {10.1111/1365-2656.13087},
#'   url = {https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2656.13087},
#'   eprint = {https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/1365-2656.13087},
#'   year = {2020}
#' }
#' 

# Edited by Benjamin Michael Marshall and Samantha Nicole Smith 2020-09-02 for publication of: 
# title = {Restricted constrictors: Space use and habitat selection of Burmese pythons in Northeast Thailand}
## authors= {Samantha Nicole Smith, Max Dolton Jones, Benjamin Michael Marshall, Surachit Waengsothorn, George A. Gale and Colin Thomas Strine}
## year = {2020}

library(dplyr)
library(raster)
library(amt)
library(ggplot2)
library(readr)

# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#install.packages("INLA", lib = .libPaths()[3], 
#                 repos=c(getOption("repos"),
#                 INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)


# Read in animal tracking data --------------------------------------------

#snake.data <- read_csv(file = "BUCA_data_ISSF_M01-M36.csv")
snake.data <- read_csv(file = "BUCA_data_complete.csv",
                       locale = locale(tz = "Asia/Bangkok"))

# Set datetime as "POSIXct" "POSIXt" format
class(snake.data$datetime)
snake.data$datetime <- as.POSIXct(snake.data$datetime, format = "%m/%d/%Y %H:%M")
class(snake.data$datetime)

# view the data
snake.data <- snake.data[snake.data$animal %in% c("M02"),]

# Since BUCA32 and BUCA35 have some locations in 47N and 48N
#subset out the different UTM zone locations
utm47data <- snake.data[snake.data$UTMzone == "47N",]
# tell R the utm zone that the data was collected in, and make a SP object
sp47 <- SpatialPoints(coords= cbind(utm47data[,3], utm47data[,4]), proj4string = CRS("+init=epsg:32647"))
# now that R knows it was collected in 47N we can convert it to 48N
sp48convert <- spTransform(sp47, CRS("+init=epsg:32648"))
# extract the newly converted locations and place them back in the orignal dataframe
snake.data$x[snake.data$UTMzone == "47N"] <- sp48convert@coords[,1]
snake.data$y[snake.data$UTMzone == "47N"] <- sp48convert@coords[,2]

##### set projection system
crs.proj <- CRS("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

dat <- snake.data

# Load in raster covariates -----------------------------------------------

## reading inverted tifs from ISSF code
covariates <- list("Settlement_dist_INVERT.tif",
                   "Road_dist_INVERT.tif",
                   "Building_dist_INVERT.tif",
                   "Agriculture_dist_INVERT.tif",
                   "Less_Distured_Areas_dist_INVERT.tif")


rp <- stack(x = covariates)


cov.names <- c("dist_settle", "dist_road", "dist_building", "dist_ag", 
               "dist_nat")
names(rp) <- cov.names 

# START OF IMPORTED (and edited) CODE from Muff et al --------------------------------------------------

# select only the information we need
dat <- dat %>% 
  select("x" = x, "y" = y,
         "t" = datetime, "id" = animal) %>% 
  group_by(id) %>% 
  arrange(t)

# separate the individuals out into a list-column of dataframes, each item an animal
dat_all <- dat %>% 
  nest(-id) 

# if you want to add meta data, like sex, do so here. We had a second dataframe
# that included metadata

#dat_all <- left_join(dat_all,
#                     snake.meta %>% 
#                       select(id = deployment.id, sex = animal.sex))

# map operates a lot like an apply or loop. It repeats the function to each item
# in the list. In this case we make all the individual dataframes into track
# objects. Check dat_all to see the list of dataframes now has a second
# dataframe/tibble for the track.
dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, x, y, t, crs = crs.proj)
  }))

# Here the summarize_sampling_rate is repeated on each track object to give you
# an individual level summary.
dat_all %>% 
  mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  select(id, sr) %>%
  unnest(cols = c(sr))

# So dat_all is set up ready for the generation of random steps and the
# extraction of raster covariate values

# Step-Selection function -------------------------------------------------

# To make the random steps we use the map function again to do this on an
# individual basis. We calculate the step lengths so we have data to draw from
# to generate the random locations. Then we drop out the zero step lengths that
# break things, or we can drop the sub-GPS error points because we cannot be
# confident they were moves (ie not a new selection choice from the animal).
# Then the random steps, can be generous here because we aren't using GPS data
# with the high resolution that often entails. We have use more randoms to pick
# up smaller changes and maximise the high res rasters we have. Final step is to
# get the covariate values. After the mutate+map bit is a few lines to compile
# the data into a single v large dataframe for the model.

# partprep <- dat_all %>% 
#   mutate(stps = map(trk, ~ .x %>%
#                       steps() %>% 
#                       filter(sl_>0) %>% # removing the non-moves, or under GPS error
#                       random_steps(n = 10) %>% 
#                       extract_covariates(rp))) %>% 
#   select(id, stps) %>%
#   unnest(cols = c(stps))
# 
# unique(partprep$id)
# as.numeric(as.factor(unique(partprep$id)))
# as.numeric(as.factor(partprep$id))
# 
# partprep %>% 
#   ungroup() %>% 
#   mutate(id = as.numeric(as.factor(id))) %>% 
#   pull(id) %>% 
#   unique()

dat_ssf <- dat_all %>% 
  mutate(stps = map(trk, ~ .x %>%
                      steps() %>% 
                      filter(sl_>0) %>% # removing the non-moves, or under GPS error
                      random_steps(n = 10) %>% 
                      extract_covariates(rp))) %>% 
  select(id, stps) %>%
  unnest(cols = c(stps)) %>% 
  group_by(id) %>% 
  arrange(t1_) %>% 
  ungroup() %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"),
    cos_ta = cos(ta_), 
    log_sl = log(sl_))

# So dat_ssf is the correct format for the modelling. We have the id, a heap of
# movement columns, and the critical case_ that describes whether they used a
# location or not.


# Running the INLA model --------------------------------------------------

# We can run the INLA model using the priors and set-up from Muff et al.
# Precision for the priors of slope coefficients
prec.beta.trls <- 1e-4

# "In the model formula for INLA, we set the stratum-specific intercept variance
# to $10^6$ (or rather: the precision to $10^{ ^'6}$) by fixing it (`fixed=T`)
# to an initial value. The other precision is given a PC(1,0.05) prior:"

###### For BUCA data
### making formulas for each habitat feature

#Natural
formula.random.n <- y ~  -1 + 
  dist_nat + # fixed covariate effect
  dist_nat:log_sl + # covar iteractions
  dist_nat:cos_ta + 
  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id, dist_nat, values = 1:length(unique(dat_ssf$id)), model="iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(1, 0.05)))) 
# settlement
formula.random.s <- y ~  -1 + 
  dist_settle + # fixed covariate effect
  dist_settle:log_sl + # covar iteractions
  dist_settle:cos_ta + 
  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id, dist_settle, values = 1:length(unique(dat_ssf$id)), model="iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(1, 0.05)))) 

# road 
formula.random.r <- y ~  -1 + 
  dist_road + # fixed covariate effect
  dist_road:log_sl + # covar iteractions
  dist_road:cos_ta + 
  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id, dist_road, values = 1:length(unique(dat_ssf$id)), model="iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(1, 0.05)))) 

# Buildings
formula.random.w <- y ~  -1 + 
  dist_building + # fixed covariate effect
  dist_building:log_sl + # covar iteractions
  dist_building:cos_ta + 
  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id, dist_building, values = 1:length(unique(dat_ssf$id)), model="iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(1, 0.05)))) 

# Agriculture
formula.random.aq <- y ~  -1 + 
  dist_ag + # fixed covariate effect
  dist_ag:log_sl + # covar iteractions
  dist_ag:cos_ta + 
  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id, dist_aq.ag, values = 1:length(unique(dat_ssf$id)), model="iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(1, 0.05)))) 

# Last two parts of the model are the gaussian processes to deal with step_id
# (each move), and the id of the animal. So steps_id is the covariate, then in
# the second one is id that is weighted by the raster cov.

# Fit the models -----------------------------------------------------------

##############################################################################################################

#########################   V      ISSUES HERE      V     ISSUES HERE       V     ISSUES HERE

##############################################################################################################

# natural model
inla.ssf.n <- inla(formula.random.n, family = "Poisson", data = dat_ssf, #verbose=TRUE,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = prec.beta.trls)))

# settlement model 
inla.ssf.s <- inla(formula.random.s, family = "Poisson", data = dat_ssf, verbose=TRUE, 
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = prec.beta.trls))
)

# road model 
inla.ssf.r <- inla(formula.random.r, family = "Poisson", data = dat_ssf,  verbose=TRUE,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = prec.beta.trls))
)

# building model 
inla.ssf.b <- inla(formula.random.b, family = "Poisson", data = dat_ssf,  verbose=TRUE,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = prec.beta.trls))
)

# agriculture
inla.ssf.a <- inla(formula.random.ag, family = "Poisson", data = dat_ssf, 
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = prec.beta.trls))
)



# Model results -----------------------------------------------------------

# We can see a dataframe of the fixed effects here, so each covar and any interaction terms
inla.ssf.n$summary.fixed
inla.ssf.s$summary.fixed
inla.ssf.r$summary.fixed
inla.ssf.b$summary.fixed
inla.ssf.a$summary.fixed


# "Since variances are parameterized and treated as precisions, the summary of
# the respective posterior distributions is given for the precisions:"
inla.ssf.n$summary.hyperpar
inla.ssf.s$summary.hyperpar
inla.ssf.r$summary.hyperpar
inla.ssf.b$summary.hyperpar
inla.ssf.a$summary.hyperpar




inla_mmarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.mmarginal(inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mode of variance"))
  results
}
inla_emarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mean of variance"))
  results
}

# "Posterior mean and mode are obtained as"
inla_emarginal(inla.ssf.n)
inla_mmarginal(inla.ssf.n)


inla_emarginal(inla.ssf.s)
inla_mmarginal(inla.ssf.s)


inla_emarginal(inla.ssf.r)
inla_mmarginal(inla.ssf.r)


inla_emarginal(inla.ssf.b)
inla_mmarginal(inla.ssf.b)


inla_emarginal(inla.ssf.a)
inla_mmarginal(inla.ssf.a)



fixed.df.n <- inla.ssf.n$summary.fixed
fixed.df.s <- inla.ssf.s$summary.fixed
fixed.df.r <- inla.ssf.r$summary.fixed
fixed.df.b <- inla.ssf.b$summary.fixed
fixed.df.a <- inla.ssf.a$summary.fixed




fixed.df.n <- inla.ssf.n$summary.fixed
fixed.df.n$term <- row.names(fixed.df.n)
names(fixed.df.n) <- c("mean", "sd", "q025", "q50", "q975",
                       "mode", "kld", "term")


fixed.df.s <- inla.ssf.s$summary.fixed
fixed.df.s$term <- row.names(fixed.df.s)
names(fixed.df.s) <- c("mean", "sd", "q025", "q50", "q975",
                       "mode", "kld", "term")


fixed.df.r <- inla.ssf.r$summary.fixed
fixed.df.r$term <- row.names(fixed.df.r)
names(fixed.df.r) <- c("mean", "sd", "q025", "q50", "q975",
                       "mode", "kld", "term")


fixed.df.b <- inla.ssf.b$summary.fixed
fixed.df.b$term <- row.names(fixed.df.b)
names(fixed.df.b) <- c("mean", "sd", "q025", "q50", "q975",
                       "mode", "kld", "term")

fixed.df.a <- inla.ssf.a$summary.fixed
fixed.df.a$term <- row.names(fixed.df.a)
names(fixed.df.a) <- c("mean", "sd", "q025", "q50", "q975",
                       "mode", "kld", "term")


allmods <- rbind(fixed.df.n, fixed.df.s, fixed.df.r, fixed.df.b, fixed.df.a)


write_csv(allmods, "Pop_ssf_mods.csv")



##### figures for population level ssf
popssf <- read_csv(file = "Pop_ssf_mods.csv")


#### Distance_feature plot

popssf %>% 
  filter(term %in% grep("log_sl|cos_ta", term, value = TRUE,
                        invert = TRUE)) %>% 
  #mutate(term = factor(term, levels = c("dist_aq.ag", "dist_forest", 
  #                                     "dist_road", "dist_settle", "dist_terr.ag", "dist_water"))) %>%
  arrange(term) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5) +
  geom_point(aes(x = term, y = mean, colour = term)) +
  geom_errorbar(aes(x = term, ymin = q025, ymax = q975,
                    colour = term),
                width = 0.25) +
  #facet_wrap(.~id) +
  scale_shape_manual(values = c(3, 16)) +
  #  scale_x_discrete(labels = c("Aq.\nagri.", "For.", "Road", "Settle.",
  #                                "Ter.\nagri.", "Water")
  #  ) +
  labs(x = "Feature", y = expression(beta), colour = "Distance to\nfeature:") +
  scale_colour_manual(values = c("#56B4E9",
                                 "#009E73",
                                 "#000000",
                                 "#E69F00",
                                 "#0072B2")) + #,
  #    labels = c("Aquatic agriculture",
  #               "Forest",
  #               "Road",
  #               "Settlement",
  #               "Terrestrial agriculture",
  #               "Water")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = 4, hjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_line(colour = "grey65", linetype = 1),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = 2),
        legend.text = element_text(lineheight = 1),
        legend.background = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 1, face = 2),
        axis.title.x = element_text(hjust = 0.5, face = 2, margin = margin(10,0,0,0)),
        plot.title = element_text(face = 4),
        strip.text.y = element_blank()
  ) +
  guides(shape = guide_none())


### interaction between step length and feature plot

popssf %>%
  filter(term %in% grep("log_sl", term, value = TRUE,
                        invert = FALSE)) %>%
  arrange(term) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5) +
  geom_errorbar(aes(x = term, ymin = q025, ymax = q975,
                    colour = term),
                width = 0.25) +
  geom_point(aes(x = term, y = mean, colour = term)) +
  #facet_wrap(.~id, scales = "free_y") +
  #scale_shape_manual(values = c(3, 16)) +
  #  scale_x_discrete(labels = c("Aq.\nagri.", "For.", "Road", "Settle.",
  #                              "Ter.\nagri.", "Water")) +
  labs(x = "Feature", y = expression(beta), colour = "Step length\ninteraction with:") +
  scale_colour_manual(values = c("#56B4E9",
                                 "#009E73",
                                 "#000000",
                                 "#E69F00",
                                 "#0072B2")) + #,
  #   labels = c("Aquatic agriculture",
  #               "Forest",
  #              "Road",
  #             "Settlement",
  #             "Terrestrial agriculture",
  #             "Water")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = 4, hjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_line(colour = "grey65", linetype = 1),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = 2),
        legend.text = element_text(lineheight = 1),
        legend.background = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 1, face = 2),
        axis.title.x = element_text(hjust = 0.5, face = 2, margin = margin(10,0,0,0)),
        plot.title = element_text(face = 4),
        strip.text.y = element_blank())



