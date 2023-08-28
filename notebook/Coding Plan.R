
# Coding Plan -------------------------------------------------------------

# Simulate landscapes - store outputs, 3, 1, for each species
#   - only run one however to start
## [FUNC] simulate landscapes

# Simulate individuals ~ 200 ideally that can be sampled from
## [FUNC] simulate individuals

# Sub-sample data tracking regimes (keep more limited than the individual multiverse)
## [FUNC] create tracking regimes, freq and dura

# Create aKDE/KDE/MCP areas (we can drop dBBMMs now as the are conceptually
# wrong)
## [FUNC] generate areas per individual per regime

# some way of randomly selecting a number of those files to create samples
# ready for Pop level area methods - and then complete those calculations
## [FUNC] sample size generation with some way of having multiple repeats of same sample size

# Use the samples to generate type II selection polygons here if possible, so
# we'd have akde, akde_II, kde, kde_II, mcp, mcp_II, then we can generate random
# points in each for the calculations next.
## [FUNC] generate points, get landscape values under points, should result in
# dataframe ready for analysis (name the stored file something linked to the
# sample)

# At this point we should have many compiled files, 1 per indi per frequency per
# duration per area method (which is method x Type [II vs III]), that have:
# - the points used
# - the random points
# - the habitat values at each.
## [FUNC] run method eisera, compana

# SSF can branch from sub-sampling stage easily first randomly select
# individuals for a sample size, then run the SSFs resulting in a single file
# with individual SSF results (no summary at this point as the way the results
# are summarised will be a decision?).
## [FUNC] run ssf, compile sff, summarise ssf

# INLA will need to branch from the regime sub-sampling stage also, but will be
# a single step (possibly with decisions on settings)
## [FUNC] run INLA, compile/summarise INLA
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

