library(rubias)
library(strataG)
library(tidyverse)
library(ggplot2)
library(Mnov.GTseq.data)
data("ind.info")
load("data/gtypes_final.sams.no.dupes_minReads.20.rda")

strat.scheme <- 'herd'

g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)
g.stratified <- g.stratified[,,
                             filter(getNumInd(g.stratified, by.strata = TRUE), num.ind >= 9) |> 
                               pull(stratum)]

rub.dat <- right_join(
  select(g.stratified@schemes, c(id, wint.area, feed.area, herd)) |> 
    mutate(
      sample_type = ifelse(herd == 'SMx.CA.OR.WA.SBC', 'mixture', 'reference'),
      repunit = ifelse(sample_type == 'reference', wint.area, NA),
      collection = herd
    ) |> 
    left_join(
      select(ind.info, c(id, HAP))
    ) |> 
    mutate(HAP_b = NA),
  g.stratified@data |>
    mutate(new.loc = paste(locus, c(1,2), sep = '_')) |>
    select(-c(locus,stratum)) |>
    pivot_wider(names_from = new.loc, values_from = allele)
) |>
  rename(indiv = id)

# self assignment of the reference samples
herd_assignment <- self_assign(reference = filter(rub.dat, sample_type == 'reference'),
                                   gen_start_col = 8)
# herd_assignment <- self_assign(reference = filter(rub.dat, wint.area %in% c('HI', 'MnMx', 'CentAm'), feed.area %in% c('CA.OR.WA.SBC', 'SEAK.NBC'), sample_type == 'reference'),
#                                gen_start_col = 8)
herd2herd <- left_join(
  select(herd_assignment, c(indiv, collection, repunit, inferred_collection, scaled_likelihood)) |> 
    pivot_wider(names_from = inferred_collection, values_from = scaled_likelihood),
  herd_assignment %>%
    group_by(indiv) %>%
    top_n(1, scaled_likelihood) %>%
    ungroup() |> 
    select(indiv, inferred_collection),
  by = 'indiv')

herd.mat <- table(herd2herd$collection, herd2herd$inferred_collection)
herd.mat
write.csv(herd.mat, file = 'results-raw/rubias_herd2herd_assignment.csv')

herd_to_DPS <- herd_assignment %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood)) |> 
  group_by(indiv) %>%
  top_n(1, repu_scaled_like) %>%
  ungroup()
DPS_mat <- table(herd_to_DPS$collection, herd_to_DPS$inferred_repunit)
DPS_mat
write.csv(DPS_mat, file = 'results-raw/rubias_herd2DPS_assignment.csv')

herd_z <- herd_assignment %>%
  group_by(indiv) %>%
  top_n(1, z_score) %>%
  ungroup()
normo <- tibble(z_score = rnorm(1e06))
ggplot(herd_z, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")

# estimate mixing proportion of the SMex samples
mix_est <- infer_mixture(reference = filter(rub.dat, sample_type == 'reference'),
                         mixture = filter(rub.dat, sample_type == 'mixture'),
                         gen_start_col = 8)

rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit

top2 <- rep_mix_ests %>%
  arrange(desc(repprop)) %>%
  slice(1:2)

# check how many MCMC sweeps were done:
nsweeps <- max(mix_est$mix_prop_traces$sweep)

# discard the first 200 sweeps as burn-in,
# and then aggregate over reporting units
# and then keep only the top2 from above
trace_subset <- mix_est$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top2$repunit)

ggplot(trace_subset, aes(x = repprop, colour = repunit)) +
  geom_density()

