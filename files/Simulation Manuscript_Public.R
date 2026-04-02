##### Script Notes----
# This script contains the code to conduct the analyses presented in "Determining accurate and precise genetic estimates of
# mixed-stock catch for efficient sampling: An application to striped bass". It relies on two data files that can also be
# downloaded from the StripedBassMSA Github repository:
# 1. "baseline_samples.rda": The baseline samples, before bias corrections, used for simulations and classifying fishery samples
# 2. "fishery_samples.rda": The fishery samples from Massachusetts in 2018 and Connecticut in 2021

# This script is broken down into 4 sections:
# 1. Basic bias and sensitivity analyses
# 2. Perform simulations to determine the power of the panel to make mixture estimates and individual assignments
#    at various levels of sampling.
# 3. Detect the presence of individuals from rare reporting groups within mixtures at various levels of sampling.
# 4. Perform mixture estimates on real-world samples from MA (n=976) and CT (n=188).

# At the end of Sections 2-4 there are code chunks to produce tables or specific values for reports and manuscripts

# Created by B. Gahagan, MADMF on 4_2_26

setwd("") #Set to your environment

library(rubias); #use rubias v 0.3.1
library(tidyverse);
library(viridis);
library(ggpubr)

set.seed((6287612))

##################### SENSITIVITY AND BIAS TESTS -----
#Background Table for panel
GTSeq_panel %>%
  group_by(repunit) %>%
  count()

###first mixture sensitivity simulation
trial<- assess_reference_loo(reference = GTSeq_panel,
                             gen_start_col = 9, 
                             reps = 50, 
                             mixsize = 200,
)
tmp <- trial %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))

tmp$repunit<-ordered(tmp$repunit, levels = c("GoSL", "SHUB", "STJ", "KenHud", "DelChes", "Carolina"))
b1 <- ggplot(tmp, aes(x = true_repprop, y = reprop_posterior_mean, colour = repunit)) + #plots post mean against true mean
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)+
  labs(color="Reporting Unit")
b1

ggplot(tmp, aes(x = repu_n_prop, y = reprop_posterior_mean, colour = repunit)) + #plots post mean against true proporiton of samples
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit) +
  labs(color="Reporting Unit")

tmp %>%
  # filter(repunit == "DelChes") %>%
  group_by(repunit) %>%
  mutate(diff_rep = reprop_posterior_mean  - true_repprop) %>%
  summarize(mean = mean(diff_rep), min= min(diff_rep), max = max(diff_rep))

# 1:1 lines show some overassignment bias for DelChes at full ref panel size.
# Alleviate by strat random sample from collections within the repunit to level (~50) of other repunits

set.seed((6287612))
GTSeq_norm <- GTSeq_panel %>% #Strat random sample approach
  group_by(repunit) %>%
  #filter(1:n() <52) %>%
  slice_sample(n = 52) %>% #random sample of KenHud and DelChes
  ungroup()

trial<- assess_reference_loo(reference = GTSeq_norm,
                             gen_start_col = 9, 
                             reps = 50, 
                             mixsize = 200,
)
tmp <- trial %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))

tmp$repunit<-ordered(tmp$repunit, levels = c("GoSL", "SHUB", "STJ", "KenHud", "DelChes", "Carolina"))

tmp %>%
  # filter(repunit == "DelChes") %>%
  group_by(repunit) %>%
  mutate(diff_rep = reprop_posterior_mean  - true_repprop) %>%
  summarize(mean = mean(diff_rep), min= min(diff_rep), max = max(diff_rep))

###############Figures S1 and S2
b2 <- ggplot(tmp, aes(x = true_repprop, y = reprop_posterior_mean, colour = repunit)) + #plots post mean against true mean
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)+
  labs(color="Reporting Unit")
b2

ggplot(tmp, aes(x = repu_n_prop, y = reprop_posterior_mean, colour = repunit)) + #plots post mean against true proporiton of samples
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit) +
  labs(color="Reporting Unit")

############Looking for bias from "strays" 
reps=200
m=300
samp_unit="individual"
boing <- assess_reference_loo(
  reference = GTSeq_norm,
  gen_start_col = 9, 
  reps = reps, #200
  mixsize = m, #300
  resampling_unit = samp_unit, 
  return_indiv_posteriors = TRUE #sampunit=indiv
)

#mix and indiv 
boing$indiv_posteriors
boing$mixing_proportions
View(boing$indiv_posteriors)
#simulated_repunit is where the indiv or gc come from, repunit is for PofZ to!

#indiv delches
boing$indiv_posteriors %>% filter(simulated_repunit == "DelChes")
View(boing$indiv_posteriors)

#Histogram of PofZs for simulated DelChes fish assigned to DelChes 
boing$indiv_posteriors %>%
  filter(simulated_repunit == "DelChes", repunit == "DelChes") %>% 
  ggplot(aes(x = PofZ)) + geom_histogram() # some indivs that consistently had high probs

#Histogram of PofZs for simulated DelChes fish assigned to KenHud 
boing$indiv_posteriors %>%
  filter(simulated_repunit == "DelChes", repunit == "KenHud") %>%
  ggplot(aes(x = PofZ)) + geom_histogram() # some indivs consistently had high probs


###################### SCENARIO 1: TWO STOCK SIMULATIONS----
# Bias evident from strays, could "sanitize" the ref panel of strays or use "gene copies" as the resampling unit,
# Chose to use gene copies

########. Simulate accuracy and precision across sample sizes and mixing rates ----

samp_unit <-  "gene_copies" #"individual"
reps <- 200
pvals <-seq(0, 1, by = 0.05) #seq(0, 1, by = 0.05) #0.2
mixsize <- c(50, 100, 150, 200, 250, 300, 400, 500, 1000) #300
names(mixsize) <- mixsize
names(pvals) <- pvals

multi_mix <- lapply(mixsize, function(m) { #This takes a long time to run! Save results for future analyses rather than re-run
  lapply(pvals, function(x) {
    kh<-x #.2
    
    arep<-GTSeq_norm %>%
      distinct(repunit) %>%
      mutate(ppn = c(0, #Carolina
                     1-kh, #DelChes
                     0, #GoSL
                     kh, #KenHud
                     0, #SHUB
                     0)) #STJ
    
    assess_reference_loo(
      reference = GTSeq_norm,
      gen_start_col = 9, 
      reps = reps, #200
      mixsize = m, #300?
      alpha_repunit = arep,
      resampling_unit = samp_unit #, return_indiv_posteriors = TRUE
    )
  }) %>%
    bind_rows(.id = "kh_ppn")
}) %>% bind_rows(.id = "mixsize")
#make things that need to be numbers numerics
multi_mix <- multi_mix %>%
  mutate(mixsize = as.numeric(mixsize), kh_ppn = as.numeric(kh_ppn))

# Get summary data from multi-mix for plotting
PlotData <- multi_mix %>%
  filter(repunit %in% c("KenHud", "DelChes")) %>%
  group_by(repunit, kh_ppn, mixsize) %>%
  summarise(
    mean_pi = mean(post_mean_pi),
    median_pi = median(post_mean_pi),
    quant95 = quantile(post_mean_pi, probs = 0.95),
    quant05 = quantile(post_mean_pi, probs = 0.05),
    mse = sum((true_pi - post_mean_pi) ^ 2) / n(),
    sqrt_mse = sqrt(mse),
    stdev = sd(post_mean_pi),
    true_pi = mean(true_pi),
    n = n()  
  ) %>%
  mutate(mixsize = as.numeric(mixsize), kh_ppn = as.numeric(kh_ppn))

bu<-PlotData
PlotData$PerAcc <- #calculate percent accuracy across all sims
  ifelse(PlotData$true_pi-PlotData$mean_pi > 0 , 100-((PlotData$true_pi - PlotData$mean_pi)/PlotData$true_pi*100),
         ifelse(PlotData$true_pi-PlotData$mean_pi < 0, 100+((PlotData$true_pi - PlotData$mean_pi)/PlotData$true_pi*100),
                NA))

PlotData$PerAcc[is.infinite(PlotData$PerAcc)] <- NA
summary(PlotData$PerAcc) 

########. Generate binomial random error mean, median and quantiles ----
#set values to loop over
rsize <-c(50, 100, 150, 200, 250, 300, 400, 500, 1000)
names(rsize) <- rsize


#### Eric method:Expand Grid
fullMixRandomExpec <- expand_grid(
  mixsize = rsize,
  ppns = seq(0, 1, by = 0.05)
) %>%
  mutate(
    binom_variates = map2(.x = mixsize, .y = ppns, .f = function(x, y) rbinom(200, size = x, prob = y) / x),
    q95 = map_dbl(binom_variates, quantile, probs = 0.95),
    q05 = map_dbl(binom_variates, quantile, probs = 0.05),
    mean = map_dbl(binom_variates, mean),
    median = map_dbl(binom_variates, median),
    mse = map2(.x  = binom_variates, .y = ppns, .f = function(x, y) sum((y - x) ^ 2 / 200)),
    rmse = map_dbl(mse, sqrt)
  )

#######.. Plots of simulation results ----
############## Static plot 

repu_cols <- c( "#1b9e77", "#d95f02" ) # Carolina"#7570b3",, DelChes, KenHud, CANADA"#e7298a"

ggplot(filter(PlotData, kh_ppn == 0.2)) +
  geom_hline(yintercept=0.2, linetype= 4, color="darkgrey") + # 0.2 horizontal line
  geom_hline(yintercept=0.8, linetype= 4, color="darkgrey") + # 0.8 horizontal line
  geom_linerange(aes( x  = mixsize - 10, y = mean_pi, ymin = quant05, ymax = quant95), # sim ests errorbars
                 color="darkblue", linewidth=0.65) +
  geom_point(aes(x  = mixsize - 10, y = mean_pi, shape = repunit, color = repunit), # sim ests mean
             size = 2) +
  geom_point(aes(x  = mixsize - 10, y= median_pi), size = 2, shape = 4) + #sim est median
  geom_linerange(data = filter(fullMixRandomExpec, ppns == 0.2), mapping = aes(x  = mixsize + 10, y = mean, ymin = q05, ymax = q95), #random val errorbars
                 color="gray", linewidth=0.65) +
  geom_linerange(data = filter(fullMixRandomExpec, ppns == 0.8), mapping = aes(x  = mixsize + 10, y = mean, ymin = q05, ymax = q95), #random val errorbars
                 color="gray", linewidth=0.65) +
  geom_point(data = filter(fullMixRandomExpec, ppns == 0.2), mapping = aes(x = mixsize + 10, y = median), colour = "black", shape = 4) + #random vals median
  geom_point(data = filter(fullMixRandomExpec, ppns == 0.8), aes(x  = mixsize + 10, y = mean), color = "black", shape = 18) + # sim ests mean
  scale_color_manual(name = "Reporting Unit", values = repu_cols ) + #make legend uniform
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_y_continuous(breaks = seq(0, 1 , 0.2)) + #set y-axis breaks
  labs(title = "80 DelChes:20 KenHud", #labels
       x = "Mixture size",
       y = "Estimated proportion") +
  theme(axis.text.x=element_text(size=10, angle = 45, hjust = 1), # x-axis text formats
        axis.text.y=element_text(size=10), # y-axis text formats
        axis.text=element_text(colour="black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        legend.position="none",
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "lightgray"))

ggsave("sim test.png", p1, width = 4.5, height = 3.5, units = "in", dpi = 300)


############ Looped plot (Figure S3) 
repu_cols <- c( "#1b9e77", "#d95f02" ) # Carolina"#7570b3",, DelChes, KenHud, CANADA"#e7298a"

PlotData$kh_ppn <- round(PlotData$kh_ppn, 2)
for (z in round(seq(0, 1, by = 0.05), 2)){
  Plot1 <- filter(PlotData, (repunit == "DelChes" & near(true_pi, 1 - z)) | (repunit == "KenHud" & near(true_pi, z)))
  Plot2 <- filter(fullMixRandomExpec, near(ppns, z) | near(ppns, 1-z))
  p <- ggplot() +
    geom_hline(yintercept=z, linetype= 4, color="darkgrey") + # 0.2 horizontal line
    geom_hline(yintercept=1-z, linetype= 4, color="darkgrey") + # 0.8 horizontal line
    geom_linerange(data = Plot1, aes(x  = mixsize - 10, y = mean_pi, ymin = quant05, ymax = quant95), # sim ests errorbars
                   color="darkblue", linewidth=0.65) +
    geom_point(data = Plot1, aes(x  = mixsize - 10, y = mean_pi, shape = repunit, color = repunit), # sim ests mean
               size = 2) +
    geom_point(data = Plot1, aes(x  = mixsize - 10, y= median_pi), size = 2, shape = 4) + #sim est median
    geom_linerange(data = Plot2, mapping = aes(x  = mixsize + 10, y = mean, ymin = q05, ymax = q95), #random val errorbars
                   color="darkgray", linewidth=0.65) +
    geom_point(data = Plot2, mapping = aes(x = mixsize + 10, y = median), colour = "black", shape = 4) + #random vals median
    geom_point(data = Plot2, aes(x  = mixsize + 10, y = mean), color = "black", shape = 18) + # sim ests mean
    scale_color_manual(name = "Reporting Unit", values = repu_cols) + #make legend uniform
    scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1 , 0.2)) + #set y-axis breaks
    labs(title = paste(1-z, "DelChes:", z, "KenHud"), #labels
         x = "Mixture size",
         y = "Estimated proportion") +
    theme(axis.text.x=element_text(size=10, angle = 45, hjust = 1), # x-axis text formats
          axis.text.y=element_text(size=10), # y-axis text formats
          axis.text=element_text(colour="black"), 
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.position="none",
          panel.background = element_rect(fill = "grey95",
                                          colour = "black",
                                          linewidth = 0.5, linetype = "solid"))
  print(p)
  
  ggsave(filename=paste("Ches", 1-z, "Hud", z, ".png", sep="_"), plot=last_plot(),
         width = 4.5, height = 3.5, units = "in", dpi = 300)
}


#############. Compare simulation and random data ----
#####.. Combine random expectation and simulation data for comparisons ----

rando <- fullMixRandomExpec %>%
  mutate(range90_bi = q95 - q05, .after = q05) %>%
  mutate(repunit = "Binom", .before = mixsize) %>%
  mutate(diff_mean_bi = mean - ppns, .after = mean) %>%
  mutate(diff_med_bi = median - ppns, .after = median) %>%
  rename(bi_mix=mixsize, bi_ppn = ppns, bi_q95 = q95, bi_q05 = q05, bi_RMSE = rmse) %>%
  select(-binom_variates, -mse)


CompData <- PlotData %>% 
  ungroup() %>%
  mutate(range90 = quant95 - quant05, .after = quant05) %>%
  mutate(diff_mean = mean_pi - true_pi, .after = mean_pi) %>%
  mutate(diff_med =  median_pi - true_pi, .after = median_pi) %>%
  select(-mse, -n)

CompData$true_pi <- round(CompData$true_pi, 2) #make sure all true_pi values are exact for merge

CompData <- merge(CompData, rando,
                  by.x = c('mixsize', 'true_pi'), 
                  by.y=c('bi_mix', 'bi_ppn'), all=TRUE) 

CompData <- select(CompData, mixsize, true_pi, repunit.x, diff_mean, sqrt_mse, stdev, diff_med,
                   range90, range90_bi, diff_mean_bi, mean_pi, bi_RMSE) %>%
  rename(repunit=repunit.x, RMSE = sqrt_mse) %>%
  mutate(range_diff = range90 - range90_bi) %>%
  mutate(sim_rand_mean = diff_mean - diff_mean_bi) %>%
  mutate(CV = RMSE/mean_pi*100) %>% #new
  mutate(diff_RMSE = RMSE - bi_RMSE) %>%
  ungroup() %>% 
  arrange(repunit, true_pi, as.integer(mixsize)) %>%
  mutate(deltaCV = (lag(CV)-CV)/lag(CV) * 100)
#CV stuff
CompData$deltaCV[CompData$mixsize==50]<-NA #removing fake deltaCVs where proportions changed
CompData$CV[CompData$true_pi==0]<-NA #removing CV for near zero values prior to computing mean CV across samp sizes
CompData$deltaCV[CompData$true_pi==0]<-NA #removing CV for near zero values prior to computing mean CV across samp sizes

print(n = Inf,
      filter(CompData, CV < 10) %>%
        group_by(true_pi, mixsize) %>% #repunit, 
        arrange(true_pi, as.integer(mixsize)) %>%
        count())# %>%
#filter(true_pi==0.5)

print(n = Inf,
      filter(CompData, mixsize == 100) %>%
        group_by(true_pi) %>% #repunit, 
        arrange(true_pi, as.integer(mixsize)) %>%
        select(mixsize, true_pi, repunit, CV, deltaCV))# %>%

#Means for plots
CompSums <- CompData %>%
  group_by(mixsize, repunit) %>%
  summarize(diff_mean = mean(diff_mean),
            diff_med = mean(diff_med),
            RMSE = mean(RMSE),
            CV = mean(CV, na.rm=TRUE),
            range90 = mean (range90),
            range_diff = mean(range_diff),
            sim_rand_mean = mean(sim_rand_mean),
            diff_RMSE = mean(diff_RMSE))

#RMSE for boxplots to show difference
cell <- rando %>% #Random binomial data to merge
  select(bi_mix, repunit, bi_RMSE) %>%
  rename(RMSE = bi_RMSE, mixsize = bi_mix)

CompRMSE <- PlotData %>% #Sim data selection and merge
  ungroup() %>%
  select(mixsize, repunit, sqrt_mse) %>%
  rename(RMSE = sqrt_mse) %>%
  bind_rows(cell)

CompRMSE$Type <- ifelse(CompRMSE$repunit == "Binom", "Binomial", "Simulation") #making columns and ordering to plot well
CompRMSE$Type <- ordered(CompRMSE$Type, levels = c("Simulation", "Binomial"))
CompRMSE$repunit <- ordered(CompRMSE$repunit, levels = c("DelChes", "KenHud", "Binomial"))

rm(cell)

###########.. Plots  ----

#### 90% CI Range
ci <- ggplot() +
  geom_hline(yintercept=0.1, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_point(data = CompData, aes(x  = mixsize, y = range90, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8, position = position_jitter(w = 5, h = 0)) +
  geom_point(data=CompSums, aes(x = mixsize, y = range90), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=CompSums, aes(x = mixsize, y = range90), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = "90% CI range") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))
ci
#### RMSE
ggplot() +
  #geom_hline(yintercept=0.1, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_point(data = CompData, aes(x  = mixsize, y = RMSE, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8, position = position_jitter(w = 5, h = 0)) +
  geom_point(data=CompSums, aes(x = mixsize, y = RMSE), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=CompSums, aes(x = mixsize, y = RMSE), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = "RMSE") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))  

#### CV
cv <- ggplot() +
  geom_hline(yintercept=10, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_point(data = CompData, aes(x  = mixsize, y = CV, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8, position = position_jitter(w = 5, h = 0)) +
  geom_point(data=CompSums, aes(x = mixsize, y = CV), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=CompSums, aes(x = mixsize, y = CV), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  #scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 10)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = "CV") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))  

#### Difference between mean proportion and true proportion
resid <- ggplot() +
  #geom_hline(yintercept=0.1, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_point(data = CompData, aes(x  = mixsize, y = diff_mean, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8, position = position_jitter(w = 5, h = 0)) +
  geom_point(data=CompSums, aes(x = mixsize, y = diff_mean), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=CompSums, aes(x = mixsize, y = diff_mean), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = "Estimated mean - true mean") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))  

#### Difference between median proportion and true proportion
ggplot() +
  #geom_hline(yintercept=0.1, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_point(data = CompData, aes(x  = mixsize, y = diff_med, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8) +
  geom_point(data=CompSums, aes(x = mixsize, y = diff_med), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=CompSums, aes(x = mixsize, y = diff_med), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = "Estimated median - true median") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))

#### Difference between accuracy of simulation mean and random error mean
ggplot() +
  geom_hline(yintercept=0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_hline(yintercept=-0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width  
  geom_point(data = CompData, aes(x  = mixsize, y = sim_rand_mean, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8) +
  geom_point(data=CompSums, aes(x = mixsize, y = sim_rand_mean), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=CompSums, aes(x = mixsize, y = sim_rand_mean), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = " sim:true  - binom:true") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))

#### Difference between simulation RMSE and random RMSE 
ggplot() +
  geom_boxplot(data = CompRMSE, aes(x = as.factor(mixsize), y = RMSE, fill = repunit)) +
  scale_fill_manual(values = c("#009E73", "#D55E00", "#CC79A7")) +
  labs(x = "Mixture Size")+
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    legend.position="none")

#### Difference between simulation 90%CI and random 90%CI
ggplot() +
  geom_hline(yintercept=0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_hline(yintercept=-0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_point(data = CompData, aes(x  = mixsize, y = range_diff, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8) +
  geom_point(data=CompSums, aes(x = mixsize, y = range_diff), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=CompSums, aes(x = mixsize, y = range_diff), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = "Simulation 90% CI  - Random 90% CI") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))

####.. Code to generate metrics for results section ----
#Overall accuracy
CompData %>%
  mutate(accuracy = (1-abs(diff_mean)/true_pi)*100) %>%
  filter(is.finite(accuracy)) %>%         
  summarize(mean = mean(accuracy), min = min(accuracy),
            max = max(accuracy), 
            med = median(accuracy), IQR = IQR(accuracy),
            Q25 = quantile(accuracy)[["25%"]],
            Q75 = quantile(accuracy)[["75%"]]) 

#Code to look at reduction in diff_mean, RMSE, and 90% CI over sample sizes
CompData %>%
  ungroup() %>%
  group_by(repunit, mixsize) %>%
  summarize(range = max(range90) - min(range90)) %>% #Swap in diagnostic of choice in this line RMSE range90 diff_mean
  mutate(per_red = (max(range) - range) / (max(range) - min(range)),
         delta_red = per_red - lag(per_red))

#Code to look at diagnostic mean, median, and diff b/t two
CompData %>%
  ungroup() %>%
  group_by(repunit, mixsize) %>%
  summarize(mean = mean(RMSE), #Swap in diagnostic of choice in this line
            median = median(RMSE)) %>%
  mutate(diff = mean - median)

#Code to find largest value of diagnostic at all mix sizes
print(n = Inf,
      CompData %>%
        ungroup() %>%
        group_by(repunit, mixsize) %>%
        select(repunit, true_pi, mixsize, diff_mean) %>% #swap last term as needed
        mutate(absvalue = abs(diff_mean)) %>% 
        filter(absvalue==max(absvalue)) %>%
        arrange(mixsize))

#Code to find smallest value of diagnostic at all mix sizes
print(n = Inf,
      CompData %>%
        ungroup() %>%
        group_by(repunit, mixsize) %>%
        select(repunit, true_pi, mixsize, diff_mean) %>% #swap last term as needed
        filter(diff_mean==min(diff_mean)) %>%
        arrange(mixsize))  

#Sim vs random metric summaries
CompData %>% 
  group_by(repunit, mixsize) %>%
  summarize(mean = mean(diff_RMSE), min = min(diff_RMSE),
            max = max(diff_RMSE), IQR = IQR(diff_RMSE)) #sim_rand_mean, range_diff, diff_RMSE

#Code to look at global mean reduction in mean_diff
CompSums %>%
  ungroup() %>%
  group_by(repunit) %>%
  arrange(repunit,mixsize) %>%
  mutate(avmean = abs(diff_mean),
         maxav = min(CompSums$diff_mean)*-1,
         per_red = (maxav-avmean) / maxav,
         delta_red = per_red - lag(per_red)) %>%
  select(repunit, mixsize, diff_mean, avmean, maxav, per_red, delta_red)

#Code to look at global mean reduction in RMSE, and 90% CI over sample sizes
CompSums %>%
  ungroup() %>%
  group_by(repunit) %>%
  arrange(repunit,mixsize) %>%
  #summarize(range = max(diff_RMSE) - min(diff_RMSE)) %>% #Swap in diagnostic of choice in this line
  mutate(absv = abs(RMSE), #RMSE   range90
         maxav = max(CompSums$RMSE),
         per_red = (maxav-avmean) / maxav,
         delta_red = per_red - lag(per_red)) %>%
  select(repunit, mixsize, diff_mean, RMSE, avmean, maxav, per_red, delta_red) #RMSE  range90

#Code to look at global mean reduction in diff_mean, RMSE, and 90% CI over sample sizes
CompSums %>%
  ungroup() %>%
  group_by(repunit) %>%
  arrange(repunit,mixsize) %>%
  #summarize(range = max(diff_RMSE) - min(diff_RMSE)) %>% #Swap in diagnostic of choice in this line
  mutate(absav = abs(range90),
         maxav = max(absav),
         minav = min(absav),
         tot_red = maxav-minav,
         per_red = (maxav-absav) / tot_red,
         delta_red = per_red - lag(per_red)) %>%
  select(repunit, mixsize, range90, absav, maxav, minav, per_red, delta_red) #RMSE  range90

#################### Thoughts and Questions ----
#Size and recruitment patterns?
#Borrow from both to increase NC or just ches? Do same with Canadian since diff allele freq?
#leads to recs for: 1. mix est; 2. Mix est plus size 3. Mix est. plus size, plus rare?
# the error \around number of individuals is captured in the variation of post_mean_pi

#################### SCENARIO 2: SMALL PROPORTION SIMULATIONS ----

########. Simulate accuracy and precision across sample sizes and mixing rates ----
#NEED TO DO THE RANDOM ERROR TOO???

samp_unit <-  "gene_copies" #"individual"
reps <- 200
pvals <-seq(0,.1,.01) #seq(0, .1, by = 0.005 to each of 2 stocks) #0.2
mixsize <- c(50, 100, 150, 200, 250, 300, 400, 500, 1000) #300
names(mixsize) <- mixsize
names(pvals) <- pvals

multi_mix_sp <- lapply(mixsize, function(m) {
  lapply(pvals, function(x) {
    kh<-x #.2
    
    arep<-GTSeq_norm %>%
      distinct(repunit) %>%
      mutate(ppn = c(x/2, #Carolina
                     0.8-(x/2), #DelChes
                     0, #GoSL
                     0.2-(x/2), #KenHud
                     0, #SHUB
                     x/2)) #STJ
    
    assess_reference_loo(
      reference = GTSeq_norm,
      gen_start_col = 9, 
      reps = reps, #200
      mixsize = m, #300?
      alpha_repunit = arep,
      resampling_unit = samp_unit #, return_indiv_posteriors = TRUE
    )
  }) %>%
    bind_rows(.id = "rare_ppn")
}) %>% bind_rows(.id = "mixsize")

#make things that need to be numbers numerics
multi_mix_sp <- multi_mix_sp %>%
  mutate(mixsize = as.numeric(mixsize), rare_ppn = as.numeric(rare_ppn))

#save(multi_mix, file = "sprinkleMixSim")

#########################. Generate binomial mean, median and quantiles ----
#set values to loop over
rsize <-c(50, 100, 150, 200, 250, 300, 400, 500, 1000)
names(rsize) <- rsize

#### Eric method:Expand Grid
sprinkleRandomExpec <- expand_grid(
  mixsize = rsize,
  ppns = seq(0, 0.05, by = 0.005)
) %>%
  mutate(
    binom_variates = map2(.x = mixsize, .y = ppns, .f = function(x, y) rbinom(200, size = x, prob = y) / x),
    q95 = map_dbl(binom_variates, quantile, probs = 0.95),
    q05 = map_dbl(binom_variates, quantile, probs = 0.05),
    mean = map_dbl(binom_variates, mean),
    median = map_dbl(binom_variates, median), 
    mse = map2(.x  = binom_variates, .y = ppns, .f = function(x, y) sum((y - x) ^ 2 / 200)),
    rmse = map_dbl(mse, sqrt)
  )

# Get summary data from multi-mix for plotting
SprinkData <- multi_mix_sp %>%
  filter(repunit %in% c("Carolina", "STJ")) %>%
  group_by(repunit, rare_ppn, mixsize) %>%
  summarise(
    mean_pi = mean(post_mean_pi),
    median_pi = median(post_mean_pi),
    quant95 = quantile(post_mean_pi, probs = 0.95),
    quant05 = quantile(post_mean_pi, probs = 0.05),
    mse = sum((true_pi - post_mean_pi) ^ 2) / n(),
    sqrt_mse = sqrt(mse),
    stdev = sd(post_mean_pi),
    true_pi = mean(true_pi),
    n = n()  
  ) %>%
  mutate(mixsize = as.numeric(mixsize), rare_ppn = as.numeric(rare_ppn))

bu<-SprinkData
SprinkData$PerAcc <- 
  ifelse(SprinkData$true_pi-SprinkData$mean_pi > 0 , 100-((SprinkData$true_pi - SprinkData$mean_pi)/SprinkData$true_pi*100),
         ifelse(SprinkData$true_pi-SprinkData$mean_pi < 0, 100+((SprinkData$true_pi - SprinkData$mean_pi)/SprinkData$true_pi*100),
                NA))

SprinkData$PerAcc[is.infinite(SprinkData$PerAcc)] <- NA
summary(SprinkData$PerAcc)


#####.. Combine random expectation and simulation data for comparisons ----

SprinkData$repunit[SprinkData$repunit == "Carolina"]<-"NCar" #correcting keypunch error on RID
SprinkData$repunit[SprinkData$repunit == "STJ"]<-"WoLSJR" #correcting keypunch error on RID

rando_sp <- sprinkleRandomExpec %>%
  mutate(range90_bi = q95 - q05, .after = q05) %>%
  mutate(repunit = "Binom", .before = mixsize) %>%
  mutate(diff_mean_bi = mean - ppns, .after = mean) %>%
  mutate(diff_med_bi = median - ppns, .after = median) %>%
  rename(bi_mix=mixsize, bi_ppn = ppns, bi_q95 = q95, bi_q05 = q05, bi_RMSE = rmse) %>%
  select(-binom_variates, -mse)

SprinkComp <- SprinkData %>%
  ungroup() %>%
  mutate(range90 = quant95 - quant05, .after = quant05) %>%
  mutate(diff_mean = mean_pi - true_pi, .after = mean_pi) %>%
  mutate(diff_med = median_pi - true_pi, .after = median_pi) %>%
  select(-mse, -n)

SprinkComp$true_pi <- round(SprinkComp$true_pi, 3) #make sure all true_pi values are exact for merge or fucks up

SprinkComp <- merge(SprinkComp, rando_sp,
                    by.x = c('mixsize', 'true_pi'), 
                    by.y=c('bi_mix', 'bi_ppn'), all=TRUE) 

SprinkComp <- select(SprinkComp, mixsize, true_pi, repunit.x, diff_mean, sqrt_mse, diff_med,
                     range90, range90_bi, diff_mean_bi, mean_pi, bi_RMSE) %>%
  rename(repunit=repunit.x, RMSE = sqrt_mse) %>%
  mutate(range_diff = range90 - range90_bi) %>%
  mutate(sim_rand_mean = diff_mean - diff_mean_bi) %>%
  mutate(CV = RMSE/mean_pi*100) %>% #new
  mutate(diff_RMSE = RMSE - bi_RMSE) %>%
  ungroup() %>% 
  arrange(repunit, true_pi, as.integer(mixsize)) %>%
  mutate(deltaCV = (lag(CV)-CV)/lag(CV) * 100)
#CV stuff
SprinkComp$deltaCV[SprinkComp$mixsize==50]<-NA #removing fake deltaCVs where proportions changed
SprinkComp$CV[SprinkComp$true_pi==0]<-NA #removing CV for near zero values prior to computing mean CV across samp sizes
SprinkComp$deltaCV[SprinkComp$true_pi==0]<-NA #removing CV for near zero values prior to computing mean CV across samp sizes

#Means for plots
SprinkSums <- SprinkComp %>%
  group_by(mixsize, repunit) %>%
  summarize(diff_mean = mean(diff_mean),
            diff_med = mean(diff_med),
            RMSE = mean(RMSE),
            CV = mean(CV, na.rm=TRUE),
            range90 = mean (range90),
            range_diff = mean(range_diff),
            sim_rand_mean = mean(sim_rand_mean),
            diff_RMSE = mean(diff_RMSE))

#RMSE for boxplots to show difference
cell <- rando_sp %>% #Random binomial data to merge
  select(bi_mix, repunit, bi_RMSE) %>%
  rename(RMSE = bi_RMSE, mixsize = bi_mix)

CompRMSE_sp <- SprinkData %>% #Sim data selection and merge
  ungroup() %>%
  select(mixsize, repunit, sqrt_mse) %>%
  rename(RMSE = sqrt_mse) %>%
  bind_rows(cell)

CompRMSE_sp$Type <- ifelse(CompRMSE_sp$repunit == "Binom", "Binomial", "Simulation") #making columns and ordering to plot well
CompRMSE_sp$Type <- ordered(CompRMSE_sp$Type, levels = c("Simulation", "Binomial"))
CompRMSE_sp$repunit <- ordered(CompRMSE_sp$repunit, levels = c("WolSJR", "NCar", "Binom"))

rm(cell)

################. Output for tables and plots ----
scenario1 <- CompData %>%
  select(mixsize, true_pi, repunit, sim_rand_mean, diff_RMSE, range_diff)
write.csv(scenario1, file = "Table 2 data 3_25.csv")

scenario2 <- SprinkComp %>%
  select(mixsize, true_pi, repunit, sim_rand_mean, diff_RMSE, range_diff)
write.csv(scenario2, file = "Table 3 data 3_25.csv")
############.. Looped plot 
repu_cols <- c("#7570b3","#e7298a") # Carolina "#7570b3", DelChes "#1b9e77", KenHud "#d95f02", CANADA "#e7298a"

z <-.02
for (z in seq(0, 0.1, by = 0.01)){
  Plot1 <- filter(SprinkData, near(rare_ppn, z))
  Plot2 <- filter(sprinkleRandomExpec, near(ppns, z/2))
  p <- ggplot() +
    geom_hline(yintercept=z/2, linetype= 4, color="darkgrey") + # 0.2 horizontal line
    #geom_hline(yintercept=1-z, linetype= 4, color="darkgrey") + # 0.8 horizontal line
    geom_linerange(data = Plot1, aes( x  = mixsize - 10, y = mean_pi, ymin = quant05, ymax = quant95), # sim ests errorbars
                   color="darkblue", linewidth=0.65) +
    geom_point(data = Plot1, aes(x  = mixsize - 10, y = mean_pi, shape = repunit, color = repunit), # sim ests mean
               size = 2) +
    geom_point(data = Plot1, aes(x  = mixsize - 10, y= median_pi), size = 2, shape = 4) + #sim est median
    geom_linerange(data = Plot2, mapping = aes(x  = mixsize + 10, y = mean, ymin = q05, ymax = q95), #random val errorbars
                   color="darkgray", linewidth=0.65) +
    geom_point(data = Plot2, mapping = aes(x = mixsize + 10, y = median), colour = "black", shape = 4) + #random vals median
    geom_point(data = Plot2, aes(x  = mixsize + 10, y = mean), color = "black", shape = 18) + # sim ests mean
    scale_color_manual(name = "Reporting Unit", values = repu_cols) + #make legend uniform
    #scale_shape_discrete(name = "Reporting Unit") + #make legend uniform
    scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
    scale_y_continuous(limits = c(-0.001,.105), breaks = seq(0, .1 , 0.02)) + #set y-axis breaks
    facet_grid(repunit~.) +
    labs(title = paste("Minor proportion = ", z/2), #labels
         x = "Mixture size",
         y = "Estimated proportion") +
    theme(axis.text.x=element_text(size=10, angle = 45, hjust = 1), # x-axis text formats
          axis.text.y=element_text(size=10), # y-axis text formats
          axis.text=element_text(colour="black"), 
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.position="none",
          panel.background = element_rect(fill = "grey95",
                                          colour = "black",
                                          linewidth = 0.5, linetype = "solid"))#,
  #panel.grid.major = element_line(size = 0.5, linetype = 'solid',
  #  colour = "lightgray"))#, 
  #  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
  #                                  colour = "lightgray"))
  print(p)
  
  ggsave(filename=paste("Minor proportion = ", z/2,".png"), plot=last_plot(),device = NULL,
         width = 4.5, height = 3.5, units = "in", dpi = 300)
}

#### Difference between accuracy of simulation mean and random error mean
ggplot() +
  geom_hline(yintercept=0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_hline(yintercept=-0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width  
  geom_point(data = SprinkComp, aes(x  = mixsize, y = sim_rand_mean, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8) +
  geom_point(data=SprinkSums, aes(x = mixsize, y = sim_rand_mean), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=SprinkSums, aes(x = mixsize, y = sim_rand_mean), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = " sim:true  - binom:true") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))

#### Difference between simulation RMSE and random RMSE
ggplot() +
  geom_boxplot(data = CompRMSE_sp, aes(x = as.factor(mixsize), y = RMSE, fill = repunit)) +
  scale_fill_manual(values = c("#009E73", "#D55E00", "#999999")) +
  labs(x = "Mixture Size")+
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    legend.position="none")

#### Difference between simulation 90%CI and random 90%CI
ggplot() +
  geom_hline(yintercept=0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_hline(yintercept=-0.005, linetype= 4, color="darkgrey") + # horizontal line for less than 10% gap in 90% CI width 
  geom_point(data = SprinkComp, aes(x  = mixsize, y = range_diff, color = true_pi), # sim ests mean , color = true_pi
             size = 2, alpha = 0.8) +
  geom_point(data=SprinkSums, aes(x = mixsize, y = range_diff), 
             size =2.8, color = "black", shape = 18) +
  geom_line(data=SprinkSums, aes(x = mixsize, y = range_diff), alpha = 0.6) +
  facet_grid(repunit~.) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_color_viridis(name = "    True \n Proportion") +
  labs(x = "Sample size of mixture", y = "Simulation 90% CI  - Random 90% CI") +
  theme(#panel.background = element_rect(fill = "gray"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text=element_text(colour="black"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    #legend.position="none",
    plot.title = element_text(hjust=0.5))

####.. Code to generate metrics for results section ----
#Accuracy
SprinkComp %>%
  mutate(accuracy = (1-abs(diff_mean)/true_pi)*100) %>%
  filter(is.finite(accuracy)) %>%         
  summarize(mean = mean(accuracy), min = min(accuracy),
            max = max(accuracy), 
            med = median(accuracy), IQR = IQR(accuracy),
            Q25 = quantile(accuracy)[["25%"]],
            Q75 = quantile(accuracy)[["75%"]]) 


#Code to look at reduction in diff_mean, RMSE, and 90% CI over sample sizes
SprinkComp %>%
  ungroup() %>%
  group_by(repunit, mixsize) %>%
  summarize(range = max(RMSE) - min(RMSE)) %>% #Swap in diagnostic of choice in this line
  mutate(per_red = (max(range) - range) / (max(range) - min(range)),
         delta_red = per_red - lag(per_red))


#Code to look at reduction in diff_mean, RMSE, and range90 90% CI over sample sizes
SprinkComp %>%
  ungroup() %>%
  group_by(repunit, mixsize) %>%
  summarize(range = max(range90) - min(range90)) %>% #Swap in diagnostic of choice in this line
  mutate(per_red = (max(range) - range) / (max(range) - min(range)),
         delta_red = per_red - lag(per_red))

#Code to look at diagnostic mean, median, and diff b/t two
SprinkComp %>%
  ungroup() %>%
  group_by(repunit, mixsize) %>%
  summarize(mean = mean(range90), #Swap in diagnostic of choice in this line
            median = median(range90)) %>%
  mutate(diff = mean - median)

#Code to find largest value of diagnostic at all mix sizes
print(n = Inf,
      SprinkComp %>%
        ungroup() %>%
        group_by(repunit, mixsize) %>%
        select(repunit, true_pi, mixsize, RMSE) %>% #swap last term as needed
        filter(RMSE==max(RMSE)) %>%
        arrange(mixsize))

#Code to find smallest value of diagnostic at all mix sizes
print(n = Inf,
      SprinkComp %>%
        ungroup() %>%
        group_by(repunit, mixsize) %>%
        select(repunit, true_pi, mixsize, RMSE) %>% #swap last term as needed
        filter(RMSE==min(RMSE)) %>%
        arrange(mixsize))  

#Sim vs random metric summaries sim_rand_mean, diff_RMSE, range_diff
SprinkComp %>% 
  group_by(repunit, mixsize) %>%
  summarize(mean = mean(diff_RMSE), min = min(diff_RMSE),
            max = max(diff_RMSE), IQR = IQR(diff_RMSE)) #sim_rand_mean, range_diff, diff_RMSE

#Code to look at reduction in sim_rand_mean, diff_RMSE, and range_diff 90% CI over sample sizes
SprinkComp %>%
  ungroup() %>%
  group_by(repunit, mixsize) %>%
  summarize(mean = mean(diff_RMSE), min = min(diff_RMSE), max = max(diff_RMSE),
            range = max(diff_RMSE) - min(diff_RMSE)) %>% #Swap in diagnostic of choice in this line
  mutate(per_red = (max(range) - range) / (max(range) - min(range)),
         delta_red = per_red - lag(per_red))


#Code to look at global mean reduction in diff_mean, RMSE, and 90% CI over sample sizes
SprinkSums %>%
  ungroup() %>%
  group_by(repunit) %>%
  arrange(repunit,mixsize) %>%
  #summarize(range = max(diff_RMSE) - min(diff_RMSE)) %>% #Swap in diagnostic of choice in this line
  mutate(absav = abs(range90),
         maxav = max(absav),
         minav = min(absav),
         tot_red = maxav-minav,
         per_red = (maxav-absav) / tot_red,
         delta_red = per_red - lag(per_red)) %>%
  select(repunit, mixsize, range90, absav, maxav, minav, per_red, delta_red) #RMSE  range90

#....................................................................................
save(multi_mix, multi_mix_sp, GTSeq_norm,
     PlotData, CompData, CompSums, CompRMSE, fullMixRandomExpec, rando,
     sprinkleRandomExpec, SprinkData, rando_sp, SprinkSums, SprinkComp, CompRMSE_sp,
     file="simulation data.rda")
#....................................................................................

##################### ESTIMATING MIXED CATCH SAMPLES ---------------------------------------------
####. Estimates across mixture sizes -----
### Create mixture dataframes
#MA
set.seed((6287612)) ##Did not use in MS analyses. Low samp size different inference the same.
mass <- slice_sample(mass, n = 977) #random reorder
ma_50 <- mass %>% #Strat random sample approach
  slice_head(n = 50) 
ma_100 <- mass %>% #Strat random sample approach
  slice_head(n = 100) 
ma_150 <- mass %>% #Strat random sample approach
  slice_head(n = 150) 
ma_200 <- mass %>% #Strat random sample approach
  slice_head(n = 200) 
ma_250 <- mass %>% #Strat random sample approach
  slice_head(n = 250)  
ma_300 <- mass %>% #Strat random sample approach
  slice_head(n = 300)
ma_400 <- mass %>% #Strat random sample approach
  slice_head(n = 400)
ma_500 <- mass %>% #Strat random sample approach
  slice_head(n = 500)
ma_1000 <- mass 

#CT
set.seed((6287612)) #Did not use in MS analyses. Low samp size different inference the same.
ct <- slice_sample(ct, n = 188) #random reorder
ct_50 <- ct %>% #Strat random sample approach
  slice_head(n = 50) 
ct_100 <- ct %>% #Strat random sample approach
  slice_head(n = 100) 
ct_150 <- ct %>% #Strat random sample approach
  slice_head(n = 150) 
ct_200 <- ct

### Generate means and credible intervals for MA samples across sample sizes
dfList <- list(ma_50, ma_100, ma_150, ma_200, ma_250, ma_300, ma_400, ma_500, ma_1000) #for loop feed
ma_mean_ci <- NULL #repository for loop gen data

for (d in 1:length(dfList)){
  SB_panel_results <- infer_mixture(GTSeq_norm,dfList[[d]], gen_start_col = 9) #calculate gen proportions and assignments  dfList[[d]]
  
  #Aggregate collections into reporting units
  # for mixing proportions
  rep_mix_ests <- SB_panel_results$mixing_proportions %>%
    group_by(mixture_collection, repunit) %>%
    summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit
  
  # find the top 3 most abundant:
  top3 <- rep_mix_ests %>%
    arrange(desc(repprop)) %>% 
    slice(1:3) %>%
    arrange(repunit)
  
  #Create data table for curves
  trace_subset <- SB_panel_results$mix_prop_traces %>%
    group_by(mixture_collection) %>%
    filter(sweep > 100) %>%
    group_by(mixture_collection, sweep, repunit) %>%
    summarise(repprop = sum(pi)) %>%
    filter(repunit %in% top3$repunit)
  
  # Credible intervals
  top3_cis <- trace_subset %>%
    group_by(mixture_collection, repunit) %>%
    summarise(loCI = quantile(repprop, probs = 0.05),
              hiCI = quantile(repprop, probs = 0.95),
              stdev = sd(repprop)) %>%
    arrange(repunit)
  
  #Combine and aggregate
  top3 <- bind_cols(top3, top3_cis[,3:5]) 
  top3$trial <- nrow(dfList[[d]]) #dfList[[d]]
  
  ma_mean_ci<-bind_rows(ma_mean_ci,top3)
}

### Generate means and credible intervals for CT samples across sample sizes
dfList <- list(ct_50, ct_100, ct_150, ct_200) #for loop feed
ct_mean_ci <- NULL #repository for loop gen data

for (d in 1:length(dfList)){
  SB_panel_results <- infer_mixture(GTSeq_norm,dfList[[d]], gen_start_col = 9) #calculate gen proportions and assignments  dfList[[d]]
  
  #Aggregate collections into reporting units
  # for mixing proportions
  rep_mix_ests <- SB_panel_results$mixing_proportions %>%
    group_by(mixture_collection, repunit) %>%
    summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit
  
  # find the top 3 most abundant:
  top3 <- rep_mix_ests %>%
    arrange(desc(repprop)) %>% 
    slice(1:3) %>%
    arrange(repunit)
  
  #Create data table for curves
  trace_subset <- SB_panel_results$mix_prop_traces %>%
    group_by(mixture_collection) %>%
    filter(sweep > 100) %>%
    group_by(mixture_collection, sweep, repunit) %>%
    summarise(repprop = sum(pi)) %>%
    filter(repunit %in% top3$repunit)
  
  # Credible intervals
  top3_cis <- trace_subset %>%
    group_by(mixture_collection, repunit) %>%
    summarise(loCI = quantile(repprop, probs = 0.05),
              hiCI = quantile(repprop, probs = 0.95),
              stdev = sd(repprop)) %>%
    arrange(repunit)
  
  #Combine and aggregate
  top3 <- bind_cols(top3, top3_cis[,3:5]) 
  top3$trial <- nrow(dfList[[d]]) #dfList[[d]]
  
  ct_mean_ci<-bind_rows(ct_mean_ci,top3)
}

ct_mean_ci <-ct_mean_ci %>%
  mutate(indivs = repprop*trial) %>%
  filter(indivs >= 1)
print(n=Inf, ct_mean_ci)

ma_mean_ci <- ma_mean_ci %>%
  mutate(indivs = repprop*trial) %>%
  filter(indivs >= 1)
print(n=Inf, ma_mean_ci)


####. Plots -----
repu_cols <- c("#7570b3", "#1b9e77", "#d95f02", "#e7298a") # Carolina, DelChes, KenHud, CANADA
pd <- position_dodge(0.4)

#MA Plot
ggplot() + 
  geom_rect(aes(xmin = 0, ymin =0.00298, xmax = 1100, ymax = 0.0115), 
            fill = "#7570b3", alpha = 0.2) + #rect for Carolina
  geom_rect(aes(xmin = 0, xmax = 1100, ymin = 0.802, ymax = 0.844),
            fill = "#1b9e77", alpha = 0.2) + #rect for DelChes
  geom_rect(aes(xmin = 0, xmax = 1100, ymin = 0.149, ymax = 0.190),
            fill = "#d95f02", alpha = 0.2) + #rect for KenHud
  geom_linerange(data = ma_mean_ci, aes( x  = trial, y = repprop, ymin = loCI, ymax = hiCI), # sim ests errorbars
                 linewidth=0.65, position = position_jitter(width = 10, seed = 123)) +
  geom_point(data = ma_mean_ci, aes(x  = trial, y = repprop, color = repunit), # sim ests mean
             size = 2, position = position_jitter(width = 10, seed = 123)) +
  scale_color_manual(name = "Reporting Unit", values = repu_cols) + #make legend uniform
  #scale_shape_discrete(name = "Reporting Unit") + #make legend uniform
  scale_x_continuous(limits = c(-0, 1100), breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_y_continuous(limits = c(-0.001, 1), breaks = seq(0, 1.0 , 0.2)) + #set y-axis breaks
  labs(x = "Mixture size",
       y = "Estimated proportion") +
  theme(axis.text.x=element_text(size=10, angle = 45, hjust = 1), # x-axis text formats
        axis.text.y=element_text(size=10), # y-axis text formats
        axis.text=element_text(colour="black"), 
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        panel.background = element_rect(fill = "grey95",
                                        colour = "black",
                                        linewidth = 0.5, linetype = "solid"))#,

#CT Plot
ggplot() + #data = ma_mean_ci, aes(group = repunit)
  geom_linerange(data = ct_mean_ci, aes( x  = trial, y = repprop, ymin = loCI, ymax = hiCI), # sim ests errorbars
                 linewidth=0.65, position = position_jitter(width = 10, seed = 123)) +
  geom_point(data = ct_mean_ci, aes(x  = trial, y = repprop, color = repunit), # sim ests mean
             size = 2, position = position_jitter(width = 10, seed = 123)) +
  #scale_color_manual(name = "Reporting Unit", values = repu_cols) + #make legend uniform
  #scale_shape_discrete(name = "Reporting Unit") + #make legend uniform
  scale_x_continuous(limits = c(-0, 250), breaks = c(50, 100, 150, 200, 250, 300, 400, 500, 1000)) + # set x-axis breaks
  scale_y_continuous(limits = c(-0.001, 1), breaks = seq(0, 1.0 , 0.2)) + #set y-axis breaks
  labs(x = "Mixture size",
       y = "Estimated proportion") +
  theme(axis.text.x=element_text(size=10, angle = 45, hjust = 1), # x-axis text formats
        axis.text.y=element_text(size=10), # y-axis text formats
        axis.text=element_text(colour="black"), 
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        panel.background = element_rect(fill = "grey95",
                                        colour = "black",
                                        linewidth = 0.5, linetype = "solid"))#,

#### Add 90% range column
ma_mean_ci$ma_range90 <- ma_mean_ci$hiCI - ma_mean_ci$loCI
ct_mean_ci$ct_range90 <- ct_mean_ci$hiCI - ct_mean_ci$loCI

#Get MA simulation comps 85:15 for potential table
sim_comps1 <- filter(CompData, repunit == "DelChes" & true_pi == 0.85)
sim_comps2 <- filter(CompData, repunit == "KenHud" & true_pi == 0.15)
sim_comps <- bind_rows(sim_comps1, sim_comps2) %>%
  select(mixsize, true_pi, repunit, mean_pi, stdev, range90) %>%
  mutate(type = "Simulation")

sim_comps_tab <- merge(sim_comps, filter(ma_mean_ci, repunit %in% c("DelChes", "KenHud")),
                       by.x = c('mixsize', "repunit"), 
                       by.y=c('trial', "repunit"), all=TRUE) %>%
  rename(stdev_sim = stdev.x, stdev_emp = stdev.y) %>%
  arrange(repunit, mixsize)

write.csv(sim_comps_tab, file = "sim_est comparison.csv")

#Make df for Figure showing emp results and sim 85:15
sim_comps <- select(sim_comps, - true_pi)
ma_comps <- filter(ma_mean_ci, repunit %in% c("DelChes", "KenHud")) %>%
  ungroup() %>%
  select(trial, repunit, repprop, stdev, ma_range90) %>%
  rename(mixsize = trial, mean_pi = repprop, range90 = ma_range90) %>%
  mutate(type = "MA Empirical")
ct_comps <- filter(ct_mean_ci, repunit %in% c("DelChes", "KenHud")) %>%
  ungroup() %>%
  select(trial, repunit, repprop, stdev, ct_range90) %>%
  rename(mixsize = trial, mean_pi = repprop, range90 = ct_range90) %>%
  mutate(type = "CT Empirical")

sim_comps <- bind_rows(sim_comps, ma_comps, ct_comps)  

rm(ma_comps, ct_comps)

####.. Code to examine diagnostics for results section ----

#Code to look at mean and reduction in 90% CI over sample sizes
print(n = Inf, ma_mean_ci %>%
        ungroup() %>%
        group_by(repunit) %>%
        select(-stdev, -mixture_collection) %>%
        arrange(repunit, trial) %>%
        mutate(maxr90 = max(ma_range90), minr90 = min(ma_range90), #Swap in diagnostic of choice in this line
               per_red = (max(ma_range90) - ma_range90) / (max(ma_range90) - min(ma_range90)),
               delta_red = per_red - lag(per_red)))#,
ifelse(repprop > 0.136 & repprop < 0.174, "YES",
       "NO")),
diff_mean = ifelse(repunit == "DelChes", 0.84 - repprop,
                   ifelse(repunit == "KenHud", 0.155 - repprop,
                          "ERROR")))

#Code to look at diagnostic mean, median, and diff b/t two
ma_mean_ci %>%
  ungroup() %>%
  group_by(repunit) %>%
  summarize(mean = mean(repprop), #Swap in diagnostic of choice in this line
            median = median(repprop),
            min = min(repprop),
            max = max(repprop)) %>%
  mutate(diff = mean - median) 

#Code to look at diagnostic mean, median, and diff b/t two
ct_mean_ci %>%
  ungroup() %>%
  group_by(repunit) %>%
  summarize(mean = mean(repprop), #Swap in diagnostic of choice in this line
            median = median(repprop)) %>%
  mutate(diff = mean - median)

ct_mean_ci %>%
  ungroup() %>%
  filter(repunit %in% c("DelChes", "KenHud")) %>%
  group_by(repunit) %>%
  select(-stdev, -mixture_collection) %>%
  arrange(repunit, trial) %>%
  mutate(per_red = (max(ct_range90) - ct_range90) / (max(ct_range90) - min(ct_range90)), #Swap in diagnostic of choice in this line
         delta_red = per_red - lag(per_red))

#Code to find largest value of diagnostic at all mix sizes
print(n = Inf,
      sim_comps %>%
        ungroup() %>%
        group_by(repunit, mixsize) %>%
        select(repunit, mixsize, type, mean_pi, range90) %>% #swap last term as needed
        #filter(RMSE==max(RMSE)) %>%
        arrange(repunit,type,mixsize))

#Code to find smallest value of diagnostic at all mix sizes
print(n = Inf,
      SprinkComp %>%
        ungroup() %>%
        group_by(repunit, mixsize) %>%
        select(repunit, true_pi, mixsize, RMSE) %>% #swap last term as needed
        filter(RMSE==min(RMSE)) %>%
        arrange(mixsize))  

#Sim vs random metric summaries sim_rand_mean, diff_RMSE, range_diff
SprinkComp %>% 
  group_by(repunit, mixsize) %>%
  summarize(mean = mean(range_diff), min = min(range_diff),
            max = max(range_diff), IQR = IQR(range_diff)) #sim_rand_mean, range_diff, diff_RMSE

#Code to look at reduction in mean_diff, diff_RMSE, and range_diff 90% CI over sample sizes
SprinkComp %>%
  ungroup() %>%
  group_by(repunit, mixsize) %>%
  summarize(range = max(range_diff) - min(range_diff)) %>% #Swap in diagnostic of choice in this line
  mutate(per_red = (max(range) - range) / (max(range) - min(range)),
         delta_red = per_red - lag(per_red))

####### .Massachusetts Regional Analysis -----------------------------------------------------------------
GTSeq_area <- mass
table(GTSeq_area$collection,GTSeq_area$Size) #more of that

table(GTSeq_area$collection)
ssizes <- as.data.frame(table(GTSeq_area$collection)) #for min prop filter below
ssizes <- rename(ssizes, mixture_collection = Var1)

SB_panel_results <- infer_mixture(reference = GTSeq_norm, #This will not run if you have rubias 0.4.0 loaded! PB function throws error
                                  mixture = GTSeq_area,
                                  gen_start_col = 9,
                                  method = "PB") #calculate gen proportions and assignments

#Look at heads of all tables
lapply(SB_panel_results, head)

#Aggregate collections into reporting units
# for mixing proportions
rep_mix_ests <- SB_panel_results$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit
#> `summarise()` regrouping output by 'mixture_collection' (override with `.groups` argument)
rep_mix_ests
filter(rep_mix_ests, repprop >0.0005)  

ggplot(filter(rep_mix_ests, repprop >0.0005), aes(x=mixture_collection, y=repprop, fill = repunit),
       color = "black")+
  geom_bar(stat="identity", position=position_dodge(width = .95)) +
  scale_fill_viridis_d(name = "Reporting Unit") +
  labs(title = "Mixture estimates by aggregation",
       x = "Capture location",
       y = "Proportion of sample")

############# Density curves
# find the top 6 most abundant:
top3 <- rep_mix_ests %>%
  left_join(ssizes,by = "mixture_collection") %>%
  mutate(indivs = repprop*Freq) %>%
  #group_by(ordcoll) %>%
  filter(indivs >= 1)
print(n=Inf, top3)


# check how many MCMC sweeps were done:
max(SB_panel_results$mix_prop_traces$sweep) # n=1999
#Create data table for curves
trace_subset <- SB_panel_results$mix_prop_traces %>%
  group_by(mixture_collection) %>%
  filter(sweep > 100) %>%
  group_by(mixture_collection, sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>%
  filter(repunit %in% top3$repunit)

trace_subset <- trace_subset %>%
  mutate(mix_rep = paste(mixture_collection, repunit)) %>%
  filter(mix_rep != "CT Carolina") %>%
  select(-5)

# Plot density curves by agg (need to figure out how to get rid of Carolina in CT)
aggs <- c("NMA", "CC", "BBVS")
for(k in aggs){
  PlotData <- filter(trace_subset, mixture_collection == k)
  print(ggplot(PlotData, aes(x = repprop, colour = repunit)) +
          geom_density() +
          labs(title = paste(k, "Mixture estimate density curves"))
  )
}

# Credible intervals
top3_cis <- trace_subset %>%
  group_by(mixture_collection, repunit) %>%
  summarise(loCI = quantile(repprop, probs = 0.05),
            hiCI = quantile(repprop, probs = 0.95))
top3_cis #look
top3_cis <- filter(top3_cis, loCI > 1.0e-5) #get rid of Carolina x CC
top3_cis
top3 <- bind_cols(top3, top3_cis[,3:4])

rep_indiv_ests <- SB_panel_results$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))
rep_indiv_ests

#Get highest likelihood of pop assignment for each individual
MaxPost <- SB_panel_results$indiv_posteriors %>%
  group_by(indiv) %>%
  filter(PofZ == max(PofZ)) %>%
  arrange(repunit,indiv)

iba_indivs <- MaxPost %>%
  filter(PofZ > .899) %>%
  group_by(mixture_collection, repunit) %>%
  count()

table(GTSeq_area$collection)
ggplot(top3, aes(x=mixture_collection, y=repprop, ymin = loCI , ymax = hiCI, fill = repunit))+
  geom_bar(stat="identity", position=position_dodge(width = .95)) +
  geom_linerange(position=position_dodge(width=0.95), colour="black", alpha=0.9, linewidth=1) +
  scale_fill_viridis_d(name = "Reporting Unit") +
  labs(title = "Mixture estimates by aggregation",
       x = "Capture location",
       y = "Proportion of sample")

########## Comparing mixing proportions between areas 
PopDiffTest <- function(mixres){
  A <- mixres$mix_prop_traces 
  
  skeleton <- expand_grid(
    mixc1 = unique(A$mixture_collection),
    mixc2 = unique(A$mixture_collection),
    sweep = unique(A$sweep),
    repunit = unique(A$repunit)
  ) %>%
    filter(mixc1 != mixc2)
  
  pw_comps <<- A %>%
    select(-collection) %>%
    left_join(skeleton, ., by = join_by(mixc1 == mixture_collection, sweep, repunit)) %>%
    rename(mixc1_pi = pi) %>%
    left_join(A %>% select(-collection), by = join_by(mixc2 == mixture_collection, sweep, repunit)) %>%
    rename(mixc2_pi = pi) %>%
    mutate(diff = mixc2_pi - mixc1_pi)
  
  #visualize differences in posterior prob distributions
  print(pw_comps %>%
          filter(sweep > 200, repunit %in% c("DelChes", "KenHud")) %>%
          ggplot(aes(x = diff)) +
          geom_density(colour = "blue") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          facet_grid(mixc1 ~ mixc2 + repunit))
  
  #df to summarize differences among groupings
  diffs <<- pw_comps %>%
    filter(sweep > 200, repunit %in% c("DelChes", "KenHud")) %>%
    group_by(repunit, mixc1, mixc2) %>%
    summarise(
      diff_gt0 = mean(diff > 0) # asks what is the mean of the occurrences this is > 0. These are the mirror images b/c everything is done 2x
    ) #%>% 
  View(diffs)
  
}
PopDiffTest(SB_panel_results)

