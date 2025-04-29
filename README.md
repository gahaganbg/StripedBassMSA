# StripedBassMSA
This is a repository for striped bass genomic data to conduct mixed stock analysis, as presented in Gahagan et al. 2025, *Accurate, precise, and efficient estimates of striped bass mixed-stock catch* (in review).

**Data provided currently includes:**
  1. The 233 SNP GTSeq baseline reference panel (Reference_Genotypes_Genepop NEW.txt) and associated metadata.
  2. The 233 SNP GTSeq genotypes for Massachusetts fisheries dependent striped bass samples (n = 977).
  3. The 233 SNP GTSeq genotypes for Connecticut fisheries dependent striped bass samples (n = 188).

After the manuscript is published R code to conduct the analyses presented will be added to the repository.

**Abstract:**
Mixed stock fisheries present a persistent dilemma to fisheries managers as they can lead to issues from overharvest of rarer populations to fundamental mismatches in biological and management units. Rapid advancements in genomics and genetic stock identification have provided powerful solutions when applied, but there has been inconsistent use of methods for assessing the quality of expected results. We developed an adaptable approach to assess the accuracy and precision of a reference panel for estimating the composition of mixed catch and extended this approach to the examination of actual mixed stock samples to determine if simulation results corresponded to empirical patterns. Applied to a newly developed GT-Seq panel for migratory striped bass, simulations demonstrated that accuracy and precision improved with sample size. Potential gains were mostly realized by 300 individuals and were minimally different than the amount of bias expected under a scenario of perfect genetic assignment. Although accuracy was very high (> 97%) across all simulations, it was best in well-mixed samples and worst when a reporting unit contributed almost all or very few fish to a mixture. Precision, as measured by RMSE and 90% CI width, was high and quickly improved as sample size increased. Results from mixture estimates for samples from coastal fisheries aligned with the trends seen in the simulations and further demonstrated the panel can be reliably used for reporting-unit-specific management of mixed stock striped bass fisheries, as long as they are sampled appropriately.

<img src="https://github.com/user-attachments/assets/58c19a3a-65d6-4f37-a7ac-adf8b0414688" alt="" width="500" height="714">

**Figure 1.** The mean estimated proportion minus the true proportion (A), root mean squared error (B), and 90% credible interval size (C) of the DelChes and KenHud reporting units across all combinations of sample sizes and proportions. Point color represents the mean value at each proportion estimated and the points are jittered horizontally for clarity. Solid black points connected by a black line represent the global mean across proportions at each sample size. In panel C, the horizontal dashed gray line represents a 10% threshold for interval size.

<img src="https://github.com/user-attachments/assets/009920df-e2f6-4dd1-8335-c0af3e7a25a6" alt="" width="500" height="500">

**Figure 2.** Two stock GT-Seq panel results as compared to random samples from the binomial distribution across all sample sizes at two combinations of proportions, 0.8 to 0.2 (A) and 0.6 to 0.4 (B), and boxplots of RMSE (C). In A and B, DelChes (teal) and KenHud (orange) means are represented by colored points, medians by black x’s and 90% credible intervals by dark blue lines. For the random binomial samples, which correspond to the result expected under perfect genetic identification of reporting unit for every fish, the means are represented by black points, medians by black x’s, and 90% confidence intervals by light grey lines. In C, boxplots of Delches (teal), KenHud (orange), and binomial (grey) represent the median with a solid black line, the 25th and 75th quartiles at the hinges, and the whiskers are to the maximum value or 1.5 times the interquartile range. Outliers are represented by black dots. 

<img src="https://github.com/user-attachments/assets/e0a0b0fe-f7a4-4656-9487-6d9471121181" alt="Striped bass mixed stock assignment at multiple smaple sizes" width="500" height="500">

**Figure 5.** Estimated reporting unit proportions of mixed origin striped bass collected during 2018 in coastal fisheries off Massachusetts across nine sample sizes. Colored points represent the mean proportions of the DelChes (teal), KenHud (orange), and NCar (purple) reporting units and black bars represent the 90% credible intervals. Colored rectangles represent the proportions spanned by the 90% credible interval of the results based on the maximum sample size of 977.

