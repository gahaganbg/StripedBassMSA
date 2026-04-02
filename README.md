# StripedBassMSA
This is a repository for striped bass genomic data to conduct mixed stock analysis, as presented in Gahagan et al., *Determining accurate and precise genetic estimates of mixed-stock catch for efficient sampling: An application to striped bass*.

**Data provided currently includes:**
  1. The 233 SNP GTSeq baseline reference panel (Reference_Genotypes_Genepop NEW.txt) and associated metadata, as presented in LeBlanc et al. 2025 (in press).
  2. The 233 SNP GTSeq genotypes for Massachusetts fisheries dependent striped bass samples (n = 977).
  3. The 233 SNP GTSeq genotypes for Connecticut fisheries dependent striped bass samples (n = 188).

Code to perform analyses can be found in the R script *insert name here*

**Abstract:**
Rapid advancements in genomics and genetic stock identification have provided powerful solutions for resolving the dilemmas presented by mixed stock fisheries, but explicit methods to fully test a panel and determine appropriate sample sizes to collect are still needed. We developed a simulation approach to assess the accuracy, precision, and source of error of any genetic reference panel for mixing proportion estimates and applied it to a GT-Seq panel for striped bass. We also compared the simulation results to empirical estimates from mixed stock fishery samples collected in Massachusetts and Connecticut. Accuracy was very high (> 97%) across all simulations, but it was best in well-mixed samples and worst when a reporting unit contributed almost all or very few fish to a mixture. Precision was also high; 90% credible intervals were lower than 0.090 (on a scale of 0.000 to 1.000) and more than 50% of the total improvement had been realized at 200 samples. Fishery mixture estimates typically had better precision than simulations but improved more slowly with sample size, demonstrating that the simulation approach provided valuable guidance when employing new panels. The Delaware-Chesapeake reporting unit contributed 79.0% to Connecticut samples and 82.3% to Massachusetts and was significantly higher in the southernmost portion of that state (89.1%), indicating spatial variability in catch composition. The Kennebec-Hudson reporting unit made up the remainder in Connecticut but not in Massachusetts, where a small proportion of North Carolina origin fish were detected (0.7%, 90% CI = 0.3 – 1.1%). 

<img src="https://github.com/user-attachments/assets/58c19a3a-65d6-4f37-a7ac-adf8b0414688" alt="" width="500" height="714">

**Figure 2.** The mean estimated proportion minus the true proportion (A), root mean squared error (B), and 90% credible interval size (C) of the DelChes and KenHud reporting units across all combinations of sample sizes and proportions. Point color represents the mean value at each proportion estimated and the points are jittered horizontally for clarity. Solid black points connected by a black line represent the global mean across proportions at each sample size. In panel C, the horizontal dashed gray line represents a 10% threshold for interval size.

<img src="https://github.com/user-attachments/assets/009920df-e2f6-4dd1-8335-c0af3e7a25a6" alt="" width="500" height="500">

**Figure 3.** Two stock GT-Seq panel results as compared to random samples from the binomial distribution across all sample sizes at two combinations of proportions, 0.8 to 0.2 (A) and 0.6 to 0.4 (B), and boxplots of RMSE (C). In A and B, DelChes (teal) and KenHud (orange) means are represented by colored points, medians by black x’s and 90% credible intervals by dark blue lines. For the random binomial samples, which correspond to the result expected under perfect genetic identification of reporting unit for every fish, the means are represented by black points, medians by black x’s, and 90% confidence intervals by light grey lines. In C, boxplots of Delches (teal), KenHud (orange), and binomial (grey) represent the median with a solid black line, the 25th and 75th quartiles at the hinges, and the whiskers are to the maximum value or 1.5 times the interquartile range. Outliers are represented by black dots. 

<img src="https://github.com/user-attachments/assets/e0a0b0fe-f7a4-4656-9487-6d9471121181" alt="Striped bass mixed stock assignment at multiple smaple sizes" width="500" height="500">

**Figure 6.** Estimated reporting unit proportions of mixed origin striped bass collected during 2018 in coastal fisheries off Massachusetts across nine sample sizes. Colored points represent the mean proportions of the DelChes (teal), KenHud (orange), and NCar (purple) reporting units and black bars represent the 90% credible intervals. Colored rectangles represent the proportions spanned by the 90% credible interval of the results based on the maximum sample size of 977.

