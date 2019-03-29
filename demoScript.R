#load in example binary depletion matrix
data(exampleDepMat)

# Generate the profiles of number of fitness genes across number of cell lines from observed data and
# corresponding comulative sums.
pprofile<-ADAM.panessprofile(depMat=exampleDepMat)

# Generate a set of random profiles of number of genes depleted for a number of cell lines and corresponding
# cumulative sums by perturbing observed data.
nullmodel<-ADAM.generateNullModel(depMat=exampleDepMat,ntrials = 1000)

#load a reference set of essential genes
data(curated_BAGEL_essential)

# Calculate log10 odd ratios of observed/expected profiles of cumulative number of fitness genes in fixed number of cell lines
# Observed values are from the ADAM.panessprofile function and expected are the average of random set from ADAM.generateNullModle
EO<-ADAM.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )

# Calculate True positive rates for fitness genes in at least n cell lines in the observed dependency matrix,
# with positive cases from a reference set of essential genes
TPR<-ADAM.truePositiveRate(exampleDepMat,curated_BAGEL_essential)


# Calculate minimum number of cell lines a gene needs to be a fitness gene in order to be considered
# as a core-fitness gene
crossoverpoint<-ADAM.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'curated BAGEL essential')

#coreFitnessGenes is the list of genes predicted as core-fitness by AdAM.
coreFitnessGenes<-rownames(exampleDepMat)[rowSums(exampleDepMat)>=crossoverpoint]
