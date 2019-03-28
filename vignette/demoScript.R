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





#Calculate log10 odds ratio of observed/expected depletion profiles. Observed values from ADAM.panessprofile and expected as average of random set from ADAM.generateNullModle
EO<-ADAM.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )

#Calculate True positive rates for genes in observed depletion matrix with true positives being in the known set of essential genes.
TPR<-ADAM.truePositiveRate(depMat,BAGEL_essential)
#Calculate minimum number of cell lines a gene needs to be depleted in to be classes as essential
crossoverpoint<-ADAM.tradeoffEO_TPR(EO,TPR$TPR,filename=filename)
#essential genes is the list of genes classed as essential by AdAM.
essentialgenes<-rownames(depMat)[rowSums(depMat)>=crossoverpoint]
