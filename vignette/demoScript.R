#load in example binary depletion matrix
data(exampleDepMat)




#load in set of known essential genes
load(BAGEL_essential)
#Generate profiles of number of genes depleted for a number of cell lines from observed data.
pprofile<-ADAM.panessprofile(depMat=depMat)
#Generate a set of random profiles of number of genes depleted for a number of cell lines by perturbing observed data.
nullmodel<-ADAM.generateNullModel(depMat=depMat)
#Calculate log10 odds ratio of observed/expected depletion profiles. Observed values from ADAM.panessprofile and expected as average of random set from ADAM.generateNullModle
EO<-ADAM.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )

#Calculate True positive rates for genes in observed depletion matrix with true positives being in the known set of essential genes.
TPR<-ADAM.truePositiveRate(depMat,BAGEL_essential)
#Calculate minimum number of cell lines a gene needs to be depleted in to be classes as essential
crossoverpoint<-ADAM.tradeoffEO_TPR(EO,TPR$TPR,filename=filename)
#essential genes is the list of genes classed as essential by AdAM.
essentialgenes<-rownames(depMat)[rowSums(depMat)>=crossoverpoint]
