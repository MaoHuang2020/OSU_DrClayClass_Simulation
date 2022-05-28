##### Other scripts
# ?select 10 individuals from GS, randomly cross them and then set phenotypic data on them. We then compare the correlation between GEBVs and the true BVs and between GEBVs and the phenotypes
# ?find the best ones based on estimated breeding values _ebv_ and the true breeding values _gv_0


```{r compare GS selection accuracy and the phenotypic selection accuracy, echo=FALSE}
## Select 10 Inbreds based on PS, 
## Select 10 Inbreds based on GS

```


## NOT SURE !!!!! ???? * Next, we can check on the homozygousity

## GS vs PS

## !!!!!!!!! Ended Here
#### Genotypic value versus breeding value  
By default, the additive genetic variance among founders will be 1.  
```{r Compare breeding value to genotypic value}
# Create haplotypes for founder population of outbred individuals
founderHaps <- AlphaSimR::runMacs(nInd=nFounders, nChr=nChr, segSites=segSites)
# Setup the genotype to phenotype mapping
SP <- AlphaSimR::SimParam$new(founderHaps)
meanDD <- 0.4
varDD <- 0.04
# addToREADME(c("Chunk Compare breeding value to genotypic value",
#               paste("Dominance degree mean", meanDD),
#               paste("Dominance degree variance", varDD), "")
#             )
SP$addTraitAD(nQtlPerChr=nQTL, meanDD=meanDD, varDD=varDD)

# Create a new population of founders
founders <- AlphaSimR::newPop(founderHaps, simParam=SP)


```

> Make a plot of the correlation between genotypic value and breeding value at difference levels of dominance deviation.  
> Does the variance in the dominance deviation have an impact on this relationship?  
  
  ### Estimated versus analytical breeding value  
  Estimating the breeding value by a progeny test.  
```{r Estimate breeding value}
# Error variance for phenotypic evaluations
varE <- 1
# Number of progeny for breeding value estimation
nProgeny1 <- 5
nProgeny2 <- 50

# addToREADME(c("Chunk Estimate breeding value",
#               paste("Phenotypic evaluation with varE", varE),
#               paste("Number of progeny for first BV estimate", nProgeny1),
#               paste("Number of progeny for first BV estimate", nProgeny2), ""
#               )
#             )

# Estimate breeding values  
# ind is the individual whose breeding value you want to estimate  
# pop is the population that individual is in  
# nProgeny is the number of progeny for the test  
# varE is the error variance with which phenotypes are evaluated
#      if the genotypic variance is 1 then varE=1 will give h2 = 0.5
estimateBV <- function(ind, pop, nProgeny, varE=1){
  # Set up crossPlan to cross ind to random others nProgeny times
  crossPlan <- cbind(ind, sample(AlphaSimR::nInd(pop), nProgeny, replace=T))
  progeny <- AlphaSimR::makeCross(founders, crossPlan)
  progPheno <- AlphaSimR::setPheno(progeny, varE=varE, onlyPheno=T)
  return(2*mean(progPheno))
}

estimatedBV <- sapply(1:AlphaSimR::nInd(founders), estimateBV, pop=founders, nProgeny=nProgeny1, varE=varE)
# Compare estimated and analytical breeding values
plot(AlphaSimR::bv(founders), estimatedBV, pch=16, xlab="Analytical value", ylab="Estimated value", main=paste("Breeding value estimated from", nProgeny1, "Progeny"))

estimatedBV <- sapply(1:AlphaSimR::nInd(founders), estimateBV, pop=founders, nProgeny=nProgeny2, varE=varE)
# Compare estimated and analytical breeding values
plot(AlphaSimR::bv(founders), estimatedBV, pch=16, xlab="Analytical value", ylab="Estimated value", main=paste("Breeding value estimated from", nProgeny2, "Progeny"))
```

## Inbreeding depression  
Inbreeding depression should increase as the degree of dominance increases. It has to be directional dominance: that the alleles shifting the phenotype in the same direction are always the dominant ones. Here, inbreeding depression is the difference in the genotypic value between an individual and the mean genotypic value of its progeny from selfing.  
```{r Inbreeding depression}
# Setup a new genotype to phenotype mapping
SP3 <- AlphaSimR::SimParam$new(founderHaps)
# Try different values of degree of dominance
meanDD <- c(0.0, 0.4, 0.8)
varDD <- 0.04
# addToREADME(c("Chunk Inbreeding depression",
#               paste("Dominance degree mean", paste(meanDD, collapse=" ")),
#               paste("Dominance degree variance", varDD), "")
#             )
# Make three genotype to phenotype mappings, one for each dominance deviation
SP3$addTraitAD(nQtlPerChr=nQTL, mean=rep(0, 3), var=rep(1, 3), meanDD=meanDD, varDD=varDD)
# Create a new population of founders
founders3 <- AlphaSimR::newPop(founderHaps, simParam=SP3)

# Estimate individual inbreeding depression
# ind is the individual whose inbreeding depression you want to estimate
# pop is the population that individual is in
# nProgeny is the number of selfed progeny for the test
estimateInbDep <- function(ind, pop, nProgeny){
  # Set up crossPlan to self ind nProgeny times
  crossPlan <- matrix(rep(ind, 2*nProgeny), ncol=2)
  progeny <- AlphaSimR::makeCross(pop, crossPlan, simParam=SP3)
  return(AlphaSimR::gv(pop[ind]) - colMeans(AlphaSimR::gv(progeny)))
}

nProgeny <- 50
estInbDepFndr <- t(sapply(1:AlphaSimR::nInd(founders3), estimateInbDep, pop=founders3, nProgeny=nProgeny))
# Standardize by the genotypic standard deviation
sigmaG <- sqrt(diag(AlphaSimR::varG(founders3)))
for (trait in 1:ncol(estInbDepFndr)){
  estInbDepFndr[,trait] <- estInbDepFndr[,trait]/sigmaG[trait]
}
# Make nice plot
boxplot(estInbDepFndr, xaxt="n", xlab="Mean dominance degree", ylab="Inbreeding depression")
axis(side=1, at=1:3, labels=F)
mtext(meanDD, at=1:3, side=1, line=1, cex=1.3)
```

## `records` data structure  
`AlphaSimR` populations only retain phenotypes of the most recent evaluation. In plant breeding, it is common to evaluate the same line more than once, and it makes sense to include all of those phenotypes in analyses.  Here, I propose a simple tibble to retain phenotypic records.  

### Set up records
Number of cycles, how to weight different traits, selection intensity, etc.
```{r Set up records}
nIndStage1 <- 400
nIndStage2 <- 100
varEstage1 <- 4
varEstage2 <- 1

# addToREADME(c("Chunk Set up records",
#               paste("Number of Stage 1 lines", nIndStage1),
#               paste("Number of Stage 2 lines", nIndStage2),
#               paste("Error variance for Stage 1", varEstage1),
#               paste("Error variance for Stage 2", varEstage2), ""
#               )
#             )
# Function to make a simple data structure out of a population
# AlphaSimR doesn't retain varE once you have setPheno, so supply it
makeRecFromPop <- function(pop, varE=1){
  return(dplyr::tibble(id=pop@id, 
                       mother=pop@mother, 
                       father=pop@father, 
                       pheno=AlphaSimR::pheno(pop), 
                       varE=varE
  )
  )
}

# Make an empty set of records
records <- dplyr::tibble()

# The production pipeline starts with a bunch of new lines
exptLines <- AlphaSimR::randCross(founders, nCrosses=nIndStage1)

# Phenotypic evaluation of experimental lines
exptLines <- AlphaSimR::setPheno(exptLines, varE=varEstage1, simParam=SP)
records <- dplyr::bind_rows(records, makeRecFromPop(exptLines, varE=4))

# Select among lines to advance to Stage 2
keep <- sort(order(AlphaSimR::pheno(exptLines), decreasing=T)[1:nIndStage2])

# Phenotypic evaluation of Stage 2 lines
stage2Lines <- exptLines[keep]
stage2Lines <- AlphaSimR::setPheno(stage2Lines, varE=varEstage2, simParam=SP)
records <- dplyr::bind_rows(records, makeRecFromPop(stage2Lines, varE=1))

str(records)
print(paste("Gain from selection", round(mean(AlphaSimR::bv(stage2Lines)) - mean(AlphaSimR::bv(exptLines)), 2)))
print(paste("Gain from selection", round(mean(AlphaSimR::bv(exptLines)[keep]) - mean(AlphaSimR::bv(exptLines)), 2)))
```

> Similar to the correlation between genotypic value and breeding value that you calculated earlier, evaluate the relationship between the gain from selection and different levels of dominance deviation.  
