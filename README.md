# Description

evalAdmix allows to evaluate the results of an admixture analysis. It takes as input the genotype data 
(either called genotypes in plink files or genotype likelihoods beagle files) used in the admixture analysis and the frequency  
and admixture propotions (P and Q files) generated.

The output is a pairwise correlation of residuals matrix between individuals The correlation will be 0 in case of a good fit of 
the data to the admixture model. When something is wrong, individuals from the same population will be positively correlated; and 
individuals from different populationts but that share one or more ancestral populations as admixture sources will have a 
negative correlation. Positive correlation between a pair of individuals might also be due to relatedness.

# Installation

# Usage

# Visualization
(maybe add R functions I use to plot?)
