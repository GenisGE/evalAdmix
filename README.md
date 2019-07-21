# Description

evalAdmix allows to evaluate the results of an admixture analysis. It takes as input the genotype data 
(either called genotypes in plink files or genotype likelihoods beagle files) used in the admixture analysis and the frequency  
and admixture propotions (P and Q files) generated.

The output is a pairwise correlation of residuals matrix between individuals The correlation will be 0 in case of a good fit of 
the data to the admixture model. When something is wrong, individuals from the same population will be positively correlated; and 
individuals from different populationts but that share one or more ancestral populations as admixture sources will have a 
negative correlation. Positive correlation between a pair of individuals might also be due to relatedness.

More detailed documentation can be found [here](http://www.popgen.dk/software/index.php/EvalAdmix).

# Installation

```
git clone https://github.com/GenisGE/evalAdmix.git
cd evalAdmix
make
```

# Usage

```
./evalAdmix
```

```
Arguments:
	Required:
		-plink path to binary plink file (excluding the .bed)
		or
	      	-beagle path to beagle file containing genotype likelihoods (alternative to -plink)
		
		-fname path to ancestral population frequencies file
	       	-qname path to admixture proportions file
		
	Optional:       
	
	       -o name of the output file
	       
	 Setup (optional):
	 
	       -P 1 number of threads
	       -autosomeMax 23	 autosome ends with this chromsome
	       -nIts 5	 number of iterations to do for frequency correction; if set to 0 calculates correlation without correction (fast but biased)
	       -useSites 1.0	 proportion of sites to use to calculate correlation of residuals
	       -useInds filename     path to tab delimited file with first column containing all individuals ID and second column with 1 to include individual in analysis and 0 otherwise (default all individuals are included)
	       -misTol 0.05 	 tolerance for considering site as missing when using genotype likelihoods. Use same value as used in NGSadmix to keep compatibility when using genotype likelihoods (-beagle)
	       -minMaf 0.05 	 minimum minor allele frequency to keep site. Use same value as used in NGSadmix to keep compatibility when using genotype likelihoods (-beagle)
```

# Visualization

In R

```
source("NicePlotCorRes.R")

pop <- as.vector(read.table("plink.fam")$V1) # N length character vector with each individual population assignment
r <- as.matrix(read.table("output.corres.txt"))

plotCorRes(cor_mat = r, pop = pop, title = "Admixture evaluation as correlation of residuals", max_z=0.25, min_z=-0.25)

```

# Citation

evalAdmix has a preprint

[Evaluation of Model Fit of Inferred Admixture Proportions](https://doi.org/10.1101/708883)