# evalAdmix

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
	       -autosomeMax 23	 autosome ends with this chromsome (needed only if genotype (plink) input)
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

Usage

```
plotCorRes(cor_mat, pop=rep(" ", nrow(cor_mat)), superpop=NULL,
                       title="Correlation of residuals", min_z=NULL,max_z=NULL, 
                       is.ord=F, cex.main=1.5, cex.lab=1.5, cex.legend=1.5, color_palette=c("#001260", "#EAEDE9", "#601200"),
                       pop_labels = c(T,T), plot_legend = T, adjlab = 0.1, rotatelab=0){



Arguments

Required
cor_mat	Correlation of residuals matrix.

Optional
pop		Vector with individuals population criteria; will be used to group individuals and for mean correlation per population. Must have same length and order as there are individuals in correlation matrix.
superpop	Vector with individual superpopulation criteria, to group populations. Must be consisten with pop, only used for labelling.
title		Plot title.
min_z		Minimum value in scale (values below will be plotted as dark blue).
max_z		Maximum value in color scale (values above will be plotted as dark red).
is.ord	Logical indicating if cor_mat and pop are in desired order. If F individuals are grouped by population alphabetical order.
cex.main cex.lab cex.legend	   Text size for title, population labels and legend.
color_palette	Vector of length 3 indicating (min, mid, max) color scale.
pop_labels	Logical vector of length 2 indicating whether to write population labels in (y axis, x axis).
plot_legend	Logical indicating whether legend with color scale should be plotted.
adjlab		Number to adjust distance of labels to axes.
rotatelab	Degrees to rotate labels.

```


# Citation

evalAdmix has a preprint

[Evaluation of Model Fit of Inferred Admixture Proportions](https://doi.org/10.1101/708883)