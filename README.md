# evalAdmix

**IMPORTANT**: version 0.95 (updated on 30/06/2021) fixes a bug in the implementation for genotype data, which caused displacement of genotypes between samples when a site had missing data. When all sites have some missingness this would result in the last samples from the analyses having a correlation of nan with all other samples; but might have some more subtle effects whenever there is some level of missingness. If you have analyses based on genotype data with any missingness might be a good idea to re-run them after updating. The bug did not affect the genotype likelihoods implementation so if you based the analyses on genotype likelihoods you do not need to worry. If you applied it to gentoype data without any missingness you also do not need to worry.

evalAdmix allows to evaluate the results of an admixture analysis (i.e. the result of applying [ADMIXTURE](https://genome.cshlp.org/content/19/9/1655.long), [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html), [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) and similar). It only needs the input genotype data used for the previous admixture analysis and the output of that analysis (admixture proportions and ancestral population frequencies). The genotype input data can either be called genotypes in [binary plink format](https://www.cog-genomics.org/plink/1.9/formats#bed) or genotype likelihoods in [beagle format](http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Beagle_format).

The output is a pairwise correlation of residuals matrix between individuals. The correlation will be close to 0 in case of a good fit of the data to the admixture model. When individuals do not fit the model, individuals with similar demographic histories (i.e. usually individuals from the same population) will be positively correlated; and individuals with different histories but that are modelled as sharing one or more ancestral populations as admixture sources will have a negative correlation. Positive correlation between a pair of individuals might also be due to relatedness.

More detailed documentation can be found in [the wiki page](http://www.popgen.dk/software/index.php/EvalAdmix).

# Installation

```
git clone https://github.com/GenisGE/evalAdmix.git
cd evalAdmix
make
```

# Usage

With genotype data

```
./evalAdmix -plink inputPlinkPrefix -fname inputPlinkPrefix.K.P -qname inputPlinkPrefix.K.Q -P 10 -o output.corres.txt

```

With genotype likelihoods

```
./evalAdmix -beagle inputBeagleFile.gz -fname myoutfiles.fopt.gz -qname myoutfiles.qopt -P 10 -o output.corres.txt
```

You can find a tutorial with example input data in [the wiki page](http://www.popgen.dk/software/index.php/EvalAdmix#Run_command_example).

```
Arguments:
	Required:
		-plink path to binary plink file (excluding the .bed)
		or
	      	-beagle path to beagle file containing genotype likelihoods (alternative to -plink)
		
		-fname path to ancestral population frequencies file (space delimited matrix where rows are sites and columns ancestral populations)
	       	-qname path to admixture proportions file (space delimited matrix where rows are individuals and columns ancestral populations)
		
	Optional:       
	
	       -o name of the output file
	       
	 Setup (optional):
	 
	       -P 1 number of threads
	       -autosomeMax 23	 autosome ends with this chromsome (needed only if using genotype input (-plink))
	       -nIts 5	 number of iterations to do for frequency correction; if set to 0 calculates correlation without correction (fast but biased)
	       -useSites 1.0	 proportion of sites to use to calculate correlation of residuals
	       -useInds filename     path to tab delimited file with first column containing all individuals ID and second column with 1 to include individual in analysis and 0 otherwise (default all individuals are included)
	       -misTol 0.05 	 tolerance for considering site as missing when using genotype likelihoods. Use same value as used in NGSadmix to keep compatibility when using genotype likelihoods (-beagle)
	       -minMaf 0.05 	 minimum minor allele frequency to keep site. Use same value as used in NGSadmix to keep compatibility when using genotype likelihoods (-beagle)
```

# Visualization

In R

```
source("visFuns.R")

pop <- as.vector(read.table("inputPlinkPrefix.fam")$V1) # N length character vector with each individual population assignment
q <- as.matrix(read.table("inputPlinkPrefix.K.Q")) # admixture porpotions q is optional for visualization but if used for ordering plot might look better
r <- as.matrix(read.table("output.corres.txt"))

ord <- orderInds(pop=pop, q=q) # ord is optional but this make it easy that admixture and correlation of residuals plots will have individuals in same order

plotAdmix(q=q, pop=pop, ord=ord)
plotCorRes(cor_mat = r, pop = pop, ord=ord, title = "Admixture evaluation as correlation of residuals", max_z=0.25, min_z=-0.25)
```

Usage

```
plotCorRes(cor_mat = cor_mat, pop=NULL, ord=NULL, superpop=NULL,
                       	  title="Correlation of residuals", min_z=NA,max_z=NA, 
                      	  cex.main=1.5, cex.lab=1.5, cex.legend=1.5, color_palette=c("#001260", "#EAEDE9", "#601200"),
                     	  pop_labels = c(T,T), plot_legend = T, adjlab = 0.1, rotatelabpop=0, rotatelabsuperpop=0,lineswidth=1, lineswidthsuperpop=2,
                       adjlabsuperpop=0.16,cex.lab.2 = 1.5)


Arguments

Required
cor_mat	Correlation of residuals matrix.

Optional
pop		Vector with individuals population criteria; will be used to group individuals and for mean correlation per population. Must have same length and order as there are individuals in correlation matrix.
ord		Vector with individual indices specifying order to plot (can be generated with orderInds(q=q, pop=pop) or just order(pop) from base R. 
superpop	Vector with individual superpopulation criteria, to group populations. Must be consisten with pop, only used for labelling.
title		Plot title.
min_z		Minimum value in scale (values below will be plotted as dark blue).
max_z		Maximum value in color scale (values above will be plotted as dark red).
cex.main cex.lab cex.lab.2 cex.legend	   Text size for title, population and superpopulation labels and legend.
color_palette	Vector of length 3 indicating (min, mid, max) color scale.
pop_labels	Logical vector of length 2 indicating whether to write population labels in (y axis, x axis).
plot_legend	Logical indicating whether legend with color scale should be plotted.
adjlab lineswidthsuperpop		Number to adjust distance of population / superpopulation labels to axes.
rotatelabpop rotatelabsuperpop	 Degrees to rotate labels.
lineswidth lineswidthsuperpop	 Width of lines separating populations and superpopulations

```

visFuns.R includes some other functions I find useful (order individuals, order ancestral populations, plot admixture proportions...).



# Citation


[Evaluation of Model Fit of Inferred Admixture Proportions](https://doi.org/10.1111/1755-0998.13171)



# Free pre-print


[Evaluation of Model Fit of Inferred Admixture Proportions](https://doi.org/10.1101/708883)
