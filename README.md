# User Manual for ANCOM-BC v1.0
The current code implements ANCOM-BC in cross-sectional datasets for comparing absolute abundance changes of taxa among different experimental groups. 

## R-package dependencies
The following libraries need to be included for the R code to run:

```r
library(dplyr)
library(nloptr)
```

## Instructions for use

### Data preprocess

#### Usage

* ```feature_table_pre_process(feature.table, meta.data, sample.var, group.var, zero.cut, lib.cut, neg.lb)```

#### Arguments

*	```feature.table```: Data frame or matrix representing observed OTU table with samples in rows and OTUs (or taxa) in columns.
*	```meta.data```: Data frame or matrix of all variables and covariates of interest.
*	```sample.var```: Character. The name of column storing sample IDs.
*	```group.var```: Character. The name of the main variable of interest. ANCOM-BC v1.0 only supports discrete ```group.var``` and aims to compare the change of absolute abundance across different levels of ```group.var```.
*	```zero.cut```: Numeric fraction. Taxa with proportion of zeroes greater than ```zero.cut``` are not included in the analysis.
* ```lib.cut```: Numeric. Samples with library size less than ```lib.cut``` are not included in the analysis.
*	```neg.lb```: Logical. TRUE indicates a taxon would be classified as a structural zero in an experimental group using its asymptotic lower bound.

#### Value

* ```feature.table```: A data frame of pre-processed OTU table.
*	```library.size```: A numeric vector of library sizes after pre-processing.
*	```group.name```: A character vector of levels of ```group.var```.
*	```group.ind```: A numeric vector. Each sample is assigned a number to indicate its group label for better internal process.
*	```structure.zeros```: A matrix consists of 0 and 1s with 1 indicating the taxon is identified as a structural zero in the corresponding group.

### ANCOM-BC main function

#### Usage:

*	```ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero, adj.method, tol.EM, max.iterNum, perNum, alpha)```

#### Arguments:

*	```feature.table```: Data frame or matrix representing the pre-processed OTU table with samples in rows and OTUs (or taxa) in columns. 
*	```grp.name```: A character vector indicating the levels of group. 
*	```grp.ind```: A numeric vector indicating group assignment for each sample. 1 corresponds to the 1st level of ```grp.name```, 2 corresponds to the 2nd level of ```grp.name``` etc.
*	```struc.zero```: A matrix consists of 0 and 1s with 1 indicating the taxon is identified as a structural zero in the corresponding group.
*	```adj.method```: Character. Returns p-values adjusted using the specified method, including ```"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"```.
*	```tol.EM```: Numeric. The iteration convergence tolerance for E-M algorithm.
*	```max.iterNum```: Numeric. The maximum number of iterations for E-M algorithm.
* ```perNum```: Numeric. The maximum number of permutations. It is active only if there are more than 2 groups.
*	```alpha```: Numeric. Level of significance.

#### Value:
*	```feature.table```: Data frame or matrix. Return the input feature.table.
*	```res```: Data frame. The primary result of ANCOM-BC consisting of: 
  * ```mean.difference```: Numeric. The estimated mean difference of absolute abundance between groups in log scale;
  * ```se```: Numeric. The standard error of ```mean.difference```;
  * ```W```: Numeric. ```mean.difference/se```, which is the test statistic of ANCOM-BC.
  * ```p.val```: Numeric. P-value obtained from two-sided asymptotic Z-test using the test statistic ```W```.
  * ```q.val```. Numeric. Q-value obtained by applying adj.method to ```p-val```.
  * ```diff.abn```. Logical. TRUE if the taxon has ```q.val``` less than ```alpha```.
*	```d```: A numeric vector. Estimated sampling fractions in log scale.
*	```mu```: A numeric vector. Estimated log mean absolute abundance for each group.
*	```bias.est```: Numeric. Estimated mean difference of log sampling fractions between groups through E-M algorithm.

## Flowchart of ANCOM-BC

![Flow Chart](/demos/flowchart.png)

## Examples

```r
# Load example data
data(dietswap)
pseq = dietswap
n_taxa = ntaxa(pseq)
n_samp = nsamples(pseq)
# Metadata
meta_data = meta(pseq)
# Taxonomy table
taxonomy = tax_table(pseq)
# Absolute abundances
otu_absolute = abundances(pseq)
```

### Two-group comparison

```r
# Pre-processing
feature.table = otu_absolute; sample.var = "sample"; group.var = "nationality"; 
zero.cut = 0.90; lib.cut = 1000; neg.lb = TRUE
pre.process = feature_table_pre_process(feature.table, meta_data, sample.var, 
                                        group.var, zero.cut, lib.cut, neg.lb)
feature.table = pre.process$feature.table
group.name = pre.process$group.name
group.ind = pre.process$group.ind
struc.zero = pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name = group.name; grp.ind = group.ind; adj.method = "bonferroni"
tol.EM = 1e-5; max.iterNum = 100; perNum = 1000; alpha = 0.05

out = ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero,
               adj.method, tol.EM, max.iterNum, perNum, alpha)
res = cbind(taxon = rownames(out$feature.table), out$res)
write_csv(res, "demo_two_group.csv")
```

### Multi-group comparison

```r
# Pre-processing
feature.table = otu_absolute; sample.var = "sample"; group.var = "bmi_group"; 
zero.cut = 0.90; lib.cut = 1000; neg.lb = TRUE
pre.process = feature_table_pre_process(feature.table, meta_data, sample.var, 
                                        group.var, zero.cut, lib.cut, neg.lb)
feature.table = pre.process$feature.table
group.name = pre.process$group.name
group.ind = pre.process$group.ind
struc.zero = pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name = group.name; grp.ind = group.ind; adj.method = "bonferroni"
tol.EM = 1e-5; max.iterNum = 100; perNum = 1000; alpha = 0.05

out = ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero,
               adj.method, tol.EM, max.iterNum, perNum, alpha)
res = cbind(taxon = rownames(out$feature.table), out$res)
write_csv(res, "demo_multi_group.csv")
```
