# gwas_simulator
em  
November 6, 2016  



Code for a GWAS simulation

```r
n=1
loci=100
pop = 100

gwas_sim <- function(pop, loci, n, dist){
  #generate allele frequencies from uniform distribution
allele.frequencies = runif(loci, min=0, max=1) 
#effect sizes for each locus, half are 0, half are sampled from a normal distribution
beetas = matrix(c(rnorm(loci/2),rep(0,loci/2)), ncol=1, nrow=loci) 

#make a matrix with a row for each individual, column for each locus, value is the underlying allele frequency of the population
ind.freqs = t(matrix(rep(allele.frequencies, pop), ncol = pop))
#then sample to get the individual genotype
pop.genotypes = apply(ind.freqs,c(1,2), function(x){sum(sample(c(1,0),2,replace=TRUE, prob = c(x,1-x)))} )
#get phenotypes for each individual
pop.phenotypes = pop.genotypes %*% beetas
#calculate the observed allele frequency in the sample
real.freqs = colSums(pop.genotypes)/(pop*2) 


#do the GWAS
doGWAS <- function(x){
  if (real.freqs[x] > 0.05 & real.freqs[x] < 0.95){ 
    out.matrix = summary.lm(lm(pop.phenotypes~pop.genotypes[,x]))$coefficients }
  else { #not enough data to do the GWAS
    out.matrix = matrix(c(rep(NA,7),1), ncol=4)}
  return(out.matrix)
}

gwas.tables = sapply(1:loci, doGWAS)
gwas.pvals = gwas.tables[8,]
gwas.betas = gwas.tables[2,]

get.maf = function(f){
if (f>0.5) {my.maf = 1-f}
else {my.maf = f}
return(my.maf)
}

initial.mafs = sapply(real.freqs, get.maf)
assoc.mafs = initial.mafs[gwas.pvals < 0.05]
causal.mafs.all = initial.mafs[0:50]
causal.mafs = causal.mafs.all[causal.mafs.all > 0.05]

maf.beta.cor = cor.test(initial.mafs,abs(gwas.betas))$estimate

   return(list(causal.mafs,assoc.mafs, maf.beta.cor))
}
```



```r
gwas.sims.200 = sapply(1:200, function(x){gwas_sim(500,100,x)})
```


```r
mean.causal.mafs = sapply(gwas.sims.200[1,], mean)
mean.assoc.mafs = sapply(gwas.sims.200[2,],mean)

plot(1:200, mean.causal.mafs, ylim = c(0,0.5))
points(1:200, mean.assoc.mafs, col = "purple")
legend('bottomright', c('causal mean maf', 'mean maf of gwas hits'),bty="n", col = c("black","purple"), pch=1)
```

![](gwas_simulator_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
hist(mean.assoc.mafs - mean.causal.mafs, col="mediumpurple3", border="white")
```

![](gwas_simulator_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```r
maf.beta.cors = unlist(gwas.sims.200[3,])
hist(maf.beta.cors, col = "darkgreen", border="white", main = "correlations between maf and abs(beta)", xlab = "")
abline(v=0, col = "navy", lwd=2)
```

![](gwas_simulator_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

OK so clearly a skew towards intermediate allele frequencies. There are a few potential issues here. First, I assume a uniform distribution of allele frequencies which is obviously inaccurate. Second, I'm not really getting at how false-positive or false-negative rates vary across allele frequencies with this summary. Third, I chose a p value cut-off of 0.05, which is likely too high.

The previous simulations used a sample size of 500. How much worse is it with a smaller sample size size (100)?

```r
gwas.sims.pop100 = sapply(1:200, function(x){gwas_sim(100,100,x)})
```


```r
mean.causal.mafs = sapply(gwas.sims.pop100[1,], mean)
mean.assoc.mafs = sapply(gwas.sims.pop100[2,],mean)
maf.beta.cors = unlist(gwas.sims.pop100[3,])

plot(1:200, mean.causal.mafs, ylim = c(0,0.5))
points(1:200, mean.assoc.mafs, col = "purple")
legend('bottomright', c('causal mean maf', 'mean maf of gwas hits'),bty="n", col = c("black","purple"), pch=1)
```

![](gwas_simulator_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
hist(mean.assoc.mafs - mean.causal.mafs, col="mediumpurple3", border="white")
abline(v=0, col = "navy", lwd=2)
```

![](gwas_simulator_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
hist(maf.beta.cors, col = "darkgreen", border="white", main = "correlations between maf and abs(beta)", xlab = "")
abline(v=0, col = "navy", lwd=2)
```

![](gwas_simulator_files/figure-html/unnamed-chunk-3-3.png)<!-- -->
