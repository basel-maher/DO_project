#calculate genotype probs and allele probs and calculate kinship

library(qtl2)
load("./results/Rdata/cross_basic.Rdata")

###############
map=cross_basic$gmap

#setting cores to 0 allows for core detection via parallel::detectCores(), but causes a crash if there isnt enough
#storage or memory (storage is used as swap for memory purposes)

#This was done on our high performance cluster
pr <- calc_genoprob(cross_basic, map, cores = 20,error_prob = 0.002,quiet = F)

#cleans genotype probabilities by setting small values to 0 and, for genotype columns where the maximum value is not large, 
#setting all values to 0. This is intended to help with the problem of unstable estimates of genotype effects 
#in scan1coef() and fit1() when there's a genotype that is largely absent.
pr = clean_genoprob(pr)

save(pr,file ="./results/Rdata/pr_basic.Rdata")

apr = genoprob_to_alleleprob(pr)#additive allele model
save(apr,file ="./results/Rdata/apr_basic.Rdata")

#calculate kinship, both Leave One Chromosome Out and overall

k = calc_kinship(apr)
save(k,file ="./results/Rdata/k_basic.Rdata")

k_loco = calc_kinship(apr, type = "loco")
save(k_loco,file ="./results/Rdata/k_loco_basic.Rdata")
