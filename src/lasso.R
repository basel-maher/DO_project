#install.packages("mctest")
#install.packages("ppcor")
#install.packages("psych")
library(mctest)
library(glmnet)

library(ppcor)
library(psych)
#lasso regression to find best predictors of max load
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")


#first get data


norm_pheno = as.data.frame(log10(cross_basic$pheno[,c(6:14,16,17,21,23:33,34,35,37:41,43:49,51,54:58,61:70,72,74,76)]))
pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(5,15,18,19,20,22,36,42,50,52,53,59,60)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}

#remove glucose
pheno_combined = pheno_combined[,-2]

#cor p vals. adjusted for multiple comparisons (Holm)
x = psych::corr.test(pheno_combined)
#use only phenotypes that are significantly correlated with max load
correlated_phenos = names(which(x$p[,60] <= 0.05))

pheno_combined = pheno_combined[,correlated_phenos]

#find multicollinear terms
#Farrar-Glauber test
#explanatory = pheno_combined[,-60]#remove all but max load
#response = pheno_combined[,60] #max load

#overall diagnostics
#od<-omcdiag(x=explanatory,y=response)

#test for multi
#id<-imcdiag(x=explanatory,y=response)

#correlation matrix
corMat = cor(pheno_combined, use = "p")

#remove body length, weight and gastroc
pheno_combined = pheno_combined[,-c(1,2,25)]
corMat = cor(pheno_combined, use = "p")

#remove Imax and Imin (correlated with ML and AP)
pheno_combined = pheno_combined[,-c(21,22)]
corMat = cor(pheno_combined, use = "p")

#remove pMOI and ttAR (correlated with ML and AP)
pheno_combined = pheno_combined[,-c(19,20)]
corMat = cor(pheno_combined, use = "p")

#remove work to yield (correlated with yield load)
pheno_combined = pheno_combined[,-c(6)]
corMat = cor(pheno_combined, use = "p")

#remove work_post_yield (corr with total work)
pheno_combined = pheno_combined[,-c(6)]
corMat = cor(pheno_combined, use = "p")

#remove histo bvtv
pheno_combined = pheno_combined[,-c(6)]
corMat = cor(pheno_combined, use = "p")


#remove histo bsbv (correlated negatiely strongly with tbth)
pheno_combined = pheno_combined[,-c(6)]
corMat = cor(pheno_combined, use = "p")

#remove histo tbsp
pheno_combined = pheno_combined[,-c(7)]
corMat = cor(pheno_combined, use = "p")

#remove uCT BMD, conn density, SMI, tbn, and tbsp (cor uct bvtv)
pheno_combined = pheno_combined[,-c(9,10,12,20,21)]
corMat = cor(pheno_combined, use = "p")

#remove uCT bsbv (cor tbth)
pheno_combined = pheno_combined[,-c(16)]
corMat = cor(pheno_combined, use = "p")


#remove bending_yield_stiffness
pheno_combined = pheno_combined[,-c(13)]
corMat = cor(pheno_combined, use = "p")

#remove frax load
pheno_combined = pheno_combined[,-c(14)]
corMat = cor(pheno_combined, use = "p")

#remove ctar/ttar
pheno_combined = pheno_combined[,-c(14)]
corMat = cor(pheno_combined, use = "p")





# #remove all non-bone traits
# pheno_combined = pheno_combined[,-c(1:8, 56)]
# corMat = cor(pheno_combined, use = "p")
# 
# 
# #remove ML and AP, highly correlated (~0.8 with uCT params)
# pheno_combined = pheno_combined[,-c(1,2)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove BMD (bvtv is better)
# pheno_combined = pheno_combined[,-c(26)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove both disp at frax and frax load
# pheno_combined = pheno_combined[,-c(4,49)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove bending yield load
# pheno_combined = pheno_combined[,-c(1)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove total work, multi with work post yield
# pheno_combined = pheno_combined[,-c(4)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove histo.bvtv
# pheno_combined = pheno_combined[,-c(6)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove histo.ovtv, ovbv
# pheno_combined = pheno_combined[,-c(6,7)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove bsbv
# pheno_combined = pheno_combined[,-c(6)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove histo tbth
# pheno_combined = pheno_combined[,-c(9)]
# corMat = cor(pheno_combined, use = "p")
# 
# 
# #remove histo osbs
# pheno_combined = pheno_combined[,-c(6)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove nobs and nocs, except nob tar and noctar
# pheno_combined = pheno_combined[,-c(10,11,13,41)]
# corMat = cor(pheno_combined, use = "p")
# 
# 
# 
# #remove obsos and histo tbn
# pheno_combined = pheno_combined[,-c(37,12)]
# corMat = cor(pheno_combined, use = "p")
# 
# 
# #remove SMI, tbsp, tbn, and conn_density (all uCT, multi with bvtv)
# pheno_combined = pheno_combined[,-c(37,15,38,13)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove yield stiffness, ctar/ttar
# pheno_combined = pheno_combined[,-c(35,32)]
# corMat = cor(pheno_combined, use = "p")
# 
# 
# #bsbv
# pheno_combined = pheno_combined[,-c(33)]
# corMat = cor(pheno_combined, use = "p")
# 
# #imin, imax, ttar, pmoi
# pheno_combined = pheno_combined[,-c(19:21, 16)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove nonzero MAT
# pheno_combined = pheno_combined[,-c(22:25)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove MAT vol 3
# pheno_combined = pheno_combined[,-c(20)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove bending disp
# pheno_combined = pheno_combined[,-c(2,3)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove obsbs and ocsbs
# pheno_combined = pheno_combined[,-c(4,5)]
# corMat = cor(pheno_combined, use = "p")
# 
# #remove ct.th as it is correlated with ct.ar, but less correlated with max load than ct.ar
# pheno_combined = pheno_combined[,-c(20)]
# corMat = cor(pheno_combined, use = "p")



explanatory = pheno_combined[,-19]#remove all but max load
response = pheno_combined[,19] #max load

#overall diagnostics
od<-omcdiag(x=explanatory,y=response)

#test for multi
id<-imcdiag(x=explanatory,y=response)

####
pheno_combined = as.data.frame(pheno_combined)
x_vars <- model.matrix(bending_max_load~. , pheno_combined)[,-1]

pheno_combined[which(rownames(pheno_combined) %in% rownames(x_vars)),]
y_var <- pheno_combined[which(rownames(pheno_combined) %in% rownames(x_vars)),"bending_max_load"]
lambda_seq <- 10^seq(2, -2, by = -.1)

# Splitting the data into test and train
set.seed(8675309)
train = sample(1:nrow(x_vars), nrow(x_vars)/2)
x_test = x_vars[-train,]
y_test = y_var[-train]

cv_output <- cv.glmnet(x_vars[train,], y_var[train], 
                       alpha = 1, lambda = lambda_seq)

# identifying best lamda
best_lam <- cv_output$lambda.min
#100

lasso_best <- glmnet(x_vars[train,], y_var[train], alpha = 1, lambda = best_lam,standardize = T,standardize.response = T)
pred <- predict(lasso_best, s = best_lam, newx = x_test)

final <- cbind(y_test, pred)
#####
actual <- y_test
preds <- pred
rss <- sum((preds - actual) ^ 2)
tss <- sum((actual - mean(actual)) ^ 2)
rsq <- 1 - rss/tss
rsq


coef(lasso_best)
