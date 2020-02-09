#install.packages("mctest")
#install.packages("ppcor")
install.packages("psych")
library(mctest)
library(ppcor)
library(psych)
library(glmnet)
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


#find multicollinear terms
#Farrar-Glauber test
explanatory = pheno_combined[,-60]#remove all but max load
response = pheno_combined[,60] #max load

#overall diagnostics
od<-omcdiag(x=explanatory,y=response)

#test for multi
id<-imcdiag(x=explanatory,y=response)

#correlation matrix
corMat = cor(pheno_combined, use = "p")

#remove all fat pads except BFP
pheno_combined = pheno_combined[,-c(2,3,5)]
corMat = cor(pheno_combined, use = "p")

#remove ML and AP, highly correlated (~0.8 with uCT params)
pheno_combined = pheno_combined[,-c(6,7)]
corMat = cor(pheno_combined, use = "p")

#remove both disp at frax
pheno_combined = pheno_combined[,-c(9)]
corMat = cor(pheno_combined, use = "p")

#remove bending yield load
pheno_combined = pheno_combined[,-c(6)]
corMat = cor(pheno_combined, use = "p")

#remove total work, multi with work post yield
pheno_combined = pheno_combined[,-c(9)]
corMat = cor(pheno_combined, use = "p")

#remove all histo
pheno_combined = pheno_combined[,-c(11:26)]
corMat = cor(pheno_combined, use = "p")

#remove SMI
pheno_combined = pheno_combined[,-c(41)]
corMat = cor(pheno_combined, use = "p")

#remove Tb.Sp, multi with BV/TV and BMD
pheno_combined = pheno_combined[,-c(15)]
corMat = cor(pheno_combined, use = "p")
#remove BMD
pheno_combined = pheno_combined[,-c(12)]
corMat = cor(pheno_combined, use = "p")

#remove BSBV and rest of histo
pheno_combined = pheno_combined[,-c(36:38)]
corMat = cor(pheno_combined, use = "p")

#remove ttAr and rest of histo
pheno_combined = pheno_combined[,-c(16)]
corMat = cor(pheno_combined, use = "p")

#remove pMOI
pheno_combined = pheno_combined[,-c(18)]
corMat = cor(pheno_combined, use = "p")

#remove frax load and yield stiffness
pheno_combined = pheno_combined[,-c(31,33)]
corMat = cor(pheno_combined, use = "p")


#remove nonzero MATs
pheno_combined = pheno_combined[,-c(24:27)]
corMat = cor(pheno_combined, use = "p")

#remove ctat/ttar
pheno_combined = pheno_combined[,-c(29)]
corMat = cor(pheno_combined, use = "p")

#still some uCT multis but left for now
explanatory = pheno_combined[,-27]#remove all but max load
response = pheno_combined[,27] #max load

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
#take out non bone traits, switch bv/tv for bmd, add some histo parameters (Nob.TAR/Noc.TAR)
#####
actual <- y_test
preds <- pred
rss <- sum((preds - actual) ^ 2)
tss <- sum((actual - mean(actual)) ^ 2)
rsq <- 1 - rss/tss
rsq


coef(lasso_best)
