#code adapted from www.rstatisticsblog.com/data-science-in-action/lasso-regression
library(glmnet)
library(car)
#library(ppcor)
library(psych)
library(tidyverse)
set.seed(8675309)
#lasso regression to find best predictors of max load
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")


#normalize
norm_pheno = as.data.frame(cross_basic$pheno)

norm_pheno$MAT_VOL1 = norm_pheno$MAT_VOL1 + 1
norm_pheno$MAT_VOL2 = norm_pheno$MAT_VOL2 + 1
norm_pheno$MAT_VOL3 = norm_pheno$MAT_VOL3 + 1
norm_pheno$MAT_VOL4 = norm_pheno$MAT_VOL4 + 1

norm_pheno$bending_work_post_yield = norm_pheno$bending_work_post_yield + 1
norm_pheno$bending_PYD = norm_pheno$bending_PYD + 1

norm_pheno = as.data.frame(log10(norm_pheno[,c(6:14,16,17,21,23:33,34,35,37:41,43:49,51,52,54:58,61:70,72,74,76)]))

pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,36,42,50,53,59,60)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}

##########


#remove glucose
pheno_combined = pheno_combined[,-2]

#cor p vals. adjusted for multiple comparisons (Holm)
x = psych::corr.test(pheno_combined)
#use only phenotypes that are significantly correlated with max load
correlated_phenos = names(which(x$p[,60] <= 0.05))

pheno_combined = pheno_combined[,correlated_phenos]


#correlation matrix
corMat = cor(pheno_combined, use = "p")

#remove body length and gastroc
pheno_combined = pheno_combined[,-c(1,2)]
corMat = cor(pheno_combined, use = "p")

#remove Imax and Imin (correlated with ML and AP)
pheno_combined = pheno_combined[,-c(21,22)]
corMat = cor(pheno_combined, use = "p")

#remove pMOI and ttAR (correlated with ML and AP)
pheno_combined = pheno_combined[,-c(19,20)]
corMat = cor(pheno_combined, use = "p")

#remove work to yield (correlated with yield load)
pheno_combined = pheno_combined[,-c(5)]
corMat = cor(pheno_combined, use = "p")

#remove work_post_yield (corr with total work)
pheno_combined = pheno_combined[,-c(5)]
corMat = cor(pheno_combined, use = "p")

#remove histo bvtv
pheno_combined = pheno_combined[,-c(5)]
corMat = cor(pheno_combined, use = "p")


#remove histo bsbv (correlated negatiely strongly with tbth)
pheno_combined = pheno_combined[,-c(5)]
corMat = cor(pheno_combined, use = "p")

#remove histo tbsp
pheno_combined = pheno_combined[,-c(6)]
corMat = cor(pheno_combined, use = "p")

#remove uCT BMD, conn density, SMI, tbn, and tbsp (cor uct bvtv)
pheno_combined = pheno_combined[,-c(8,9,10,20,12)]
corMat = cor(pheno_combined, use = "p")

#remove uCT bsbv (cor tbth)
pheno_combined = pheno_combined[,-c(15)]
corMat = cor(pheno_combined, use = "p")


#remove bending_yield_stiffness
pheno_combined = pheno_combined[,-c(12)]
corMat = cor(pheno_combined, use = "p")

#remove frax load
pheno_combined = pheno_combined[,-c(13)]
corMat = cor(pheno_combined, use = "p")

#remove ctar/ttar
pheno_combined = pheno_combined[,-c(13)]
corMat = cor(pheno_combined, use = "p")

######


car::vif(lm(bending_max_load~., data = as.data.frame(pheno_combined)))
#ct.ar and ct.th are multicollinear. They are highly correlated. remove the one thats least correlated with max load (uCT_Ct.th)
pheno_combined = pheno_combined[,-13]
car::vif(lm(bending_max_load~., data = as.data.frame(pheno_combined)))

#looks ok now
pheno_combined = as.data.frame(pheno_combined)

#put all predictors in a data matrix
x_vars =  model.matrix(bending_max_load~. , pheno_combined)[,-1]

#make y, pruned to contain only smaples in x
y_var <- pheno_combined[which(rownames(pheno_combined) %in% rownames(x_vars)),"bending_max_load"]

##calculate VIF. over 10 is very significant multicollinearity, over 5 is significant



lambda_seq = 10^seq(2,-3, by=-.1)

#split data into test and train
#80% train

n = 0.8*nrow(x_vars)
train= sample(1:nrow(x_vars), n)

x_test = x_vars[-train,]
x_train = x_vars[train,]

y_test = y_var[-train]
y_train = y_var[train]




#run model for lambdas
#alpha=1 for lasso

#cross validation
#alpha=1 for lasso
cv_output = cv.glmnet(x = x_train, y = y_train, alpha = 1, lambda = lambda_seq)
plot(cv_output)

best_lam= cv_output$lambda.min
#0.2511886


#predict values and compute R2 for the data we trained on, using best lambda
lasso_best = glmnet(x_train, y_train, alpha = 1, lambda = best_lam)

pred = predict(lasso_best, s=best_lam, newx = x_test)
final = cbind(y_test, pred)


##### get R2
actual = final[,1]
preds = final[,2]

rss = sum((preds-actual)^2)
tss = sum((actual - mean(actual))^2)

rmse_test = sqrt(rss/nrow(x_test))
rsq_test <- 1 - rss/tss

rsq_test
rmse_test






#repeat but for training set
pred = predict(lasso_best, s=best_lam, newx = x_train)
final = cbind(y_train, pred)


##### get R2
actual = final[,1]
preds = final[,2]

rss = sum((preds-actual)^2)
tss = sum((actual - mean(actual))^2)

rmse_train = sqrt(rss/nrow(x_train))

rsq_train <- 1 - rss/tss

rsq_train  
rmse_train


coef(lasso_best)
