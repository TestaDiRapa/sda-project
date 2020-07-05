library(leaps)
library(car)
library(glmnet)
library(ggplot2)
path_a="/Users/Antonio/Documents/GitHub/sda-project"
path_v="/Users/vincenzo/Documents/GitHub/sda-project"
setwd(path_a)

open_dataset <- function(iso_codes, dates=FALSE) {
  df.data <- read.csv('covid-data.csv')
  if(dates) {
    columns <- colnames(df.data)[c(1,4,9,10,15,16,20,21,22,23,24,25,26,28,29,30,31)]
  } else {
    columns <- colnames(df.data)[c(1,9,10,15,16,20,21,22,23,24,25,26,28,29,30,31)]
  }
  df.data <- df.data[which(df.data$iso_code %in% iso_codes), columns]
  df.data <- na.omit(df.data)
  df.data$new_cases <- df.data$new_cases_per_million
  df.data$new_cases_per_million <- NULL
  df.data$total_cases <- df.data$total_cases_per_million
  df.data$total_cases_per_million <- NULL
  df.data$new_tests <- df.data$new_tests_per_thousand
  df.data$new_tests_per_thousand <- NULL
  df.data$total_tests <- df.data$total_tests_per_thousand
  df.data$total_tests_per_thousand <- NULL
  actual_cases <- c()
  for(code in iso_codes) {
    tmp_actual <- df.data[df.data$iso_code == code,'total_cases']
    actual_cases <- c(actual_cases, 0, tmp_actual[1:length(tmp_actual)-1])
  }
  df.data$actual_cases <- actual_cases
  df.data$iso_code <- NULL
  df.data$total_cases <- NULL
  return(df.data)
}
############################################################################################################
###LOAD THE DATA FROM THE CSV
data_train = open_dataset(c('ITA', 'GBR', 'IND', 'JPN', 'ISR', 'LUX', 'AUS', 'AUT', 'NZL', 'ARG', 'BEL', 'CAN', 'ZAF', 'PRT','ISL','CHE'))
data_test = open_dataset( c('USA','IRN','KOR','URY'))
data_validation = open_dataset(c('RUS','TUR','DNK'))

#CONVERT DATA FROM LINEAR RESPONSE TO LOG RESPONSE
data_train.log <- data_train
data_train.log$new_cases <- log(data_train$new_cases + 0.1)
data_train.log$actual_cases = log(data_train$actual_cases + 0.1)
data_test.log <- data_test
data_test.log$new_cases <- log(data_test$new_cases + 0.1)
data_test.log$actual_cases = log(data_test$actual_cases + 0.1)
data_validation.log <- data_validation
data_validation.log$new_cases <- log(data_validation$new_cases + 0.1)
data_validation.log$actual_cases = log(data_validation$actual_cases + 0.1)
#############################################################################################################


###############################################################################################################
#CHECK DEPENDENCY
dev.new()
pairs(~new_cases+stringency_index+population_density+median_age+population,data_train)
dev.new()
pairs(~new_cases+aged_65_older+aged_70_older+gdp_per_capita+extreme_poverty,data_train)
dev.new()
pairs(~new_cases+cvd_death_rate+diabetes_prevalence+female_smokers+male_smokers,data_train)
dev.new()
pairs(~new_cases+new_tests+total_tests+actual_cases,data_train)
#############################################################################################

########################################################################################################
# FITTING LINEAR MODEL 
fit = lm(new_cases~.,data_train)
summary(fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(fit)

#REMOVE ININFLUENTIAL POINTS
cooksd <- cooks.distance(fit)

# Plot the Cook's Distance using the traditional 4/n criterion
sample_size <- nrow(data_train)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4/sample_size, col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
data_train = data_train[which(!(rownames(data_train) %in% influential)),]
# FITTING LINEAR MODEL WITHOUT ININFLUENT POINTS
fit = lm(new_cases~.,data_train)
summary(fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(fit)



# FITTING LINEAR MODEL (LOG OF RESPONSE)
log.fit = lm(new_cases~.,data_train.log)
summary(log.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(log.fit)

#REMOVE ININFLUENTIAL POINTS
cooksd <- cooks.distance(log.fit)
# Plot the Cook's Distance using the traditional 4/n criterion
sample_size <- nrow(data_train.log)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4/sample_size, col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
data_train.log = data_train.log[which(!(rownames(data_train.log) %in% influential)),]
log.fit = lm(new_cases~.,data_train.log)
summary(log.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(log.fit)
########################################################################################################

#######################################################################################################
#Refit the model with only the significant predictors
adj.fit = lm(new_cases~.-female_smokers-diabetes_prevalence-aged_70_older-cvd_death_rate,data_train)
summary(adj.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(adj.fit)

# Refit the log model with only the significant predictors
adjlog.fit = lm(new_cases~.-female_smokers-actual_cases,data_train.log)
summary(adjlog.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(adjlog.fit)
###############################################################################################################

######################################################################################################################################################################
#BESTSUBSET SELECTION
best.fit=regsubsets(new_cases~.,data_train,nvmax = 15)
reg.summary=summary(best.fit)
#PLOT TO CHOOSE THE BEST NUMBER OF PREDICTORS
dev.new()
par(mfrow=c(2,2))
plot(reg.summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
points(which.min(reg.summary$rss),min(reg.summary$rss), col="red",cex=2,pch=20)
plot(reg.summary$adjr2 ,xlab="Number of Variables ",
     ylab="Adjusted RSq",type="l")
points(which.max(reg.summary$adjr2),max(reg.summary$adjr2), col="red",cex=2,pch=20)
plot(reg.summary$cp ,xlab="Number of Variables ",ylab="Cp", type="l")
points(which.min(reg.summary$cp ),min(reg.summary$cp),col="red",cex=2,pch=20)
plot(reg.summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
points(which.min(reg.summary$bic),min(reg.summary$bic),col="red",cex=2,pch=20)
dev.new()
plot(best.fit,scale="r2")
dev.new()
plot(best.fit,scale="adjr2")
dev.new()
plot(best.fit,scale="Cp")
dev.new()
plot(best.fit,scale="bic")

#CHECK MIN ERROR ON VALIDATION SET
test.mat=model.matrix(new_cases~.,data=data_validation)
# Now we run a loop, and for each size i, we extract the coefficients 
# from regfit.best for the best model of that size, 
# multiply them into the appropriate columns of the test model matrix
# to form the predictions, and compute the test MSE.
val.errors=rep(NA,15)
for(i in 1:15){
  coefi = coef(best.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_validation$new_cases-pred)^2)
}
# The best model is the one that contains which.min(val.errors) (ten in the book) variables.
val.errors; which.min(val.errors) 
coef(regfit.best,which.min(val.errors)) # This is based on training data
# there is no predict() method for regsubsets(). 
# Since we will be using this function again, we can capture our steps above and write our own predict method.
predict.regsubsets = function(object,newdata,id,...){ # ... <-> ellipsis
  form=as.formula(object$call[[2]])
  mat=model.matrix(form, newdata)
  coefi=coef(object, id=id)
  xvars=names(coefi)
  mat[,xvars]%*%coefi
}
# We will demonstrate how we use this function while performing test MSE estimation by cross-validation.

# BACKWARD SELECTION
back.fit=regsubsets(data_train$new_cases~.,data_train,method = "backward",nvmax=15)
back.summary=summary(back.fit)
dev.new()
par(mfrow=c(2,2))
plot(back.summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
points(which.min(back.summary$rss),min(back.summary$rss), col="red",cex=2,pch=20)
plot(back.summary$adjr2 ,xlab="Number of Variables ",
     ylab="Adjusted RSq",type="l")
points(which.max(back.summary$adjr2),max(back.summary$adjr2), col="red",cex=2,pch=20)
plot(back.summary$cp ,xlab="Number of Variables ",ylab="Cp", type="l")
points(which.min(back.summary$cp ),min(back.summary$cp),col="red",cex=2,pch=20)
plot(back.summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
points(which.min(back.summary$bic),min(back.summary$bic),col="red",cex=2,pch=20)
dev.new()
plot(back.fit,scale="r2")
dev.new()
plot(back.fit,scale="adjr2")
dev.new()
plot(back.fit,scale="Cp")
dev.new()
plot(back.fit,scale="bic")

# FORWARD SELECTION
for.fit=regsubsets(data_train$new_cases~.,data_train,method = "forward",nvmax=15)
for.summary=summary(for.fit)
dev.new()
par(mfrow=c(2,2))
plot(for.summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
points(which.min(for.summary$rss),min(for.summary$rss), col="red",cex=2,pch=20)
plot(for.summary$adjr2 ,xlab="Number of Variables ",
     ylab="Adjusted RSq",type="l")
points(which.max(for.summary$adjr2),max(for.summary$adjr2), col="red",cex=2,pch=20)
plot(for.summary$cp ,xlab="Number of Variables ",ylab="Cp", type="l")
points(which.min(for.summary$cp ),min(for.summary$cp),col="red",cex=2,pch=20)
plot(for.summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
points(which.min(for.summary$bic),min(for.summary$bic),col="red",cex=2,pch=20)
dev.new()
plot(for.fit,scale="r2")
dev.new()
plot(for.fit,scale="adjr2")
dev.new()
plot(for.fit,scale="Cp")
dev.new()
plot(for.fit,scale="bic")


#BESTSUBSET SELECTION LOG MODEL
regfit.full=regsubsets(data_train$new_cases_per_million~.-aged_65_older-aged_70_older-gdp_per_capita-extreme_poverty-cvd_death_rate-diabetes_prevalence-female_smokers-male_smokers,data_train,nvmax = 6)
reg.summary=summary(regfit.full)
dev.new()
par(mfrow=c(2,2))
plot(reg.summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
points(which.min(reg.summary$rss),min(reg.summary$rss), col="red",cex=2,pch=20)
plot(reg.summary$adjr2 ,xlab="Number of Variables ",
     ylab="Adjusted RSq",type="l")
points(which.max(reg.summary$adjr2),max(reg.summary$adjr2), col="red",cex=2,pch=20)
plot(reg.summary$cp ,xlab="Number of Variables ",ylab="Cp", type="l")
points(which.min(reg.summary$cp ),min(reg.summary$cp),col="red",cex=2,pch=20)
plot(reg.summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
points(which.min(reg.summary$bic),min(reg.summary$bic),col="red",cex=2,pch=20)
dev.new()
plot(regfit.full,scale="r2")
dev.new()
plot(regfit.full,scale="adjr2")
dev.new()
plot(regfit.full,scale="Cp")
dev.new()
plot(regfit.full,scale="bic")



# BACKWARD SELECTION LOG MODEL
regfit.full=regsubsets(data_train$new_cases~.,data_train,method = "backward")
reg.summary=summary(regfit.full)
dev.new()
par(mfrow=c(2,2))
plot(reg.summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
points(which.min(reg.summary$rss),min(reg.summary$rss), col="red",cex=2,pch=20)
plot(reg.summary$adjr2 ,xlab="Number of Variables ",
     ylab="Adjusted RSq",type="l")
points(which.max(reg.summary$adjr2),max(reg.summary$adjr2), col="red",cex=2,pch=20)
plot(reg.summary$cp ,xlab="Number of Variables ",ylab="Cp", type="l")
points(which.min(reg.summary$cp ),min(reg.summary$cp),col="red",cex=2,pch=20)
plot(reg.summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
points(which.min(reg.summary$bic),min(reg.summary$bic),col="red",cex=2,pch=20)
dev.new()
plot(regfit.full,scale="r2")
dev.new()
plot(regfit.full,scale="adjr2")
dev.new()
plot(regfit.full,scale="Cp")
dev.new()
plot(regfit.full,scale="bic")

# FORWARD SELECTION LOG MODEL
regfit.full=regsubsets(data_train$new_cases_~.,data_train,method = "forward")
reg.summary=summary(regfit.full)
dev.new()
par(mfrow=c(2,2))
plot(reg.summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
points(which.min(reg.summary$rss),min(reg.summary$rss), col="red",cex=2,pch=20)
plot(reg.summary$adjr2 ,xlab="Number of Variables ",
     ylab="Adjusted RSq",type="l")
points(which.max(reg.summary$adjr2),max(reg.summary$adjr2), col="red",cex=2,pch=20)
plot(reg.summary$cp ,xlab="Number of Variables ",ylab="Cp", type="l")
points(which.min(reg.summary$cp ),min(reg.summary$cp),col="red",cex=2,pch=20)
plot(reg.summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
points(which.min(reg.summary$bic),min(reg.summary$bic),col="red",cex=2,pch=20)
dev.new()
plot(regfit.full,scale="r2")
dev.new()
plot(regfit.full,scale="adjr2")
dev.new()
plot(regfit.full,scale="Cp")
dev.new()
plot(regfit.full,scale="bic")

###model with 3 predictors
regfit.3=regsubsets(data_train$new_cases_per_million~.,data_train,nvmax = 3)
lm3 = lm(new_cases_per_million~new_tests_per_thousand+median_age+population, data_train) 
lm3.log = lm(new_cases_per_million~+new_tests_per_thousand+median_age+population, data_train.log)
dev.new()
par(mfrow=c(2,2))
#plot(lm3)
plot(lm3.log)
vif(lm3)
vif(lm3.log)

# #model with 5 predictors per 5 
# regfit.5=regsubsets(data_train$new_cases_per_million~.,data_train,nvmax = 5)
# lm5 = lm(new_cases_per_million~total_tests_per_thousand+new_tests_per_thousand+diabetes_prevalence+stringency_index+median_age, data_train) 
# lm5.log = lm(new_cases_per_million~total_tests_per_thousand+new_tests_per_thousand+diabetes_prevalence+stringency_index+median_age, data_train.log)
# dev.new()
# par(mfrow=c(2,2))
# plot(lm5)
# plot(lm5.log)
# vif(lm5)
# vif(lm5.log)

###model with 4 predictors 
# regfit.4=regsubsets(data_train$new_cases~.,data_train,nvmax = 4)
# lm4 = lm(new_cases~total_tests+cvd_death_rate+male_smokers+new_tests,data_train)
# lm4.log = lm(new_cases~total_tests+cvd_death_rate+male_smokers+new_tests,data_train.log)
# lm4.sqrt = lm(new_cases~total_tests+cvd_death_rate+male_smokers+new_tests,data_train.sqrt)
# dev.new()
# par(mfrow=c(2,2))
# plot(lm4)
# plot(lm4.log)
# plot(lm4.sqrt)
# vif(lm4)
# vif(lm4.log)
# vif(lm4.sqrt)
#model with 7 predictors
# regfit.10=regsubsets(data_train$new_cases~.,data_train,nvmax = 10)
# lm7.log = lm(new_cases~.-stringency_index-population_density-aged_65_older-cvd_death_rate-extreme_poverty-median_age-aged_70_older,data_train.log)
# lm7= lm(new_cases~.-stringency_index-population_density-aged_65_older-cvd_death_rate-extreme_poverty-median_age-aged_70_older,data_train)
# lm7.sqrt = lm(new_cases~.-stringency_index-population_density-aged_65_older-cvd_death_rate-extreme_poverty-median_age-aged_70_older,data_train.sqrt)
# dev.new()
# par(mfrow=c(2,2))
# plot(lm7)
# plot(lm7.log)
# plot(lm7.sqrt)
# vif(lm7)
# vif(lm7.log)
# vif(lm7.sqrt)
# #model with 10 predictors
# lm10.log = lm(new_cases~.-stringency_index-population_density-aged_65_older-cvd_death_rate,data_train.log)
# lm10= lm(new_cases~.-stringency_index-population_density-aged_65_older-cvd_death_rate,data_train)
# lm10.sqrt = lm(new_cases~.-stringency_index-population_density-aged_65_older-cvd_death_rate,data_train.sqrt)
# dev.new()
# par(mfrow=c(2,2))
# plot(lm10)
# plot(lm10.log)
# plot(lm10.sqrt)
# vif(lm10)
# vif(lm10.log)
# vif(lm10.sqrt)


#variable selection with Ridge Regression and the Lasso

# use the glmnet package in order to perform ridge regression and the lasso. We do not use the y ??? x syntax here, but matrix and vector.
# perform ridge regression and the lasso in order to predict Salary on the Hitters data. Missing values has to be removed
# It has an alpha argument that determines what type of model is fit.
# If alpha = 0 a ridge regression model is fit.
# If alpha = 1 (default) then a lasso model is fit. 

x = model.matrix(data_train.log$new_cases_per_million~., data_train.log)[,-1] # without 1's 
y = data_train.log$new_cases_per_million

grid=10^seq(10,-2,length=100)
ridge.mod=glmnet(x,y,alpha=0,lambda=grid)

dim(coef(ridge.mod))




############


ridge.mod$lambda[50] # grid[50] = 11497.57
coef(ridge.mod)[,50] # corresponding coefficients
sqrt(sum(coef(ridge.mod)[-1,50]^2)) # l2 norm
ridge.mod$lambda[60] # lambda = 705.48
coef(ridge.mod)[,60] # corresponding coefficients
sqrt(sum(coef(ridge.mod)[-1,60]^2)) # l2 norm > l2 for lambda[50]

predict(ridge.mod,s=50,type="coefficients")[1:15,] 
# Validation approach to estimate test error
set.seed(1)
train=sample(1:nrow(x), nrow(x)/2) # another typical approach to sample
test=(-train)
y.test=y[test]
# fit a ridge regression model on the training set, and evaluate its MSE on the test set, using lambda = 4. 
ridge.mod=glmnet(x[train,],y[train],alpha=0,lambda=grid,thresh=1e-12)
ridge.pred=predict(ridge.mod,s=4,newx=x[test,]) # Note the use of the predict() function again. This time we get predictions for a test set, by replacing type="coefficients" with the newx argument.
mean((ridge.pred-y.test)^2) # test MSE
mean((mean(y[train ])-y.test)^2) # test MSE, if we had instead simply fit a model with just an intercept, we would have predicted each test observation using the mean of the training observations.
# We could also get the same result by fitting a ridge regression model
# with a very large value of lambda, i.e. 10^10:
ridge.pred=predict(ridge.mod,s=1e10,newx=x[test,])
mean((ridge.pred-y.test)^2) # like intercept only
# Least squares is simply ridge regression with lambda=0;
ridge.pred=predict(ridge.mod,s=0,newx=x[test,],exact=T,x=x[train,],y=y[train]) # corrected according to errata (glmnet pack updated)
# In order for glmnet() to yield the exact least squares coefficients when lambda = 0, we use the argument exact=T when calling the predict() function. Otherwise, the predict() function will interpolate over the grid of lambda values used in fitting the glmnet() model, yielding approximate results. When we use exact=T, there remains a slight discrepancy in the third decimal place between the output of glmnet() when lambda = 0 and the output of lm(); this is due to numerical approximation on the part of glmnet().
mean((ridge.pred-y.test)^2)
# Comapre the results from glmnet when lambda=0 with lm()
lm(y~x, subset=train)
predict(ridge.mod,s=0,exact=T,type="coefficients",x=x[train,],y=y[train])[1:15,] # corrected according to errata 
# In general, if we want to fit a (unpenalized) least squares model, then we should use the lm() function, since that function provides more useful outputs,such as standard errors and p-values.
## CROSS-VALIDATION
# Instead of using the arbitrary value lambda=4, cv.glmnet() uses cross-validation to choose the tuning parameter
# By default, the function performs ten-fold cross-validation, 
# though this can be changed using the argument nfolds.
set.seed (1)
cv.out=cv.glmnet(x[train,],y[train],alpha=0)
dev.new()
plot(cv.out)
bestlam=cv.out$lambda.min; bestlam;log(bestlam) # the best lambda (212 on the text)
cv.out$lambda.1se
ridge.pred=predict(ridge.mod,s=bestlam ,newx=x[test,])
mean((ridge.pred-y.test)^2)

# This represents a further improvement over the test MSE when lambda=4. 
# Finally refit our ridge regression model on the full data set with the best lambda
out=glmnet(x,y,alpha=0)
predict(out,type="coefficients",s=bestlam)[1:15,]
# As expected, none of the coefficients are zero
# ridge regression does not perform variable selection!
dev.new()
plot(out,label = T, xvar = "lambda")


## CROSS-VALIDATION
# Instead of using the arbitrary value lambda=4, cv.glmnet() uses cross-validation to choose the tuning parameter
# By default, the function performs ten-fold cross-validation, 
# though this can be changed using the argument nfolds.
set.seed (1)
cv.out=cv.glmnet(x[train,],y[train],alpha=0)
dev.new()
plot(cv.out)
bestlam=cv.out$lambda.min; bestlam;log(bestlam) # the best lambda (212 on the text)
cv.out$lambda.1se
ridge.pred=predict(ridge.mod,s=bestlam ,newx=x[test,])
mean((ridge.pred-y.test)^2)
# This represents a further improvement over the test MSE when lambda=4. 
# Finally refit our ridge regression model on the full data set with the best lambda
out=glmnet(x,y,alpha=0)
predict(out,type="coefficients",s=bestlam)[1:20,]
# As expected, none of the coefficients are zero
# ridge regression does not perform variable selection!
dev.new()
plot(out,label = T, xvar = "lambda")



test_country <- function(iso_code, model) {
  df.test <- open_dataset(c(iso_code), dates = TRUE)
  print(df.test)
  to.plot <- data.frame(y_real = df.test$new_cases, date = df.test$date)
  df.test$total_cases <- NULL
  df.test$date <- NULL
  to.plot$y_pred <- predict(model, newdata=df.test)
  print(mean((to.plot$y_real-to.plot$y_pred)^2))
  plot <- ggplot(to.plot, aes(x=date)) +
    geom_line(aes(y=y_real, group=1)) +
    geom_line(aes(y=y_pred, group=2))
  print(plot)
}
