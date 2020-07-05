#LOAD LIBRARIES AND DEFINE FUNCTIONS##############################
library(leaps)
library(car)
library(glmnet)
library(ggplot2)
path_a="/Users/Antonio/Documents/GitHub/sda-project"
path_v="/Users/vincenzo/Documents/GitHub/sda-project"
setwd(path_a)
open_dataset <- function(iso_codes, dates=FALSE, iso=FALSE) {
  df.data <- read.csv('covid-data.csv')
  if(dates) {
    columns <- colnames(df.data)[c(1,4,9,10,15,16,20,21,22,23,24,25,26,28,29,30,31)]
  } else {
    columns <- colnames(df.data)[c(1,9,10,15,16,20,21,22,23,24,25,26,28,29,30,31)]
  }
  df.data <- df.data[which(df.data$iso_code %in% iso_codes), columns]
  df.data <- na.omit(df.data)
  df.data$total_cases <- df.data$total_cases_per_million
  df.data$total_cases_per_million <- NULL
  df.data$new_cases <- df.data$new_cases_per_million
  df.data$new_cases_per_million <- NULL
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
  if(!iso) {
    df.data$iso_code <- NULL
  }
  return(df.data)
}
plot_best_predictors = function(summary,fit){
  dev.new()
  par(mfrow=c(2,2))
  plot(summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
  points(which.min(summary$rss),min(summary$rss), col="red",cex=2,pch=20)
  plot(summary$adjr2 ,xlab="Number of Variables ",
       ylab="Adjusted RSq",type="l")
  points(which.max(summary$adjr2),max(summary$adjr2), col="red",cex=2,pch=20)
  plot(summary$cp ,xlab="Number of Variables ",ylab="Cp", type="l")
  points(which.min(summary$cp ),min(summary$cp),col="red",cex=2,pch=20)
  plot(summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
  points(which.min(summary$bic),min(summary$bic),col="red",cex=2,pch=20)
  dev.new()
  plot(fit,scale="r2")
  dev.new()
  plot(fit,scale="adjr2")
  dev.new()
  plot(fit,scale="Cp")
  dev.new()
  plot(fit,scale="bic")
}
test_country <- function(df.test, iso, model) {
  df.test <- df.test[df.test$iso_code == iso,]
  to.plot <- data.frame(y_real = df.test$new_cases, date = df.test$date)
  df.test$new_cases <- NULL
  df.test$date <- NULL
  to.plot$y_pred <- predict(model, newdata=df.test)
  print(mean((to.plot$y_real-to.plot$y_pred)^2))
  plot <- ggplot(to.plot, aes(x=date)) +
    geom_line(aes(y=y_real, group=1)) +
    geom_line(aes(y=y_pred, group=2), color='red')
  print(plot)
}



#LOAD DATA##################
data_train = open_dataset(c('ITA', 'GBR', 'IND', 'JPN', 'ISR', 'LUX', 'AUS', 'AUT', 'NZL', 'ARG', 'BEL', 'CAN', 'ZAF', 'PRT','ISL','CHE'))
data_test = open_dataset( c('USA','IRN','KOR','URY'),TRUE,TRUE)
data_validation = open_dataset(c('RUS','TUR','DNK'))
#CONVERT DATA TO LOG RESPONSE#############################################
data_train.log <- data_train
data_train.log$new_cases <- log(data_train$new_cases + 0.1)
data_train.log$actual_cases = log(data_train$actual_cases + 0.1)
data_test.log <- data_test
data_test.log$new_cases <- log(data_test$new_cases + 0.1)
data_test.log$actual_cases = log(data_test$actual_cases + 0.1)
data_validation.log <- data_validation
data_validation.log$new_cases <- log(data_validation$new_cases + 0.1)
data_validation.log$actual_cases = log(data_validation$actual_cases + 0.1)


#CHECK DEPENDENCIES###################################
dev.new()
pairs(~new_cases+stringency_index+population_density+median_age+population,data_train)
dev.new()
pairs(~new_cases+aged_65_older+aged_70_older+gdp_per_capita,data_train)
dev.new()
pairs(~new_cases+cvd_death_rate+diabetes_prevalence+female_smokers+male_smokers,data_train)
dev.new()
pairs(~new_cases+new_tests+total_tests+actual_cases,data_train)

#FITTING LINEAR MODEL#########################
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
#FITTING LINEAR MODEL LOG OF RESPONSE######################################
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

#REFIT LINEAR MODEL WITH ONLY SIGNFICANT PREDICTORS###################
adj.fit = lm(new_cases~.-female_smokers-diabetes_prevalence-aged_70_older-cvd_death_rate,data_train)
summary(adj.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(adj.fit)
vif(adj.fit)

adj.fit = lm(new_cases~.-female_smokers-diabetes_prevalence-aged_70_older-cvd_death_rate-median_age-aged_65_older,data_train)
summary(adj.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(adj.fit)

#REFIT LINEAR MODEL (LOG RESPONSE) WITH ONLY SIGNFICANT PREDICTORS#########
adjlog.fit = lm(new_cases~.-female_smokers-actual_cases,data_train.log)
summary(adjlog.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(adjlog.fit)
vif(adjlog.fit)

#REFIT LINEAR MODEL (LOG RESPONSE) WITH ONLY SIGNFICANT PREDICTORS#########
adjlog.fit = lm(new_cases~.-female_smokers-actual_cases-population-population_density-median_age-aged_65_older-aged_70_older-cvd_death_rate,data_train.log)
summary(adjlog.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(adjlog.fit)
#BEST SUBSET SELECTION#######
best.fit=regsubsets(new_cases~.,data_train,nvmax = 15)
reg.summary=summary(best.fit)
#PLOT TO CHOOSE THE BEST NUMBER OF PREDICTORS
plot_best_predictors(reg.summary,best.fit)
#CHECK MIN ERROR ON VALIDATION SET
test.mat=model.matrix(new_cases~.,data=data_validation)
val.errors=rep(NA,14)
for(i in 1:14){
  coefi = coef(best.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_validation$new_cases-pred)^2)
}
val.errors; which.min(val.errors) 
coef(best.fit,which.min(val.errors)) # This is based on training data

best.fit = lm(new_cases~stringency_index+population+median_age+gdp_per_capita+total_cases+new_tests+total_tests,data_train)


# BACKWARD SELECTION ((((((((((((((((((DA QUI))))))))))))))))))
back.fit=regsubsets(data_train$new_cases~.,data_train,method = "backward",nvmax=15)
back.summary=summary(back.fit)
plot_best_predictors(back.summary,back.fit)

test.mat=model.matrix(new_cases~.,data=data_validation)
val.errors=rep(NA,14)
for(i in 1:14){
  coefi = coef(back.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_validation$new_cases-pred)^2)
}

val.errors; which.min(val.errors) 
min_error.back = coef(back.fit,which.min(val.errors)) 

# FORWARD SELECTION
for.fit=regsubsets(data_train$new_cases~.,data_train,method = "forward",nvmax=15)
for.summary=summary(for.fit)
plot_best_predictors(for.summary,for.fit)

test.mat=model.matrix(new_cases~.,data=data_validation)
val.errors=rep(NA,14)
for(i in 1:14){
  coefi = coef(for.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_validation$new_cases-pred)^2)
}
# The best model is the one that contains which.min(val.errors) (ten in the book) variables.
val.errors; which.min(val.errors) 
min_error.for = coef(for.fit,which.min(val.errors)) 

######################################################################################################################################################################
#BESTSUBSET SELECTION FOR LOG MODEL
best.log.fit=regsubsets(new_cases~.,data_train.log,nvmax = 15)
reg.summary.log=summary(best.log.fit)
#PLOT TO CHOOSE THE BEST NUMBER OF PREDICTORS
plot_best_predictors(reg.summary.log,best.log.fit)

#CHECK MIN ERROR ON VALIDATION SET
test.mat=model.matrix(new_cases~.,data=data_validation.log)
val.errors=rep(NA,14)
for(i in 1:14){
  coefi = coef(best.log.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_validation$new_cases-pred)^2)
}
# The best model is the one that contains which.min(val.errors) (ten in the book) variables.
val.errors; which.min(val.errors) 
min_error.log.best = coef(best.log.fit,which.min(val.errors)) # This is based on training data

# BACKWARD SELECTION FOR LOG MODEL
back.log.fit=regsubsets(data_train$new_cases~.,data_train,method = "backward",nvmax=15)
back.summary.log=summary(back.log.fit)
plot_best_predictors(back.summary.log,back.log.fit)

test.mat=model.matrix(new_cases~.,data=data_validation)
val.errors=rep(NA,14)
for(i in 1:14){
  coefi = coef(back.log.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_validation$new_cases-pred)^2)
}
# The best model is the one that contains which.min(val.errors) (ten in the book) variables.
val.errors; which.min(val.errors) 
min_error.back.log = coef(back.log.fit,which.min(val.errors)) 

# FORWARD SELECTION LOG MODEL
for.log.fit=regsubsets(data_train$new_cases~.,data_train,method = "forward",nvmax=15)
for.summary.log=summary(for.fit)
plot_best_predictors(for.summary.log,for.log.fit)

test.mat=model.matrix(new_cases~.,data=data_validation)
val.errors=rep(NA,14)
for(i in 1:14){
  coefi = coef(for.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_validation$new_cases-pred)^2)
}
# The best model is the one that contains which.min(val.errors) (ten in the book) variables.
val.errors; which.min(val.errors) 
min_error.for.log = coef(for.fit,which.min(val.errors)) 



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

