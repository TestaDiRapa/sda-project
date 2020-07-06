#LOAD LIBRARIES AND DEFINE FUNCTIONS##############################
library(leaps)
library(car)
library(glmnet)
library(ggplot2)
path_a="/Users/Antonio/Desktop/sda-project"
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
  return(mean((to.plot$y_real-to.plot$y_pred)^2))
}.
predict.regsubsets = function(object,newdata,id,...){ 
  form=as.formula(object$call[[2]])
  mat=model.matrix(form, newdata)
  coefi=coef(object, id=id)
  xvars=names(coefi)
  mat[,xvars]%*%coefi
}
#LOAD DATA##################
data_train = open_dataset(c('ITA', 'GBR', 'IND', 'JPN', 'ISR', 'LUX', 'AUS', 'AUT', 'NZL', 'ARG', 'BEL', 'CAN', 'ZAF', 'PRT','ISL','CHE','RUS','TUR','DNK')) #-334 for validation
data_test = open_dataset( c('USA','IRN','KOR','URY'),TRUE,TRUE)
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
adj.fit = lm(new_cases~.-aged_65_older-diabetes_prevalence-male_smokers,data_train)
summary(adj.fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(adj.fit)
vif(adj.fit)

adj.fit = lm(new_cases~.-aged_65_older-diabetes_prevalence-male_smokers-median_age-aged_70_older,data_train)
summary(adj.fit)

best.fit.res = resid(adj.fit)
dev.new()
plot(data_train$new_cases, best.fit.res, ylab="Residuals", xlab="Samples", main="Residual Plot")
abline(0, 0)

#CONVERT DATA TO LOG RESPONSE#############################################
data_train$new_cases <- log(data_train$new_cases + 0.1)
data_train = na.omit(data_train.log)
data_train$actual_cases = log(data_train$actual_cases + 0.1)
data_test$new_cases <- log(data_test$new_cases + 0.1)
data_test$actual_cases = log(data_test$actual_cases + 0.1)


#FITTING LINEAR MODEL LOG OF RESPONSE######################################
fit = lm(new_cases~.,data_train)

#REMOVE ININFLUENTIAL POINTS
cooksd <- cooks.distance(fit)
# Plot the Cook's Distance using the traditional 4/n criterion
sample_size <- nrow(data_train)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4/sample_size, col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
data_train = data_train[which(!(rownames(data_train) %in% influential)),]


fit = lm(new_cases~.,data_train)
summary(fit)
# model diagnositic plots
dev.new()
par(mfrow=c(2,2))
plot(fit)
#REFIT LINEAR MODEL (LOG RESPONSE) WITH ONLY SIGNFICANT PREDICTORS
adj.fit = lm(new_cases~.-diabetes_prevalence,data_train)
summary(adj.fit)
vif(adj.fit)
adj.fit = lm(new_cases~.-diabetes_prevalence-population-population_density-median_age-aged_65_older-aged_70_older-cvd_death_rate,data_train)

#BEST SUBSET SELECTION VIA CROSS-VALIDATION#######

k=10
set.seed(1)
folds = sample(1:k,nrow(data_train),replace=TRUE)
cv.errors = matrix(NA,k,15, dimnames=list(NULL, paste(1:15)))
# write a for loop that performs cross-validation
for(j in 1:k){
  best.fit=regsubsets(new_cases~., data=data_train[folds!=j,], nvmax=15)
  print(j)
  for(i in 1:15){
    pred = predict(best.fit, data_train[folds==j,], id=i)
    print(i)
    print(length(pred))
    print(length(data_train[folds==j,]$new_cases))
    cv.errors[j,i] = mean((data_train[folds==j,]$new_cases-pred)^2)
  }
}
# This has given us a 10×19 matrix, of which the (i,j)th element corresponds to the test MSE for the i-th cross-validation fold for the best j-variable model.
mean.cv.errors=apply(cv.errors, 2, mean); mean.cv.errors# Column average
colMeans(cv.errors) # the same
par(mfrow=c(1,1))
dev.new()
plot(mean.cv.errors, type="b") #it selects an 10-variable model (11 on the textbook)
# We now perform best subset selection on the full data set to obtain the 10-variable model.
reg.best=regsubsets (Salary~., data=Hitters, nvmax=19)
# coef(reg.best, 11)
coef(reg.best, which.min(mean.cv.errors))

# BACKWARD SELECTION###############################
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
back.fit = lm(new_cases~stringency_index+population+population_density+median_age+aged_65_older+gdp_per_capita+new_cases+total_tests,data_train)

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
coef(for.fit,which.min(val.errors)) 

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
coef(best.log.fit,which.min(val.errors)) # This is based on training data
best.log.fit = lm(new_cases~stringency_index+population+population_density+median_age+aged_65_older+aged_70_older+gdp_per_capita+female_smokers+male_smokers+total_cases+new_tests+total_tests+actual_cases,data_train.log)
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
coef(back.log.fit,which.min(val.errors)) 
back.log.fit = lm(new_cases~stringency_index+population+population_density+median_age+aged_65_older+gdp_per_capita+new_tests+total_tests,data_train.log)
# BACKWARD SELECTION FOR LOG MODEL
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
coef(for.fit,which.min(val.errors)) 

#check best performance
data_validation = open_dataset(c('RUS','TUR','DNK'),TRUE,TRUE)
data_validation.log <- data_validation
data_validation.log$new_cases <- log(data_validation$new_cases + 0.1)
data_validation.log$actual_cases = log(data_validation$actual_cases + 0.1)
val.errors.best = 0
val.errors.back =0
for( iso in c('RUS','TUR','DNK')){
  val.errors.best = test_country(data_validation,iso,best.fit) + val.errors.best
  val.errors.back = test_country(data_validation,iso,back.fit) + val.errors.back
}
val.errors.best = val.errors.best/3
val.errors.back = val.errors.back/3
#the best subset selection outperforms
val.errors.log.best = 0
val.errors.log.back =0
for( iso in c('RUS','TUR','DNK')){
  val.errors.log.best = test_country(data_validation.log,iso,best.log.fit) + val.errors.log.best
  val.errors.log.back = test_country(data_validation.log,iso,back.log.fit) + val.errors.log.back
}
val.errors.log.best = val.errors.log.best/3
val.errors.log.back = val.errors.log.back/3
#RIDGE AND LASSO#######################
#The best subset selection outperforms
#variable selection with Ridge Regression and the Lasso
# use the glmnet package in order to perform ridge regression and the lasso. We do not use the y ??? x syntax here, but matrix and vector.
# perform ridge regression and the lasso in order to predict Salary on the Hitters data. Missing values has to be removed
# It has an alpha argument that determines what type of model is fit.
# If alpha = 0 a ridge regression model is fit.
# If alpha = 1 (default) then a lasso model is fit. 

x = model.matrix(new_cases~stringency_index+population+median_age+gdp_per_capita+total_cases+new_tests+total_tests, data_train)[,-1] # without 1's 
x_val = model.matrix(new_cases~stringency_index+population+median_age+gdp_per_capita+total_cases+new_tests+total_tests, data_validation)[,-1] # without 1's 
y = data_train$new_cases
y_val = data_validation$new_cases
x.log = model.matrix(best.log.fit, data_train.log)[,-1] # without 1's 
y.log = data_train$new_cases

grid=10^seq(15,-5,length=200)
ridge.mod=glmnet(x,y,alpha=0,lambda=grid)

log.grid=10^seq(3,-1,length=100)
log.ridge.mod=glmnet(x,y,alpha=0,lambda=log.grid)

ridge.pred=predict(ridge.mod,s=4,newx=x_val,exact=T,x=x,y=y) # Note the use of the predict() function again. This time we get predictions for a test set, by replacing type="coefficients" with the newx argument.
mean((ridge.pred-y_val)^2) # test MSE
x_merged = model.matrix(new_cases~stringency_index+population+median_age+gdp_per_capita+total_cases+new_tests+total_tests, data_merged)[,-1] # without 1's 
y_merged = data_merged$new_cases
cv.out=cv.glmnet(x_merged,y_merged,alpha=0)
dev.new()
plot(cv.out)
bestlam=cv.out$lambda.min; bestlam;log(bestlam) # the best lambda (212 on the text)


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

