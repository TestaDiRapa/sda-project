#### LOAD LIBRARIES ####
library(leaps)
library(car)
library(glmnet)
library(ggplot2)
#path_a="/Users/Antonio/Desktop/sda-project"
path_v="/Users/vincenzo/Documents/GitHub/sda-project"
setwd(path_v)
set.seed(42)

#### FUNCTION DEFINITION ####

#### LOAD DATA 
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
  actual_cases <- c() #inizialize the actual_cases vector which contains the daily increase for each country
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

##### PLOT FUNCTION
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

#### TEST FUNCTION, PLOT DAY BY DAY
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
}

#### BEST SUBSET SELECTION 
predict.regsubsets = function(object,newdata,id,...){ 
  form=as.formula(object$call[[2]])
  mat=model.matrix(form, newdata)
  coefi=coef(object, id=id)
  xvars=names(coefi)
  mat[,xvars]%*%coefi
}

my_predict = function(model,data_test,y_test){
  y_pred = predict(model,data_test)
  return(mean((y_test-y_pred)^2))
}

#### LOAD DATA ####
data_train = open_dataset(c('ITA', 'GBR', 'IND', 'JPN', 'ISR', 'LUX', 'AUS', 'AUT', 'NZL', 'ARG', 'BEL', 'CAN', 'ZAF', 'PRT','ISL','CHE','RUS','TUR','DNK')) #-334 for validation
data_test = open_dataset( c('USA','IRN','KOR','URY'),TRUE,TRUE) #'USA',

#### CHECK DEPENDENCIES
dev.new()
pairs(~new_cases+stringency_index+population_density+median_age+population,data_train)
dev.new()
pairs(~new_cases+aged_65_older+aged_70_older+gdp_per_capita,data_train)
dev.new()
pairs(~new_cases+cvd_death_rate+diabetes_prevalence+female_smokers+male_smokers,data_train)
dev.new()
pairs(~new_cases+new_tests+total_tests+actual_cases,data_train)

#### FITTING LINEAR MODEL ####

fit = lm(new_cases~.,data_train)
summary(fit)

#### REMOVE ININFLUENTIAL POINTS ACCORDING TO THE COOK'S DISTANCE

cooksd <- cooks.distance(fit)
# Plot the Cook's Distance using the traditional 4/n criterion
sample_size <- nrow(data_train)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4/sample_size, col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
data_train = data_train[which(!(rownames(data_train) %in% influential)),]

####FITTING LINEAR MODEL WITHOUT ININFLUENT POINTS
fit = lm(new_cases~.,data_train)
summary(fit)

####MODEL DIAGNOSTIC PLOTS
dev.new()
par(mfrow=c(2,2))
plot(fit)

#### CHECK PREDICTORS p-value
adj.fit = lm(new_cases~.-aged_65_older-diabetes_prevalence-male_smokers,data_train)
summary(adj.fit)

#### MODEL DIAGNOSTIC PLOTS
dev.new()
par(mfrow=c(2,2))
plot(adj.fit)

#### CHECK INTERACTION OF PREDICTORS
vif(adj.fit)

#### SELECTION 
adj.fit = lm(new_cases~.-aged_65_older-diabetes_prevalence-male_smokers-aged_70_older,data_train)
summary(adj.fit)

#### CHECK INTERACTION OF PREDICTORS
vif(adj.fit)

#### SELECTION
adj.fit = lm(new_cases~.-aged_65_older-diabetes_prevalence-male_smokers-aged_70_older-population_density,data_train)
summary(adj.fit)

#### CHECK INTERACTION OF PREDICTORS
vif(adj.fit)
adj.fit = lm(new_cases~.-aged_65_older-diabetes_prevalence-male_smokers-aged_70_older-population_density-female_smokers,data_train)

#### PLOT RESIDUAL
best.fit.res = resid(adj.fit)
dev.new()
plot(data_train$new_cases, best.fit.res, ylab="Residuals", xlab="Samples", main="Residual Plot")
abline(0, 0)

#### CONVERT DATA TO LOG RESPONSE ####
data_train$actual_cases = log(data_train$actual_cases + 0.00001)
data_test$actual_cases = log(data_test$actual_cases + 0.00001)
data_train = na.omit(data_train)
data_test = na.omit(data_test)

#### FITTING LINEAR MODEL LOG OF PREDICTOR ####
fit = lm(new_cases~.,data_train)

#REMOVE ININFLUENTIAL POINTS
cooksd <- cooks.distance(fit)

#### PLOT THE COOK'S DISTANCE USING THE TRADITIONAL 4/n CRITERION
sample_size <- nrow(data_train)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4/sample_size, col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
data_train = data_train[which(!(rownames(data_train) %in% influential)),]

#### FITTING LINEAR MODEL WITHOUT ININFLUENT POINTS 
fit = lm(new_cases~.,data_train)
summary(fit)

#### MODEL DIAGNOSTIC PLOTS
dev.new()
par(mfrow=c(2,2))
plot(fit)

#### CHECK INTERACTION OF PREDICTORS
vif(fit)

#### SELECTION 
adj.fit = lm(new_cases~.-aged_70_older,data_train)

#### CHECK INTERACTION OF PREDICTORS
vif(adj.fit)

#### SELECTION 
adj.fit = lm(new_cases~.-aged_70_older-median_age,data_train)

#### CHECK INTERACTION OF PREDICTORS
vif(adj.fit)

#### SELECTION 
adj.fit = lm(new_cases~.-aged_70_older-median_age-cvd_death_rate,data_train)

#### CHECK INTERACTION OF PREDICTORS
vif(adj.fit)

#### RESIDUAL PLOT
best.fit.res = resid(adj.fit)
dev.new()
plot(data_train$new_cases, best.fit.res, ylab="Residuals", xlab="Samples", main="Residual Plot")
abline(0,0)

#### FITTING BEST MODEL ACCORDING p-value AND VIF 
adj.fit = lm(new_cases~.-aged_70_older-median_age-cvd_death_rate-population_density-gdp_per_capita-female_smokers,data_train)
vif(adj.fit)

#### BEST SUBSET SELECTION VIA CROSS-VALIDATION ####

k=10
folds = sample(1:k,nrow(data_train),replace=TRUE)
cv.errors = matrix(NA,k,15, dimnames=list(NULL, paste(1:15)))
# write a for loop that performs cross-validation
for(j in 1:k){
  best.fit=regsubsets(new_cases~., data=data_train[folds!=j,], nvmax=15)
  for(i in 1:15){
    pred = predict(best.fit, data_train[folds==j,], id=i)
    cv.errors[j,i] = mean((data_train[folds==j,]$new_cases-pred)^2)
  }
}

mean.cv.errors = colMeans(cv.errors) # Column average
par(mfrow=c(1,1))
dev.new()
plot(mean.cv.errors, type="b")

# We now perform best subset selection on the full data set to obtain the 15-variable model.
reg.best=regsubsets (new_cases~., data=data_train, nvmax=15)
best.fit = lm(new_cases~.,data_train)

# BACKWARD SELECTION OF THE 9 PREDICTORS MODEL VIA VALIDATION APPROACH #####
train=sample(c(TRUE,FALSE), nrow(data_train),rep=TRUE)
test=(!train)
back.fit= regsubsets(new_cases~stringency_index+population+aged_65_older+diabetes_prevalence+male_smokers+total_cases+new_tests+total_tests+actual_cases,data_train[train,],method = "backward",nvmax=15)
back.summary= summary(back.fit)
plot_best_predictors(back.summary,back.fit)
test.mat = model.matrix(new_cases~.,data=data_train[test,])
val.errors=rep(NA,9)
for(i in 1:9){
  coefi = coef(back.fit,id=i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((data_train[test,]$new_cases-pred)^2)
}

which.min(val.errors) #find the model with minimum error
min_error.back = coef(back.fit,which.min(val.errors)) 

#### FITTING THE LINEAR MODEL WITH 7 PREDICTORS (similar error with respect to 9 predictors model)
back.fit = lm(new_cases~stringency_index+aged_65_older+male_smokers+new_tests+total_tests+actual_cases+total_cases,data_train)

##### RIDGE WITH 9 PREDICTORS MODEL #### 

train=sample(c(TRUE,FALSE), nrow(data_train),rep=TRUE)
test=(!train)

x = model.matrix(new_cases~stringency_index+population+aged_65_older+diabetes_prevalence+male_smokers+new_tests+total_tests+actual_cases+total_cases,data_train)[,-1] # without 1's 
x_val = model.matrix(new_cases~stringency_index+population+aged_65_older+diabetes_prevalence+male_smokers+new_tests+total_tests+actual_cases+total_cases,data_train)[,-1] # without 1's 
y = data_train$new_cases

grid=10^seq(10,-2,length=500)

#finding the best lambda through the cross-validation
cv.out=cv.glmnet(x,y,alpha=0)
dev.new()
plot(cv.out)
bestlamr=cv.out$lambda.min; bestlamr
cv.out$lambda.1se
ridge.mod=glmnet(x[train,],y[train],alpha=0,lambda=bestlamr,thresh=1e-12)
ridge.pred=predict(ridge.mod,newx=x[test,])
mean((ridge.pred-y[test])^2)

# REFITTING THE RIDGE REGRESSION MODEL WITH THE BEST LAMBDA 
out=glmnet(x,y,alpha=0)
predict(out,type="coefficients",s=bestlamr)
# As expected, none of the coefficients are zero
# ridge regression does not perform variable selection!
dev.new()
plot(out,label = T, xvar = "lambda")

ridge.mod = glmnet(x,y,alpha=0, lambda=bestlamr)

out=glmnet(x,y,alpha=1,lambda=grid)
ridge.coef=predict(out,type="coefficients",s=bestlamr)
ridge.coef


##### LASSO WITH 9 PREDICTORS MODEL #####

#finding the best lambda through the cross-validation
cv.out=cv.glmnet(x[train,],y[train],alpha=1)
dev.new()
plot(cv.out)
bestlaml=cv.out$lambda.min
print(cv.out$lambda.1se)
lasso.mod = glmnet(x[train,], y[train], alpha=1, lambda=bestlaml)
lasso.pred=predict(lasso.mod,s=bestlaml ,newx=x[test,])
mean((lasso.pred-y[test])^2) # slighly larger than ridge
out=glmnet(x,y,alpha=1,lambda=grid)
lasso.coef=predict(out,type="coefficients",s=bestlaml)


#### RIDGE WITH 7 PREDICTORS MODEL ####
x = model.matrix(new_cases~stringency_index+aged_65_older+male_smokers+new_tests+total_tests+actual_cases+total_cases,data_train)[,-1] # without 1's 
x_val = model.matrix(new_cases~stringency_index+aged_65_older+male_smokers+new_tests+total_tests+actual_cases+total_cases,data_train)[,-1] # without 1's 
y = data_train$new_cases
grid=10^seq(10,-2,length=500)

#finding the best lambda through the cross-validation
cv.out=cv.glmnet(x,y,alpha=0)
dev.new()
plot(cv.out)
bestlamr7=cv.out$lambda.min;
cv.out$lambda.1se
ridge7.mod=glmnet(x[train,],y[train],alpha=0,lambda=bestlamr7,thresh=1e-12)
ridge7.pred=predict(ridge7.mod,newx=x[test,])
mean((ridge7.pred-y[test])^2)

# REFIT RIDGE REGRESSION MODEL ON FULL DATA SET WITH THE BEST LAMBDA
out=glmnet(x,y,alpha=0)
predict(out,type="coefficients",s=bestlamr7)
# As expected, none of the coefficients are zero
# ridge regression does not perform variable selection!
dev.new()
plot(out,label = T, xvar = "lambda")

ridge7.mod = glmnet(x,y,alpha=0, lambda=bestlamr7)

out=glmnet(x,y,alpha=1,lambda=grid)
ridge7.coef=predict(out,type="coefficients",s=bestlamr7)

#### LASSO WITH 7 PREDICTORS MODEL ####
#finding the best lambda through the cross-validation
cv.out=cv.glmnet(x[train,],y[train],alpha=1)
dev.new()
plot(cv.out)
bestlaml7=cv.out$lambda.min
print(cv.out$lambda.1se)
lasso7.mod = glmnet(x[train,], y[train], alpha=1, lambda=bestlaml7)
lasso7.pred=predict(lasso7.mod,s=bestlaml7 ,newx=x[test,])
mean((lasso7.pred-y[test])^2) # slighly larger than ridge
out=glmnet(x,y,alpha=1,lambda=grid)
lasso7.coef=predict(out,type="coefficients",s=bestlaml7)

####### TEST PERFORMANCE #####
print(my_predict(adj.fit,data_test,data_test$new_cases))

print(my_predict(best.fit,data_test,data_test$new_cases))

print(my_predict(back.fit,data_test,data_test$new_cases))

data_test = data_test[,-c(1,2)]
x_test = model.matrix(new_cases~stringency_index+population+aged_65_older+diabetes_prevalence+male_smokers+new_tests+total_tests+actual_cases+total_cases,data_test)[,-1] # without 1's 
y_test = data_test$new_cases

ridge.pred=predict(ridge.mod,s=bestlamr,x_test)
mean((ridge.pred-y_test)^2)

lasso.pred=predict(lasso.mod,s=bestlaml,x_test)
mean((lasso.pred-y_test)^2)

x_test = model.matrix(new_cases~stringency_index+aged_65_older+male_smokers+new_tests+total_tests+actual_cases+total_cases,data_test)[,-1] # without 1's 
y_test = data_test$new_cases

ridge7.pred=predict(ridge7.mod,s=bestlamr,x_test)
mean((ridge7.pred-y_test)^2)

lasso7.pred=predict(lasso7.mod,s=bestlaml,x_test)
mean((lasso7.pred-y_test)^2)


