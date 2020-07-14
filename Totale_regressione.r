library(car)
library(ggplot2)
library(glmnet)
library(leaps)
library(MASS)
library(tidyverse)
library(hrbrthemes)
setwd("D:/Documenti HDD/Uni/Magistrale/Statistical Data Analysis/Esercizi/Progetto")
set.seed(42)
########################################################################################################
### FUNCTIONS
open_dataset_past_cases <- function(iso_codes, dates=FALSE, iso=FALSE) {
  df.data <- read.csv('covid-data.csv')
  if(dates) {
    columns <- colnames(df.data)[c(1,4,9,15,16,20,21,22,23,24,25,26,28,29,30,31)]
  } else {
    columns <- colnames(df.data)[c(1,9,15,16,20,21,22,23,24,25,26,28,29,30,31)]
  }
  df.data <- df.data[which(df.data$iso_code %in% iso_codes), columns]
  df.data <- na.omit(df.data)
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
  if(!iso) {
    df.data$iso_code <- NULL
  }
  return(df.data)
}

plot_best_predictors = function(summary,fit){
  #dev.new()
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
}

predict.regsubsets = function(object,newdata,id,...){ 
  form=as.formula(object$call[[2]])
  mat=model.matrix(form, newdata)
  coefi=coef(object, id=id)
  xvars=names(coefi)
  mat[,xvars]%*%coefi
}

ridge_regression <- function(x, y) {
  grid=10^seq(10,-2,length=100)
  ridge.mod=glmnet(x,y,alpha=0,lambda=grid)
}
########################################################################################################
### DATA IMPORT
iso_codes_train <- c('ITA', 'GBR', 'IND', 'JPN', 'ISR', 'LUX', 'AUS', 'AUT', 'NZL', 'ARG', 'BEL', 'ZAF', 'PRT', 'CHE', 'ISL', 'RUS','TUR','DNK')
iso_codes_test <- c('USA','IRN','KOR','URY')
covid.data <- open_dataset_past_cases(iso_codes_train)
covid.data.test <- open_dataset_past_cases(iso_codes_test, dates = TRUE, iso = TRUE)
########################################################################################################
### CHECKING NONLINEAR RELATIONSHIP BETWEEN VARIABLES AND RESPONSE
dev.new()
pairs(~total_cases+stringency_index+population_density+median_age+population,covid.data)
dev.new()
pairs(~total_cases+aged_65_older+aged_70_older+gdp_per_capita,covid.data)
dev.new()
pairs(~total_cases+cvd_death_rate+diabetes_prevalence+female_smokers,covid.data)
dev.new()
pairs(~total_cases+new_tests+total_tests+actual_cases,covid.data)
########################################################################################################
### FITTING LINEAR MODEL WITH ALL PREDICTORS
fit <- lm(total_cases~., data=covid.data)
summary(fit)

# Diagnostic plots
dev.new()
par(mfrow=c(2,2))
plot(fit)

# Removing outliers
dev.new()
cooksd <- cooks.distance(fit)
sample_size <- nrow(covid.data)
plot(cooksd, pch="*", cex=2, main="Influential Points by Cooks distance")
abline(h = 4/sample_size, col="red")
cook_th <- 4/sample_size
influential <- names(cooksd)[(cooksd > cook_th)]
covid.data <- covid.data[which(!(rownames(covid.data) %in% influential)),]

# Model re-fit
fit <- lm(total_cases~., data=covid.data)
vif(fit)

#Adjustinhg model using p-value and vif
adj.fit <- lm(total_cases~.-aged_70_older-aged_65_older-population, data=covid.data)
vif(adj.fit)
summary(adj.fit)

adj.fit <- lm(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, data=covid.data)
summary(adj.fit)

# Diagnostic plots (again)
dev.new()
par(mfrow=c(2,2))
plot(fit)

# Residual plot
fit.res = resid(adj.fit)
dev.new()
plot(covid.data$total_cases, fit.res, ylab="Residuals", xlab="Samples", main="Residual Plot")
abline(0, 0)

########################################################################################################
### BEST SUBSET SELECTION USING K-FOLD
k <- 10

predictors <- 14
folds <- sample(1:k,nrow(covid.data),replace=TRUE)
cv.errors <- matrix(NA,k,predictors, dimnames=list(NULL, paste(1:predictors)))
# write a for loop that performs cross-validation
for(j in 1:k){
  best.fit <- regsubsets(total_cases~., data=covid.data[folds!=j,], nvmax=predictors)
  for(i in 1:predictors){
    pred <- predict(best.fit, covid.data[folds==j,], id=i)
    cv.errors[j,i] <- mean((covid.data[folds==j,]$total_cases-pred)^2)
  }
}
# MSE for each model using K-fold
mean.cv.errors <- apply(cv.errors, 2, mean)
par(mfrow=c(1,1))
dev.new()
plot(mean.cv.errors, type="b") 

best.subsets <- regsubsets(total_cases~., data = covid.data, nvmax = 14)
summary(best.subsets)
best.fit <- lm(total_cases~.-gdp_per_capita-diabetes_prevalence-new_tests-actual_cases, data=covid.data)
plot_best_predictors(summary(best.subsets), best.subsets)

### BEST SUBSET STARTING FROM ADJUSTED MODEL

predictors <- 9
folds <- sample(1:k,nrow(covid.data),replace=TRUE)
cv.errors <- matrix(NA,k,predictors, dimnames=list(NULL, paste(1:predictors)))
# write a for loop that performs cross-validation
for(j in 1:k){
  best.fit.adj <- regsubsets(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, data=covid.data[folds!=j,], nvmax=predictors)
  for(i in 1:predictors){
    pred <- predict(best.fit.adj, covid.data[folds==j,], id=i)
    cv.errors[j,i] <- mean((covid.data[folds==j,]$total_cases-pred)^2)
  }
}
# MSE for each model using K-fold
mean.cv.errors <- apply(cv.errors, 2, mean)
par(mfrow=c(1,1))
dev.new()
plot(mean.cv.errors, type="b", xlab="#Predictors", ylab="MSE") 

best.subsets <- regsubsets(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, data = covid.data, nvmax = 14)
summary(best.subsets)
best.fit.adj <- lm(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate-new_tests, data=covid.data)
plot_best_predictors(summary(best.subsets), best.subsets)

### BACKWARDS SELECTION
train <- sample(c(TRUE,FALSE), nrow(covid.data), rep=TRUE)
test <- (!train)
back.fit <- regsubsets(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, covid.data[train,], method = "backward", nvmax=14)
back.summary <- summary(back.fit)
plot_best_predictors(back.summary, back.fit)

test.mat <- model.matrix(total_cases~.,data=covid.data[test,])
val.errors=rep(NA,predictors)
for(i in 1:predictors){
  coefi <- coef(back.fit, id=i)
  pred <- test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((covid.data[test,]$total_cases-pred)^2)
}

val.errors
which.min(val.errors) 
back.fit = lm(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate-new_tests, data = covid.data)

### FORWARD SELECTION
train <- sample(c(TRUE,FALSE), nrow(covid.data), rep=TRUE)
test <- (!train)
for.fit <- regsubsets(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, covid.data[train,], method = "forward", nvmax=14)
for.summary <- summary(for.fit)
plot_best_predictors(for.summary, for.fit)

test.mat <- model.matrix(total_cases~.,data=covid.data[test,])
val.errors=rep(NA,predictors)
for(i in 1:predictors){
  coefi <- coef(for.fit, id=i)
  pred <- test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((covid.data[test,]$total_cases-pred)^2)
}

val.errors
which.min(val.errors) 
for.fit = lm(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate-new_tests-new_tests, data = covid.data)
########################################################################################################
### RIDGE REGRESSION
train <- sample(c(TRUE,FALSE), nrow(covid.data),rep=TRUE)
test <- (!train)

x <- model.matrix(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, covid.data)[,-1]
y <- covid.data$total_cases
grid <- 10^seq(10,-2,length=100)

cv.out <- cv.glmnet(x, y, alpha=0, lambda=grid)
dev.new()
plot(cv.out)

bestlam <- cv.out$lambda.min
cv.out$lambda.1se
ridge.mod <- glmnet(x[train,], y[train], alpha=0, lambda=bestlam, thresh=1e-12)
ridge.pred <- predict(ridge.mod, newx=x[test,])
mean((ridge.pred-y[test])^2)

out <- glmnet(x, y, alpha=0, lambda=grid)
predict(out, type="coefficients", s=bestlam)[1:9,]
dev.new()
plot(out,label = T, xvar = "lambda")
abline(v = log(bestlam), col="red", lwd=2, lty=2)

ridge.mod <- glmnet(x, y, alpha=0, lambda=bestlam, thresh=1e-12)

### LASSO REGRESSION
cv.out <- cv.glmnet(x, y, alpha=1, lamda=grid)
dev.new()
plot(cv.out)

bestlam <- cv.out$lambda.1se

out <- glmnet(x, y, alpha=1, lambda=grid)
predict(out, type="coefficients", s=bestlam)[1:9,]
dev.new()
plot(out,label = T, xvar = "lambda")
abline(v = log(bestlam), col="red", lwd=2, lty=2)

lasso.mod <- glmnet(x, y, alpha=1, lambda=bestlam, thresh=1e-12)
########################################################################################################

x.test <- model.matrix(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, covid.data.test[,-c(1, 2)])[,-1]
x.train <- model.matrix(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate, covid.data)[,-1]
y.test <- covid.data.test$total_cases
y.train <- covid.data$total_cases

type <- rep(c("Test", "Train"), 6)
labels <- c("all_predictors", "all_predictors", "vif_selected", "vif_selected", "best_subset_all", "best_subset_all", "best_subset_9", "best_subset_9", "ridge", "ridge", "lasso", "lasso")
test.mse <- c(1:12)
df.plot <- covid.data.test
#df.plot$date <- substr(df.plot$date, 6, 10)

y.fit <- predict(fit, covid.data.test[,-c(1, 2)])
df.plot$all_predictors <- y.fit
test.mse[1] <- mean((y.test - y.fit)^2)
test.mse[2] <- mean((y.train - predict(fit, covid.data))^2)

y.adj.fit = predict(adj.fit, covid.data.test[,-c(1, 2)])
df.plot$vif_selected <- y.adj.fit
test.mse[3] <- mean((y.test - y.adj.fit)^2)
test.mse[4] <- mean((y.train - predict(adj.fit, covid.data))^2)

y.best.fit = predict(best.fit, covid.data.test[,-c(1, 2)])
df.plot$best_subset_all <- y.best.fit
test.mse[5] <- mean((y.test - y.best.fit)^2)
test.mse[6] <- mean((y.train - predict(best.fit, covid.data))^2)

y.best.fit.adj = predict(best.fit.adj, covid.data.test[,-c(1, 2)])
df.plot$best_subset_9 <- y.best.fit.adj
test.mse[7] <- mean((y.test - y.best.fit.adj)^2)
test.mse[8] <- mean((y.train - predict(best.fit.adj, covid.data))^2)

y.back.fit = predict(back.fit, covid.data.test[,-c(1, 2)])
mean((y.test - y.back.fit)^2)

y.for.fit = predict(for.fit, covid.data.test[,-c(1, 2)])
mean((y.test - y.for.fit)^2)

y.ridge <- predict(ridge.mod, newx=x.test)
df.plot$ridge <- y.ridge
test.mse[9] <- mean((y.test - y.ridge)^2)
test.mse[10] <- mean((y.train - predict(ridge.mod, newx=x.train))^2)

y.lasso <- predict(lasso.mod, newx=x.test)
df.plot$lasso <- y.lasso
test.mse[11] <- mean((y.test - y.lasso)^2)
test.mse[12] <- mean((y.train - predict(lasso.mod, newx=x.train))^2)

result.df <- data.frame(models = labels, mse = test.mse, type = type)
p <- ggplot(result.df, aes(x = models, y = mse, fill = type)) +
      geom_bar(stat="identity", position=position_dodge()) +
      #geom_text(aes(label=test.mse), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5) +
      scale_fill_brewer(palette="Paired") +
      theme_minimal()
p


test_multiple_country <- function(df.test, iso) {
  line.width <- 1.05
  to.plot <- df.test[df.test$iso_code == iso,]
  to.plot$date  <- as.Date(to.plot$date)
  plot <- ggplot(to.plot, aes(x=date)) +
    geom_line(aes(y=total_cases, colour='black', group = 1), size=line.width, ) +
    geom_line(aes(y=all_predictors, colour='red', group = 2), size=line.width) +
    geom_line(aes(y=best_subset_9, colour='cyan', group = 3), size=line.width) +
    geom_line(aes(y=best_subset_all, colour='blue', group = 4), size=line.width) +
    geom_line(aes(y=vif_selected, colour='green', group = 5), size=line.width) +
    geom_line(aes(y=ridge, colour='magenta', group = 6), size=line.width) +
    geom_line(aes(y=lasso, colour='yellow', group = 7), size=line.width) + 
    scale_color_discrete(name = "Legend", labels = c("real_cases", "all_predictors", "vif_selected", "best_subset_all", "best_subset_9","ridge", "lasso")) +
    scale_x_date("Days", breaks = "10 days") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(plot)
}


test_multiple_country(df.plot, 'KOR')
test_multiple_country(df.plot, 'USA')
test_multiple_country(df.plot, 'URY')
test_multiple_country(df.plot, 'IRN')


single_state_prediction <- function(country, split) {
  s.train <- open_dataset_past_cases(c(country))[c(1:split),]
  s.test <- open_dataset_past_cases(c(country), iso = TRUE, dates = TRUE)
  plot.s <- data.frame(date=as.Date(s.test$date), total_cases=s.test$total_cases)
  adj.fit.s <- lm(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate-population_density-gdp_per_capita-diabetes_prevalence-female_smokers-male_smokers-new_tests, data=s.train)
  summary(adj.fit.s)
  pred <- rep(0, length(s.test$total_cases))
  pred[1:split] <- predict(adj.fit.s, s.test[c(1:split),-c(1, 2)])
  # s.test$total_tests[(split+1):length(pred)] <- s.test$total_tests[split] + s.test$new_tests[split]*c(1:(length(pred)-split))
  # s.test$stringency_index[(split+1):length(pred)] <- 0
  for(i in (split+1):length(pred)) {
    s.test[i,]$actual_cases <- pred[i-1]
    pred[i] <- predict(adj.fit.s, s.test[i,-c(1, 2)])
  }
  plot.s$prediction <- pred
  line.width <- 1.05
  plot <- ggplot(plot.s, aes(x=date)) +
    geom_line(aes(y=total_cases, colour='black', group = 1), size=line.width, ) +
    geom_line(aes(y=prediction, colour='green', group = 2), size=line.width) +
    scale_color_discrete(name = "Legend", labels = c("real_cases", "prediction")) +
    geom_vline(xintercept = plot.s$date[split], linetype="dashed", color = "blue", size=line.width) +
    scale_x_date("Days", breaks = "5 days") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(plot)
  return(mean((pred[(split+1):length(pred)] - plot.s$total_cases[(split+1):length(pred)])^2))
}

s.mse <- rep(0, 4)
s.mse[1] <- single_state_prediction('JPN', 50)
s.mse[2] <- single_state_prediction('USA', 50)
s.mse[3] <- single_state_prediction('URY', 50)
s.mse[4] <- single_state_prediction('IRN', 25)
single.state.mse <- mean(s.mse)


single_state_prediction.2 <- function(country, split, model) {
  s.train <- open_dataset_past_cases(c(country))[c(1:split),]
  s.test <- open_dataset_past_cases(c(country), iso = TRUE, dates = TRUE)
  s.test.t <- open_dataset_past_cases(c(country), iso = TRUE, dates = TRUE)
  plot.s <- data.frame(date=as.Date(s.test$date), total_cases=s.test$total_cases)
  adj.fit.s <- lm(total_cases~.-aged_70_older-aged_65_older-population-median_age-cvd_death_rate-population_density-gdp_per_capita-diabetes_prevalence-female_smokers-male_smokers-new_tests, data=s.train)
  summary(adj.fit.s)
  pred <- rep(0, length(s.test$total_cases))
  pred_t <- rep(0, length(s.test$total_cases))
  pred[1:split] <- predict(adj.fit.s, s.test[c(1:split),-c(1, 2)])
  pred_t[1:split] <- predict(model, s.test[c(1:split),-c(1, 2)])
  # s.test$total_tests[(split+1):length(pred)] <- s.test$total_tests[split] + s.test$new_tests[split]*c(1:(length(pred)-split))
  # s.test$stringency_index[(split+1):length(pred)] <- 0
  for(i in (split+1):length(pred)) {
    s.test[i,]$actual_cases <- pred[i-1]
    pred[i] <- predict(adj.fit.s, s.test[i,-c(1, 2)])
    s.test.t[i,]$actual_cases <- pred_t[i-1]
    pred_t[i] <- predict(model, s.test.t[i, -c(1, 2)])
  }
  plot.s$prediction <- pred
  plot.s$prediction_all <- pred_t
  line.width <- 1.05
  plot <- ggplot(plot.s, aes(x=date)) +
    geom_line(aes(y=total_cases, colour='black', group = 1), size=line.width, ) +
    geom_line(aes(y=prediction, colour='green', group = 2), size=line.width) +
    geom_line(aes(y=prediction_all, colour='magenta', group = 3), size=line.width) +
    scale_color_discrete(name = "Legend", labels = c("real_cases", "prediction", "prediction_all")) +
    geom_vline(xintercept = plot.s$date[split], linetype="dashed", color = "blue", size=line.width) +
    scale_x_date("Days", breaks = "5 days") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(plot)
  return(mean((pred_t[(split+1):length(pred)] - plot.s$total_cases[(split+1):length(pred)])^2))
}
t.mse <- rep(0, 4)
t.mse[1] <- single_state_prediction.2('KOR', 50, adj.fit)
t.mse[2] <- single_state_prediction.2('USA', 50, adj.fit)
t.mse[3] <- single_state_prediction.2('URY', 50, adj.fit)
t.mse[4] <- single_state_prediction.2('IRN', 25, adj.fit)
all.states.mse <- mean(t.mse)

mse.plot <- data.frame(model=c("single state", "multi state"), mse=c(single.state.mse, all.states.mse))
p <- ggplot(mse.plot, aes(x = model, y = mse)) +
  geom_bar(stat="identity", position=position_dodge(), fill="lightblue") +
  scale_fill_brewer(palette="Paired") +
  theme_minimal()
p

single_state_prediction.3 <- function(country, split, model) {
  s.train <- open_dataset_past_cases(c(country))[c(1:split),]
  s.test <- open_dataset_past_cases(c(country), iso = TRUE, dates = TRUE)
  s.test.t <- open_dataset_past_cases(c(country), iso = TRUE, dates = TRUE)
  s.test.i <- open_dataset_past_cases(c(country), iso = TRUE, dates = TRUE)
  plot.s <- data.frame(date=as.Date(s.test$date), total_cases=s.test$total_cases)
  pred <- rep(0, length(s.test$total_cases))
  pred_t <- rep(0, length(s.test$total_cases))
  pred_i <- rep(0, length(s.test$total_cases))
  pred[1:split] <- predict(model, s.test[c(1:split),-c(1, 2)])
  pred_t[1:split] <- predict(model, s.test[c(1:split),-c(1, 2)])
  pred_i[1:split] <- predict(model, s.test[c(1:split),-c(1, 2)])
  # s.test$total_tests[(split+1):length(pred)] <- s.test$total_tests[split] # + s.test$new_tests[split]*c(1:(length(pred)-split))
  # s.test.i$total_tests[(split+1):length(pred)] <- 2*s.test.t$total_tests[(split+1):length(pred)]
  s.test$stringency_index[(split+1):length(pred)] <- seq(s.test$stringency_index[split], 0, length=length(pred)-split)
  s.test.i$stringency_index[(split+1):length(pred)] <- seq(s.test$stringency_index[split], 100, length=length(pred)-split)
  for(i in (split+1):length(pred)) {
    s.test[i,]$actual_cases <- pred[i-1]
    pred[i] <- predict(model, s.test[i,-c(1, 2)])
    s.test.t[i,]$actual_cases <- pred_t[i-1]
    pred_t[i] <- predict(model, s.test.t[i, -c(1, 2)])
    s.test.i[i,]$actual_cases <- pred_i[i-1]
    pred_i[i] <- predict(model, s.test.i[i, -c(1, 2)])
  }
  plot.s$prediction <- pred
  plot.s$prediction_all <- pred_t
  plot.s$prediction_inc <- pred_i
  line.width <- 1.05
  plot <- ggplot(plot.s, aes(x=date)) +
    geom_line(aes(y=total_cases, colour='black', group = 1), size=line.width, ) +
    geom_line(aes(y=prediction, colour='green', group = 2), size=line.width) +
    geom_line(aes(y=prediction_all, colour='magenta', group = 3), size=line.width) +
    geom_line(aes(y=prediction_inc, colour='pink', group = 4), size=line.width) +
    scale_color_discrete(name = "Legend", labels = c("real_cases", "reduced stringency", "real stringency", "increased stringency")) +
    geom_vline(xintercept = plot.s$date[split], linetype="dashed", color = "blue", size=line.width) +
    scale_x_date("Days", breaks = "5 days") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(plot)
}

single_state_prediction.3('KOR', 50, adj.fit)
single_state_prediction.3('USA', 50, adj.fit)
single_state_prediction.3('URY', 50, adj.fit)
single_state_prediction.3('IRN', 25, adj.fit)
