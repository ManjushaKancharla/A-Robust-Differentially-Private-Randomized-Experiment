######### A Robust, Differentially Private Randomized Experiment ####
######## for Evaluating Online Educational Programs With Sensitive Student Data ####

#### Authors: Manjusha Kancharla and Hyunseung Kang

#### This .R file contains the code needed to reproduce results in Section 4 of our paper

#### Required R libraries

library(dplyr); library(plyr); library(forcats); library(utils); library(boot)


#### Read in cleaned data

set.seed(213) # For reproducibility

# setwd("~/Box/Research/Classroom Experiment/Data Analysis/Full Data May 16, 2021")
full_d <- read.csv("RPRCT_Data.csv", 
                   header = TRUE, na.strings=c(""," ","NA"))

d <- full_d[full_d$Incomplete.in.Outcomes == 'No',] # 72 usable observations

## Treatments : Trt1: Control (Slides only), Trt2: Treatment (Slides + Face)

d$Treatment <- revalue(d$Treatment, c("Trt1"="Control Group", "Trt2"="Treatment Group"))


## Create Numeric Columns of Ys (Binary Outcomes of Interest)

# Y1 <- "I found it hard to pay attention to the video". (Attention)
# Y2 <- "I was unable to recall concepts covered in the video while attempting the followup quiz". (Retention)
# Y3 <- “I don't feel that I learnt a great deal by watching the video”. (Self-judgement of learning)
# Y4 <- "I found the topic covered in the video to be hard".(Comprehension)

d$Y1_num <- NULL
d$Y1_num <- ifelse(d$Y1=="Yes",1,0)

d$Y2_num <- NULL
d$Y2_num <- ifelse(d$Y2=="Yes",1,0)

d$Y3_num <- NULL
d$Y3_num <- ifelse(d$Y3=="Yes",1,0)

d$Y4_num <- NULL
d$Y4_num <- ifelse(d$Y4=="Yes",1,0)

d$Trt <- NULL
d$Trt <- ifelse(d$Treatment=="Treatment Group",1,0)

## Randomization Probabilities used in Experiment

r0 <- r0_p <- 0
r1_p <- 0.104
r1 <- 0.1667

d <- d[d$Incomplete.in.Outcomes == 'No',]
delta <- mean(d$Trt)
n <- nrow(d)

## Clean Column names

colnames(d)[2] <- "Gender"
colnames(d)[4] <- "Race"
colnames(d)[5] <- "High School"
colnames(d)[6] <- "Statistics in HS"
colnames(d)[7] <- "College Year"
colnames(d)[8] <- "Major"
colnames(d)[10] <- "Prior Knowledge"
colnames(d)[11] <- "Number of MSC Courses"
colnames(d)[12] <- "Math Grade"
colnames(d)[13] <- "Stats Grade"
colnames(d)[14] <- "Comp. Sci. Grade"
colnames(d)[15] <- "Enjoy Math/Stat"
colnames(d)[16] <- "Hours Per Week"
colnames(d)[17] <- "GPA"
colnames(d)[18] <- "English FL"
colnames(d)[19] <- "Online Courses"
colnames(d)[20] <- "Online Courses - Face"
colnames(d)[21] <- "Hours of Video Lecs"
colnames(d)[22] <- "Preference"

## Collapse factor levels into fewer categories

d$StatGrade <- fct_collapse(d$`Stats Grade`,
                            "A-A+" = c("A", "A+"),
                            other = c("A-", "AB", "B", "B-", "B+", "D", "Other")
)

d$MathGrade <- fct_collapse(d$`Math Grade`,
                            "A-A+" = c("A", "A+"),
                            other = c("A-", "AB", "B", "B-", "C", "D", "Other")
)

d$`Online Courses - Face` <- fct_collapse(d$`Online Courses - Face`,
                                          "0-2" = c("None", "1", "2"),
                                          "3+" = c("3 or more")
)

d$`Hours of Video Lecs` <- fct_collapse(d$`Hours of Video Lecs`,
                                        "0-4" = c("0 hours", "< 2 hours", "2-4 hours"),
                                        "> 4" = c("> 4 hours")
)

d$`Enjoy Math/Stat` <- fct_collapse(d$`Enjoy Math/Stat`,
                                    "A lot" = c("A lot", "A great deal"),
                                    "Little to None" = c("A little", "None at all")
)

d$`High School` <- fct_collapse(d$`High School`,
                                "Other" = c("Private", "International")
)

d$`College Year` <- fct_collapse(d$`College Year`,
                                 "Senior" = c("Senior", "Other")
)

d[sapply(d, is.character)] <- lapply(d[sapply(d, is.character)], 
                                     as.factor)

###################### Y1 (Attention) Analysis ##############################

## Lambda hat estimation

se_func <- function(data,indices){
  d <- data[indices,]
  ll <- function(theta) { #theta = c(pi, gamma, lambda_0, lambda_1)
    pi <- theta[1]
    gamma <- theta[2]
    lambda_0 <- theta[3]
    lambda_1 <- theta[4]
    
    p <- lambda_1 + 0.5*(r1 + r1_p)*(gamma) +
      0.5*(2-r0-r0_p)*(pi)
    q <- lambda_0 + 0.5 * (r0+r0_p)* pi + 0.5 * (2-r1-r1_p)*gamma
    l <- sum(d$Y1_num) * log(p) + 
      sum(1-d$Y1_num) * log(q)
    #l <- -1*l
    return(l)
  }
  
  # Set pi = 0 
  
  
  grid <- expand.grid(pi = seq(0,1, 0.02), gamma= seq(0,1, 0.02), lambda_0 = seq(0,1, 0.02),
                      lambda_1 = seq(0,1, 0.02))
  grid$sum1 <- grid$pi + grid$gamma + grid$lambda_0 + grid$lambda_1
  grid <- grid[grid$sum1 ==1,]
  liklihoods <- apply(grid[,1:4], 1, ll)
  max <- which.max(liklihoods)
  lambda0_hat1 <- grid[max,3] 
  lambda1_hat1 <- grid[max,4]
  
  lambda_hat1 <- lambda0_hat1 + lambda1_hat1
  return(lambda_hat1)
}

results1 <- boot(d, se_func, R=5000, parallel="multicore", ncpus = 8) #se = 0.1296489
lambda_hat1 <- results1$t0
lambda1_var <- (sd(results1$t))^2 

## Estimate Tau_H,Diff
ybar_t <- with(d, mean(Y1_num*Trt)/mean(Trt))
ybar_c <- with(d, mean(Y1_num*(1-Trt))/mean(1-Trt))
denom <- (1-lambda_hat1)*(1/2)*(2-(r1+r1_p))
tau_hat1 <- (ybar_t-ybar_c)/denom

## Postulated model for A = 1, S = 1
fit1_trt1 <- glm(Y1_num ~ RD +`Enjoy Math/Stat` + MathGrade + `High School` + 
                   `Statistics in HS` +  `Hours of Video Lecs`+ Preference + `English FL`, 
                 data = d[(d$Trt == 1),], family = 'binomial')
fit1_ctl1 <- glm(Y1_num ~  RD +`Enjoy Math/Stat` + MathGrade + `High School` + 
                   `Statistics in HS` +  `Hours of Video Lecs`+ Preference + `English FL`, 
                 data = d[(d$Trt == 0),], family = 'binomial')

d$mu1_trt <- predict(fit1_trt1,newdata = d, type = 'response') 
d$mu1_ctl <- predict(fit1_ctl1, newdata = d, type = 'response') 

treated.dr1 = with(d, ((Y1_num - mu1_trt) * Trt)/delta + mu1_trt)
control.dr1 = with(d, ((Y1_num - mu1_ctl) *(1-Trt))/(1 - delta) + mu1_ctl)

## Estimate Tau_H,Cov

tau_hat_dr1 <-(mean(treated.dr1, na.rm = TRUE)-mean(control.dr1, na.rm = TRUE))/denom



## SE of Tau_hats

### SE of Tau_H,Diff

delta <- mean(d$Trt)
mult1 <- ((0.5)*(2-r1-r1_p))^{-2}*(1-lambda_hat1)^{-2}
sum1 <- (delta)^{-1} * with(d[d$Trt ==1,], var(Y1_num)) +
  (1-delta)^{-1} * with(d[d$Trt ==0,], var(Y1_num))
mult2 <- tau_hat1^2
sum2 <- 2+  lambda1_var*(1-lambda_hat1)^{-2}
var_tau1 <- mult1*sum1 + mult2 * sum2
var_tau1 <- var_tau1/n
se_tau1 <- sqrt(var_tau1)

### SE of Tau_H,Cov

sum1 <- with(d, mean((mu1_trt-mu1_ctl)^2,na.rm = TRUE)) + with(d[d$Trt ==1,], mean((Y1_num - mu1_trt)^2,na.rm = TRUE)*delta^{-1}) +
  with(d[d$Trt ==0,], mean((Y1_num - mu1_ctl)^2,na.rm = TRUE)*(1-delta)^{-1})
sum2 <- 1 + (lambda1_var)*(1-lambda_hat1)^{-2}
mult2 <- tau_hat_dr1^2

var_tau_dr1 <- mult1*sum1 + mult2 * sum2
var_tau_dr1 <- var_tau_dr1/n

se_tau_dr1 <- sqrt(var_tau_dr1)

######################### Y2 (Retention) ###############################


# Lambda hat estimation

se_func <- function(data,indices){
  d <- data[indices,]
  ll <- function(theta) { #theta = c(pi, gamma, lambda_0, lambda_1)
    pi <- theta[1]
    gamma <- theta[2]
    lambda_0 <- theta[3]
    lambda_1 <- theta[4]
    
    p <- lambda_1 + 0.5*(r1 + r1_p)*(gamma) +
      0.5*(2-r0-r0_p)*(pi)
    q <- lambda_0 + 0.5 * (r0+r0_p)* pi + 0.5 * (2-r1-r1_p)*gamma
    l <- sum(d$Y2_num) * log(p) + 
      sum(1-d$Y2_num) * log(q)
    #l <- -1*l
    return(l)
  }
  
  # Set pi = 0 
  
  
  grid <- expand.grid(pi = seq(0,1, 0.02), gamma= seq(0,1, 0.02), lambda_0 = seq(0,1, 0.02),
                      lambda_1 = seq(0,1, 0.02))
  grid$sum1 <- grid$pi + grid$gamma + grid$lambda_0 + grid$lambda_1
  grid <- grid[grid$sum1 ==1,]
  liklihoods <- apply(grid[,1:4], 1, ll)
  max <- which.max(liklihoods)
  lambda0_hat2 <- grid[max,3] 
  lambda1_hat2 <- grid[max,4]
  
  lambda_hat2 <- lambda0_hat2 + lambda1_hat2
  return(lambda_hat2)
}

results2 <- boot(d, se_func, R=5000, parallel="multicore", ncpus = 8) #se = 0.1804607

lambda_hat2 <- results2$t0
lambda2_var <- (sd(results2$t))^2 

## Estimate Tau_H,Diff
ybar_t <- with(d, mean(Y2_num*Trt)/mean(Trt))
ybar_c <- with(d, mean(Y2_num*(1-Trt))/mean(1-Trt))
denom <- (1-lambda_hat2)*(1/2)*(2-(r1+r1_p))
tau_hat2 <- (ybar_t-ybar_c)/denom

## Postulated model for A = 1, S = 1
fit2_trt1 <- glm(Y2_num ~ RD + MathGrade + StatGrade + `Enjoy Math/Stat` + 
                   GPA + Gender + Preference  + `English FL` + `Hours of Video Lecs`, 
                 data = d[(d$Trt == 1),], family = 'binomial')

fit2_ctl1 <- glm(Y2_num ~  RD + MathGrade + StatGrade + `Enjoy Math/Stat` +  
                   GPA + Gender + Preference+ `English FL` + `Hours of Video Lecs`, 
                 data = d[(d$Trt == 0),], family = 'binomial')

d$mu2_trt <- predict(fit2_trt1,newdata = d, type = 'response') 
d$mu2_ctl <- predict(fit2_ctl1, newdata = d, type = 'response') 

treated.dr2 = with(d, ((Y2_num - mu2_trt) * Trt)/delta + mu2_trt)
control.dr2 = with(d, ((Y2_num - mu2_ctl) *(1-Trt))/(1 - delta) + mu2_ctl)

## Estimate Tau_H,Cov

tau_hat_dr2 <-(mean(treated.dr2, na.rm = TRUE)-mean(control.dr2, na.rm = TRUE))/denom



## SE of Tau_hats
### SE of Tau_H,Diff

mult1 <- ((0.5)*(2-r1-r1_p))^{-2}*(1-lambda_hat2)^{-2}
sum1 <- (delta)^{-1} * with(d[d$Trt ==1,], var(Y2_num)) +
  (1-delta)^{-1} * with(d[d$Trt ==0,], var(Y2_num))
mult2 <- tau_hat2^2
sum2 <- 2+  lambda2_var*(1-lambda_hat2)^{-2}
var_tau2 <- mult1*sum1 + mult2 * sum2
var_tau2 <- var_tau2/n
se_tau2 <- sqrt(var_tau2)

### SE of Tau_H,Cov

sum1 <- with(d, mean((mu2_trt-mu2_ctl)^2,na.rm = TRUE)) + with(d[d$Trt ==1,], mean((Y2_num - mu2_trt)^2,na.rm = TRUE)*delta^{-1}) +
  with(d[d$Trt ==0,], mean((Y2_num - mu2_ctl)^2,na.rm = TRUE)*(1-delta)^{-1})
sum2 <- 1 + (lambda2_var)*(1-lambda_hat2)^{-2}
mult2 <- tau_hat_dr2^2

var_tau_dr2 <- mult1*sum1 + mult2 * sum2
var_tau_dr2 <- var_tau_dr2/n
se_tau_dr2 <- sqrt(var_tau_dr2)


###################### Y3 (JOL) ###########################

# Lambda hat estimation

se_func <- function(data,indices){
  d <- data[indices,]
  ll <- function(theta) { #theta = c(pi, gamma, lambda_0, lambda_1)
    pi <- theta[1]
    gamma <- theta[2]
    lambda_0 <- theta[3]
    lambda_1 <- theta[4]
    
    p <- lambda_1 + 0.5*(r1 + r1_p)*(gamma) +
      0.5*(2-r0-r0_p)*(pi)
    q <- lambda_0 + 0.5 * (r0+r0_p)* pi + 0.5 * (2-r1-r1_p)*gamma
    l <- sum(d$Y3_num) * log(p) + 
      sum(1-d$Y3_num) * log(q)
    #l <- -1*l
    return(l)
  }
  
  grid <- expand.grid(pi = seq(0,1, 0.02), gamma= seq(0,1, 0.02), lambda_0 = seq(0,1, 0.02),
                      lambda_1 = seq(0,1, 0.02))
  grid$sum1 <- grid$pi + grid$gamma + grid$lambda_0 + grid$lambda_1
  grid <- grid[grid$sum1 ==1,]
  liklihoods <- apply(grid[,1:4], 1, ll)
  max <- which.max(liklihoods)
  lambda0_hat3 <- grid[max,3] 
  lambda1_hat3 <- grid[max,4]
  
  lambda_hat3 <- lambda0_hat3 + lambda1_hat3
  return(lambda_hat3)
}

results3 <- boot(d, se_func, R=5000, parallel="multicore", ncpus = 8) #se =  0.1486595
lambda_hat3 <- results3$t0
lambda3_var <- (sd(results3$t))^2 

## Estimate Tau_H,Diff

ybar_t <- with(d, mean(Y3_num*Trt)/mean(Trt))
ybar_c <- with(d, mean(Y3_num*(1-Trt))/mean(1-Trt))

denom <- (1-lambda_hat3)*(1/2)*(2-(r1+r1_p))
tau_hat3 <- (ybar_t-ybar_c)/denom



## Postulated model for A = 1, S = 1
fit1_trt1 <- glm(Y3_num ~ RD + MathGrade + StatGrade + `Enjoy Math/Stat` +
                   GPA + Gender + Preference + Major, 
                 data = d[(d$Trt == 1),], family = 'binomial')

fit1_ctl1 <- glm(Y3_num ~  RD + MathGrade + StatGrade + `Enjoy Math/Stat` +
                   GPA + Gender + Preference + Major, 
                 data = d[(d$Trt == 0),], family = 'binomial')


d$mu1_trt <- predict(fit1_trt1,newdata = d, type = 'response') 
d$mu1_ctl <- predict(fit1_ctl1, newdata = d, type = 'response')

treated.dr3 = with(d, ((Y3_num - mu1_trt) * Trt)/delta + mu1_trt)
control.dr3 = with(d, ((Y3_num - mu1_ctl) *(1-Trt))/(1 - delta) + mu1_ctl)

## Estimate Tau_H,Cov

tau_hat_dr3 <-(mean(treated.dr3, na.rm = TRUE)-mean(control.dr3, na.rm = TRUE))/denom



### SE of Tau_H,Diff

mult1 <- ((0.5)*(2-r1-r1_p))^{-2}*(1-lambda_hat3)^{-2}
sum1 <- (delta)^{-1} * with(d[d$Trt ==1,], var(Y3_num)) +
  (1-delta)^{-1} * with(d[d$Trt ==0,], var(Y3_num))
mult2 <- tau_hat3^2
sum2 <- 2+  lambda3_var*(1-lambda_hat3)^{-2}
var_tau3 <- mult1*sum1 + mult2 * sum2
var_tau3 <- var_tau3/n
se_tau3 <- sqrt(var_tau3)

### SE of Tau_H,Cov

sum1 <- with(d, mean((mu1_trt-mu1_ctl)^2,na.rm = TRUE)) + with(d[d$Trt ==1,], mean((Y3_num - mu1_trt)^2,na.rm = TRUE)*delta^{-1}) +
  with(d[d$Trt ==0,], mean((Y3_num - mu1_ctl)^2,na.rm = TRUE)*(1-delta)^{-1})
sum2 <- 1 + (lambda3_var)*(1-lambda_hat3)^{-2}
mult2 <- tau_hat_dr3^2

var_tau_dr3 <- mult1*sum1 + mult2 * sum2
var_tau_dr3 <- var_tau_dr3/n

se_tau_dr3 <- sqrt(var_tau_dr3)



############################ Y4 (Comprehension) #########################

# Lambda hat estimation

se_func <- function(data,indices){
  d <- data[indices,]
  ll <- function(theta) { #theta = c(pi, gamma, lambda_0, lambda_1)
    pi <- theta[1]
    gamma <- theta[2]
    lambda_0 <- theta[3]
    lambda_1 <- theta[4]
    
    p <- lambda_1 + 0.5*(r1 + r1_p)*(gamma) +
      0.5*(2-r0-r0_p)*(pi)
    q <- lambda_0 + 0.5 * (r0+r0_p)* pi + 0.5 * (2-r1-r1_p)*gamma
    l <- sum(d$Y4_num) * log(p) + 
      sum(1-d$Y4_num) * log(q)
    #l <- -1*l
    return(l)
  }
  
  # Set pi = 0 
  
  
  grid <- expand.grid(pi = seq(0,1, 0.02), gamma= seq(0,1, 0.02), lambda_0 = seq(0,1, 0.02),
                      lambda_1 = seq(0,1, 0.02))
  grid$sum1 <- grid$pi + grid$gamma + grid$lambda_0 + grid$lambda_1
  grid <- grid[grid$sum1 ==1,]
  liklihoods <- apply(grid[,1:4], 1, ll)
  max <- which.max(liklihoods)
  lambda0_hat4 <- grid[max,3] 
  lambda1_hat4 <- grid[max,4]
  
  lambda_hat4 <- lambda0_hat4 + lambda1_hat4
  return(lambda_hat4)
}

results4 <- boot(d, se_func, R=5000, parallel="multicore", ncpus = 8) #se = 0.1998648
lambda_hat4 <- results4$t0
lambda4_var <- (sd(results4$t))^2 

## Estimate Tau_H,Diff

ybar_t <- with(d, mean(Y4_num*Trt)/mean(Trt))
ybar_c <- with(d, mean(Y4_num*(1-Trt))/mean(1-Trt))

denom <- (1-lambda_hat4)*(1/2)*(2-(r1+r1_p))
tau_hat4 <- (ybar_t-ybar_c)/denom



## Postulated model for A = 1, S = 1
fit1_trt1 <- glm(Y4_num ~ RD + MathGrade  + `Enjoy Math/Stat` + `College Year`, 
                 data = d[(d$Trt == 1),], family = 'binomial')

fit1_ctl1 <- glm(Y4_num ~  RD + MathGrade  + `Enjoy Math/Stat` + `College Year`, 
                 data = d[(d$Trt == 0),], family = 'binomial')


d$mu1_trt <- predict(fit1_trt1,newdata = d, type = 'response') 
d$mu1_ctl <- predict(fit1_ctl1, newdata = d, type = 'response')

treated.dr4 = with(d, ((Y4_num - mu1_trt) * Trt)/delta + mu1_trt)
control.dr4 = with(d, ((Y4_num - mu1_ctl) *(1-Trt))/(1 - delta) + mu1_ctl)

## Estimate Tau_H,Cov
tau_hat_dr4 <-(mean(treated.dr4, na.rm = TRUE)-mean(control.dr4, na.rm = TRUE))/denom



### SE of Tau_H,Diff

mult1 <- ((0.5)*(2-r1-r1_p))^{-2}*(1-lambda_hat4)^{-2}
sum1 <- (delta)^{-1} * with(d[d$Trt ==1,], var(Y4_num)) +
  (1-delta)^{-1} * with(d[d$Trt ==0,], var(Y4_num))
mult2 <- tau_hat4^2
sum2 <- 2+  lambda4_var*(1-lambda_hat4)^{-2}
var_tau4 <- mult1*sum1 + mult2 * sum2
var_tau4 <- var_tau4/n
se_tau4 <- sqrt(var_tau4)

### SE of Tau_H,Cov

sum1 <- with(d, mean((mu1_trt-mu1_ctl)^2,na.rm = TRUE)) + with(d[d$Trt ==1,], mean((Y4_num - mu1_trt)^2,na.rm = TRUE)*delta^{-1}) +
  with(d[d$Trt ==0,], mean((Y4_num - mu1_ctl)^2,na.rm = TRUE)*(1-delta)^{-1})
sum2 <- 1 + (lambda4_var)*(1-lambda_hat4)^{-2}
mult2 <- tau_hat_dr4^2

var_tau_dr4 <- mult1*sum1 + mult2 * sum2
var_tau_dr4 <- var_tau_dr4/n

se_tau_dr4 <- sqrt(var_tau_dr4)

