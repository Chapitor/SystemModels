rm(list = ls())

#' Load dependencies
pacman::p_load(tidyverse, optimx, caret, lubridate)


#' Initiate an example data frame
example_data <- read.csv(file = "./data/example.csv") %>%
  mutate(datetime = as_date(datetime))

# parameters tests --------------------------------------------------------

data = example_data %>%
  mutate(rest_days = rep(1,32))
y = "PERF"
var = c("Charge", "datetime", "PERF", "rest_days")
vars <- list("input" = example_data$Charge, "time" = example_data$rest_days)
TL = example_data$Charge
days = example_data$rest_days
k1 = 0.1
k3 = 0.01
tau1 = 40
tau2 = 20
tau3 = 5
P0 = init_perf(data = data, y = y)
theta <- c(P0, k1,k3,tau1,tau2,tau3)

# to1 = 40
# to2= 20
# to3= 5

# adaptation -------------------------------------------------------------
#' Calculate adaptation from each training loads input
#'@param k1 is the amplitude constant
#'@param tau1 is the time constant for the adaptation to decline
#'@param vars a list that contains input and time vectors. Input relates daily training loads and time is the time between two consecutive inputs.
adaptation_fn <- function(data, k1, tau1, vars){

  adapt_val <- vector(length=nrow(data))

  # Make the function return aberration if the time constant takes a negative or 0 value
  if (is.na(tau1)) adapt_val <- rep(-9000, nrow(data))
  else if(tau1 == 0) adapt_val <- rep(-9000, nrow(data))
  else {
    adapt_val[1] <- 0
    for(i in 2:nrow(data)){
      adapt_val[i] <- (k1*vars[["input"]][i-1] + adapt_val[i-1])*exp(-vars[["time"]][i]/(tau1))
    }
  }
  return(adapt_val)
}

# # DO NOT RUN
adaptation_fn(data = data, k1 = k1,tau1 = tau1,vars = list("input" = example_data$Charge, "time" = example_data$rest_days))

# Fatigue -----------------------------------------------------------------

#' Calculate the k2i variable required for the fatigue calculation
#' @param k3 is the amplitude constant
#' @param tau3 is the time constant for the fatigue to decline
k2i_fn <- function(data, k3, tau3, vars){

  k2i_val <- vector(length=nrow(data))

  if (is.na(tau3)) k2i_val <- rep(-9000, nrow(data))
  else if(tau3==0) k2i_val <- rep(-9000, nrow(data))
  else {
    k2i_val[1] <- 0
    for(i in 2:nrow(data)){
      k2i_val[i] <- (k3*vars[["input"]][i-1] + k2i_val[i-1])*exp(-vars[["time"]][i]/(tau3))
    }
  }
  return(k2i_val)
}


#' Calculate the fatigue component for each training imput
#' @param tau2 is the time constant for the fatigue remanence
fatigue_fn <- function(data, k3, tau3, tau2, vars){
  fat <- vector(length=nrow(data))
  if (is.na(tau3) | is.na(tau2)) apt <- rep(-9000, nrow(data))
  else if (tau3==0 | tau2==0) fatigue <- rep(-9000, nrow(data))
  else {
    fat[1] <- 0
    k2i <- k2i_fn(data, k3, tau3, vars)
    for(i in 2:nrow(data)){
      fat[i] <- (k2i[i-1]*vars[["input"]][i-1]+fat[i-1])*exp(-vars[["time"]][i]/(tau2))
    }
  }
  return(fat)
}

#' # DO NOT RUN
#' example_data <- fatigue_fn(data = example_data, k3 = 0.1, tau3 = 5, tau2 = 10, TL = example_data$training_load, days = example_data$rest)
#' plot(example_data$fatigue)


# Model -------------------------------------------------------------------

#' Retrieve the first performance further used in basic level of performance
#' @param y is a character vector that indicates the performances column name.
init_perf <- function(data, y){
  return(data %>% dplyr::filter(y != 0) %>%
           data.table::first() %>%
           dplyr::select(y) %>%
           as.numeric())
}


#' Real performances used in RSS
#' @param y is a character vector that indicates the performances column name.
real_perf <- function(data, y){
  res <- NULL
  res <- data[,y]
  res[is.na(res)] <- 0
  return(res)
}


#' @param data is the data frame that contains training loads andp erformances time ordered
#' @param P0, @param k1, @param k2, @param tau1, @param tau2, @param tau3 are basic level of performance, gains terms and time constants to be optimized through RSS
perf_model <- function(data, P0, k1, k3, tau1, tau2, tau3, vars, y){

  apt <- adaptation_fn(data, k1, tau1, vars)
  fat <- fatigue_fn(data, k3, tau3, tau2, vars)
  res <- vector(length = length(fat))
  P0 <- P0
  obs <- real_perf(data, y)

  for(i in 1:length(fat)){
    ifelse(obs[i] != 0,
           res[i] <- P0 + apt[i] - fat[i],
           res[i] <- 0)
  }
  return(res)
}



#' Function used for parameter estimate, called for optimization
#' @param theta is a vector of parameters to optimize
#' @param data is a data frame that contains training loads and performances
#' @param data_ini is the data frame thats contains the initial performance (First folder of time series CV)
RSS <- function(theta, data, y, vars){
  obs <- real_perf(data, y)
  pred <- perf_model(data, P0=theta[1], k1=theta[2], k3=theta[3], tau1=theta[4], tau2=theta[5], tau3=theta[6], vars, y)
  diff <- rep(0, length=length(y))

  for(i in 1:length(obs)){
    if(obs[i]!=0){
      diff[i] <- obs[i]-pred[i]
    }
  }
  rss <- sum((diff)^2)
  return(rss)
}



# DEV ---------------------------------------------------------------------


#' @param var is a character vector that indicates the names of input, performance and datetime object in the data frame
#' @param CV indicates whether the model is evaluated within a time series cross-validation. TRUE / FALSE is needed. Default is FALSE.
#' @param date_obj is a datetime class object from data frame to be ordered. Default is NULL
#' @param x is a character which defines the name of input variable
#' @param y is the character vector which specifies the name of the performance column in data frame
#' @param z is a character which defines the time (days) between each input
#' @param init_param is a numeric vector that indicates initial values of parameters (k1, k3, tau1, tau2 and tau3) to be optimized
#' @param specify default is NULL. A list of numeric vector of initial values for P0, k1, k3, tau1, tau2, tau3 parameters, lower and upper bounds can be specified. PO is the initial level of performance and can be exctrated through the function init_perf
sysmod <- function(data, vars, y, specify = NULL){

  df <- data

  if (is.null(specify) == FALSE) {

    P0_init <- init_perf(data, y)
    theta <- specify[["theta"]]
    res_optim <- optimx::optimx(par = theta, fn = RSS, data = df, method = "nlm",
                                lower = specify[["lower"]],
                                upper = specify[["upper"]])

  } else {

    P0_init <- init_perf(data, y)
    k1_init <- 1e-06
    k3_init <- 1e-08
    tau1_init <- 41
    tau2_init <- 26
    tau3_init <- 5

    theta <- c(P0_init, k1_init, k3_init, tau1_init, tau2_init, tau3_init)


    res_optim <- optimx::optimx(par = theta, fn = RSS,data = data, y=y, vars=vars, method = "nlm",
                                lower = c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1),
                                upper = c(P0_init, 1, 1, 80, 40, 10))

    res_optim <- optim(par = theta, fn = RSS, data = df)

  }



  P0 <- res_optim$p1
  k1 <- res_optim$p2
  k3 <- res_optim$p3
  tau1 <- res_optim$p4
  tau2 <- res_optim$p5
  tau3 <- res_optim$p6

  df$perf <- real_perf(data, y)
  df$adaptation <- adaptation_fn(df, k1, tau1)
  df$k2i <- k2i_fn(df,k3,tau3)
  df$fatigue <- fatigue_fn(df, k3, tau3, tau2)
  df$perf_mod <- perf_model(data=df, data_ini = df_ini, P0=P0, k1=k1, k3=k3, tau1=tau1, tau2=tau2, tau3=tau3)

  theta <- c(P0, k1, k3, tau1, tau2, tau3)

  list <- list("df"=df, "theta"=theta,)

  return(list)
}

# sysmod(data = example_data, y = "PERF", var = c("Charge", "datetime", "PERF"), date_obj = "datetime")



