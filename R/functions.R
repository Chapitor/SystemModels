rm(list = ls())
library(usethis)
library(devtools)
library(roxygen2)


#' Load dependencies
pacman::p_load(tidyverse, optimx, caret, lubridate)


#' Initiate an example data frame
example_data <- data.frame("training_load" = c(rnorm(100, mean = 1000, sd = 150),
                                            rep(0,50)),
                        "rest" = c(rnorm(100, mean= 2, sd=1),
                                   rep(1,50)),
                        "perf" = c(seq(from = 10, to = 30, length.out = 100),
                                   rep(0,50)))
# usethis::use_data(example_data)

# adaptation -------------------------------------------------------------
#' To calculate adaptation from each training loads input
#'@param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#'@param k1 denotes the gain term of adaptations induced by training
#'@param tau1 denotes the time constant for the adaptation to decline
#'@param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#'
#'@return a vector of numerical values
#'
#' @examples
#'adaptation_fn(data = example_data, k1 = 0.5, tau1 = 40, vars = list("input" = example_data$training_load, "time" = example_data$rest))
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
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

# Fatigue -----------------------------------------------------------------

#' To make the fatigue impulse response to vary over time and function to accumulation of training sessions
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param k3 denotes gain term of fatigue induced by training
#' @param tau3 denotes the time constant for the fatigue remanence
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#'
#' @return a vector of numerical values
#'
#' @examples
#'k2i_fn(data = example_data, k3=5, tau3 = 5, vars = list("input" = example_data$training_load, "time" = example_data$rest))
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
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


#' To calculate the variable fatigue component for each training input
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param k3 denotes gain term of fatigue induced by training
#' @param tau2 denotes the time constant for the fatigue remanence
#' @param tau3 denotes the time constant for the fatigue to decline
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#'
#' @return a vector of numerical values
#'
#' @examples
#'fatigue_fn(data = example_data, k3=5, tau2 = 20, tau3 = 5, vars = list("input" = example_data$training_load, "time" = example_data$rest))
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
fatigue_fn <- function(data, k3, tau2, tau3, vars){
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

# Model -------------------------------------------------------------------

#' Extract the first performance
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param target is a character that indicates the performances column name.
#'
#' @return a numeric value
#'
#' @examples
#'init_perf(data = example_data, target = "perf")
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
init_perf <- function(data, target){
  return(data %>% dplyr::filter(target != 0) %>%
           data.table::first() %>%
           dplyr::select(target) %>%
           as.numeric())
}


#' Extract known performances
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param target is a character vector that indicates the performances column name.
#'
#' @return a vector of numerical values
#'
#' @examples
#'real_perf(data = example_data, target = "perf")
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
real_perf <- function(data, target){
  res <- NULL
  res <- data[,target]
  res[is.na(res)] <- 0
  return(res)
}

#' Function to model the performance following Busso Variable dose-response model (Busso, 2003)
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param P0, @param k1, @param k2, @param tau1, @param tau2, @param tau3 denote the basic level of performance, gains terms and time constants used in the model
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#' @param target is a character vector that indicates the performances column name.
#'
#' @return a vector of numerical values
#'
#' @examples
#'perf_model(data = example_data, PO = 10, k1 = 0.1, k3 = 0.01, tau1 = 40, tau2 = 20, tau3 = 5, vars = list("input" = example_data$training_load, "time" = example_data$rest),
#' target = "perf")
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
perf_model <- function(data, P0, k1, k3, tau1, tau2, tau3, vars, target){

  apt <- adaptation_fn(data, k1, tau1, vars)
  fat <- fatigue_fn(data, k3, tau3, tau2, vars)
  res <- vector(length = length(fat))
  P0 <- P0
  obs <- real_perf(data, target)

  for(i in 1:length(fat)){
    ifelse(obs[i] != 0,
           res[i] <- P0 + apt[i] - fat[i],
           res[i] <- 0)
  }
  return(res)
}



#' Function used for parameter estimate, called for optimization
#' @param data is a data frame that contains training loads and performances
#' @param theta is a vector of parameters to optimize
#' @param target is a character vector that indicates the performances column name.
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#'
#' @return a numerical value
#'
#' @examples
#' RSS(data = example_data, theta= c(10, 0.1, 0.01, 40, 20, 5), target = "perf", vars = list("input" = example_data$training_load, "time" = example_data$rest))
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
RSS <- function(data, theta, target, vars){
  y <- real_perf(data, target)
  y_hat <- perf_model(data, P0=theta[1], k1=theta[2], k3=theta[3], tau1=theta[4], tau2=theta[5], tau3=theta[6], vars, target)
  diff <- rep(0, length=length(y))

  for(i in 1:length(y)){
    if(y[i]!=0){
      diff[i] <- y[i]-y_hat[i]
    }
  }
  rss <- sum((diff)^2)
  return(rss)
}


# parameters tests --------------------------------------------------------

data = example_data %>%
  mutate(rest_days = rep(1,150))
# y = "PERF"
target="perf"
# var = c("Charge", "datetime", "PERF", "rest_days")
# vars <- list("input" = example_data$Charge, "time" = example_data$rest_days)
vars <- list("input" = example_data$training_load, "time" = example_data$rest)
k1 = 0.1
k3 = 0.01
tau1 = 40
tau2 = 20
tau3 = 5
P0 = init_perf(data = data, target = target)
theta <- c(P0, k1,k3,tau1,tau2,tau3)
optim.method = "nlm"


# DEV ---------------------------------------------------------------------

#' The main function that optimises parameters used in the Variable dose-response model (Busso, 2003), predicts values and calculates model performance.
#' @param data is a data frame that contains training loads and performances
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#' @param target is a character vector that indicates the performances column name.
#' @param specify default is NULL. If TRUE, a list of numeric vector of initial values for P0, k1, k3, tau1, tau2, tau3 parameters, lower and upper bounds has to be specified.
#' PO is the initial level of performance and can be extracted through the function init_perf
#' @param validation default is NULL. A list that contains the method used for validation (see details), initial window numerical value, horizon numerical values and
#' a logical for the fixed window.
#' @param optm.method is a character that indicates which method has to be used for parameter optimisation (see details)
#'
#' @details validation can be used to learn model and evaluate within a cross-validation procedure. Methods available are c("TS_CV") for time-series cross-validation.
#' For "TS-CV" method, specify numeric values for "initialWindow", "horizon" and logical term for the fixed window within a list.
#' Each of the optimization algorithm with constraints and used by the optimx package can be specified in @param optim.method.
#'
#'

#' @param CV indicates whether the model is evaluated within a time series cross-validation. TRUE / FALSE is needed. Default is FALSE.
#' @param date_obj is a datetime class object from data frame to be ordered. Default is NULL
#' @param x is a character which defines the name of input variable
#' @param y is the character vector which specifies the name of the performance column in data frame
#' @param z is a character which defines the time (days) between each input
#' @param init_param is a numeric vector that indicates initial values of parameters (k1, k3, tau1, tau2 and tau3) to be optimized
#' @param specify default is NULL. A list of numeric vector of initial values for P0, k1, k3, tau1, tau2, tau3 parameters, lower and upper bounds can be specified. PO is the initial level of performance and can be exctrated through the function init_perf
sysmod <- function(data, vars, target, specify = NULL, validation = NULL, optim.method){

  df <- data

  if (is.null(specify) == FALSE) {

    P0_init <- init_perf(data, target)
    theta <- specify[["theta"]]
    res_optim <- optimx::optimx(par = theta, fn = RSS, data = df, method = optim.method,
                                lower = specify[["lower"]],
                                upper = specify[["upper"]])

  } else {

    P0_init <- init_perf(data, target)
    k1_init <- 0.5
    k3_init <- 0.1
    tau1_init <- 40
    tau2_init <- 20
    tau3_init <- 5

    theta <- c(P0_init, k1_init, k3_init, tau1_init, tau2_init, tau3_init)


    res_optim <- optimx::optimx(par = theta, fn = RSS,data = data, target=target, vars=vars, method = optim.method,
                                lower = c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1),
                                upper = c(P0_init, 1, 1, 80, 40, 10))
  }

  P0 <- res_optim$p1
  k1 <- res_optim$p2
  k3 <- res_optim$p3
  tau1 <- res_optim$p4
  tau2 <- res_optim$p5
  tau3 <- res_optim$p6

  df$perf <- real_perf(data, y)
  df$adaptation <- adaptation_fn(data = df, k1 = k1, tau1 = tau1, vars = vars)
  df$k2i <- k2i_fn(data = df, k3 = k3, tau3 = tau3, vars = vars)
  df$fatigue <- fatigue_fn(data = df, k3 = k3, tau3 = tau3, tau2 = tau2, vars = vars)
  df$predicted <- perf_model(data=df, P0=P0, k1=k1, k3=k3, tau1=tau1, tau2=tau2, tau3=tau3, vars = vars, y = y)

  theta <- c(P0, k1, k3, tau1, tau2, tau3)

  list <- list("df"=df, "theta"=theta)

  return(list)
}

# sysmod(data = example_data, y = "PERF", var = c("Charge", "datetime", "PERF"), date_obj = "datetime")



