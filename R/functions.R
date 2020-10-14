rm(list = ls())
library(usethis)
library(devtools)
library(roxygen2)

# Initiate an example data frame
load("./data/example_data.rda")

#' @title  Adaptations response modeling
#' @description To model adaptation from each training session
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param k1 denotes the gain term for adaptations occurrences.
#' @param tau1 denotes the time constant of the adaptations exponential decrease.
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#'
#' @return a vector of numerical values
#'
#' @examples
#'adaptation_fn(data = example_data, k1 = 0.5, tau1 = 40, vars = list("input" = example_data$training_load, "time" = example_data$rest))
#'
#' @author Frank Imbach <frankimbach@gmail.com>
#' @export
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


#' @title Variable fatigue
#' @description To make the fatigue impulse response varying over time and function to accumulation of training sessions
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param k3 denotes the gain term for fatigue occurrences.
#' @param tau3 denotes the time constant of the fatigue remanence exponential decrease.
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#'
#' @return a vector of numerical values
#'
#' @examples
#'k2i_fn(data = example_data, k3=5, tau3 = 5, vars = list("input" = example_data$training_load, "time" = example_data$rest))
#'
#' @author Frank Imbach <frankimbach@gmail.com>
#' @export
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


#' @title Fatigue response modeling
#' @description To calculate the variable fatigue component for each training input
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions
#' @param k3 denotes the gain term for fatigue occurrences.
#' @param tau2 denotes the time constant of fatigue exponential decrease.
#' @param tau3 denotes the time constant of the fatigue remanence exponential decrease.
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


#' @title Initial performance
#' @description This function extracts the first performance of a time ordered data frame object.
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
  require(tidyverse)
  require(data.table)

  return(data %>% dplyr::filter(target != 0) %>%
           data.table::first() %>%
           dplyr::select(target) %>%
           as.numeric())
}


#' @title Extract observations
#' @description This function extracts the known performances from a data frame object.
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

#' @title Variable dose-response modeling
#' @description The function models the performance according to Busso Variable dose-response model (Busso, 2003) and based on previously defined parameters.
#' @param data is a data frame object that contains training loads, performances and time between two consecutive sessions.
#' @param P0 denotes the basic level of performance.
#' @param k1 denotes the gain term for adaptations occurrences.
#' @param k3 denotes the gain term for fatigue occurrences.
#' @param tau1 denotes the time constant of the adaptations exponential decrease.
#' @param tau2 denotes the time constant of fatigue exponential decrease.
#' @param tau3 denotes the time constant of the fatigue remanence exponential decrease.
#' @param vars denotes a list that contains input and time numeric vectors. Input refers to session training loads and time refers to the time between two consecutive inputs.
#' @param target is a character vector that indicates the performances column name.
#'
#' @return a vector of numerical values
#'
#' @examples
#'perf_model(data = example_data, P0 = 10, k1 = 0.1, k3 = 0.01, tau1 = 40, tau2 = 20, tau3 = 5, vars = list("input" = example_data$training_load, "time" = example_data$rest),
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



#' @title Residual sum of squares calculation
#' @description This function calculates the residual sum of squares from observations and predicted performances.
#' @param data is a data frame that contains training loads and performances
#' @param theta is a vector of parameters P0, k1, k3, tau1, tau2, tau3.
#' @param target is a character that indicates the performance column name.
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


#' @title Dose-response modeling
#' @description The functions train a dose-response model with or without cross-validation.
#' @param data is a data frame that contains training loads, performances, a date time values and the time between two consecutive sessions.
#' @param vars denotes a list in which "input" (e.g. training loads) and "time" between two consecutive inputs are displayed.
#' @param target is a character that indicates the performances column name.
#' @param date_ID is a character that indicates the date time object name.
#' @param specify default is NULL. If TRUE, a list of numeric vector of initial values for "P0", "k1", "k3", "tau1", "tau2", "tau3" parameters, "lower" and "upper" bounds has to be specified.
#' "PO" is the initial level of performance and can be extracted through the function [init_perf].
#' @param validation default is NULL. A list that contains the method used for validation (see details), initial window numerical value, horizon numerical values and
#' a logical for the fixed window.
#' @param optm.method is a character that indicates which method has to be used for parameter optimisation (see details).
#'
#' @details validation can be used to learn model and evaluate within a cross-validation procedure. Methods available are \code{"TS-CV"} for time-series cross-validation.
#' For \code{"TS-CV"} method, specify numeric values for \code{"initialWindow"}, \code{"horizon"} and logical term for the \code{"fixedWindow"} within a list.
#' Each of the optimization algorithm with constraints and used by the optimx package can be specified in optim.method.
#'
#' @return A list describing the model output and its performances.
#'
#' @note Model performances (RMSE, MAE and R squared) are calculated on test data for validation = c("simple", "TS-CV").
#'
#' @examples
#' sysmod(data = example_data, vars = list("input" = example_data$training_load, "time" = example_data$rest), target = "perf", specify = NULL,
#' validation = list("initialWindow" = 50, "horizon" = 15, "fixedWindow" = FALSE), optim.method = "nlm", date_ID = "datetime")
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@export
sysmod <-
  function(data,
           vars,
           target,
           date_ID,
           specify = NULL,
           validation = validation,
           optim.method) {

    require(tidyverse)
    require(optimx)
    require(caret)

    df <-
      data %>% dplyr::slice(-c(which(data[, target] == 0)))   # A reduced data frame that allow to model with no performance days. Used for time slices.

    # Initiate and optimize parameters
    if (is.null(specify) == FALSE) {
      P0_init <- init_perf(data, target)
      theta_init <- specify[["theta"]]


    } else {
      P0_init <- init_perf(data, target)
      k1_init <- 0.5
      k3_init <- 0.1
      tau1_init <- 40
      tau2_init <- 20
      tau3_init <- 5

      theta_init <-
        c(P0_init, k1_init, k3_init, tau1_init, tau2_init, tau3_init)
    }

    # Initiate empty objects for further saves
    theta_df <- data.frame()   # To save each optim
    dfs <- data.frame()
    rmse_vec <- c()
    MAE_vec <- c()
    Rsq_vec <- c()

    # Create timeslices
    if (is.null(validation) == FALSE) {
      time_slice <- caret::createTimeSlices(
        df[, target],
        initialWindow = validation[["initialWindow"]],
        horizon = validation[["horizon"]],
        fixedWindow = validation[["fixedWindow"]]
      )

      for (k in 1:length(time_slice$train)) {
        # Compute the model for each folder. The last folder will return the full model

        # split performance dataframe
        folder_train <- data[unlist(time_slice$train[k]),]
        folder_test <- data[unlist(time_slice$test[k]),]

        # retrieve date information of split
        datetime_train_min <- min(folder_train[, date_ID])
        datetime_train_max <- max(folder_train[, date_ID])

        datetime_test_min <- min(folder_test[, date_ID])
        datetime_test_max <- max(folder_test[, date_ID])

        # split datetime dataframe
        folder_train <-
          data %>% dplyr::filter(datetime <= datetime_train_max &
                            datetime >= datetime_train_min) %>%
          mutate("base" = "train")
        folder_test <-
          data %>% dplyr::filter(datetime <= datetime_test_max &
                            datetime >= datetime_test_min)  %>%
          mutate("base" = "test")
        folder_test <- rbind(folder_train, folder_test)


        # Model training
        if (is.null(specify) == FALSE) {
          res_optim <-
            optimx::optimx(
              par = theta_init,
              fn = RSS,
              data = folder_train,
              method = optim.method,
              lower = specify[["lower"]],
              upper = specify[["upper"]]
            )
        } else {
          res_optim <-
            optimx::optimx(
              par = theta_init,
              fn = RSS,
              data = folder_train,
              target = target,
              vars = vars,
              method = optim.method,
              lower = c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1),
              upper = c(P0_init, 1, 1, 80, 40, 10)
            )
        }

        P0 <- res_optim$p1
        k1 <- res_optim$p2
        k3 <- res_optim$p3
        tau1 <- res_optim$p4
        tau2 <- res_optim$p5
        tau3 <- res_optim$p6

        theta <-
          data.frame(
            P0 = P0,
            k1 = k1,
            k3 = k3,
            tau1 = tau1,
            tau2 = tau2,
            tau3 = tau3
          )

        folder_test$perf <-
          real_perf(data = folder_test, target = target)
        folder_test$adaptation <-
          adaptation_fn(
            data = folder_test,
            k1 = k1,
            tau1 = tau1,
            vars = vars
          )
        folder_test$k2i <-
          k2i_fn(
            data = folder_test,
            k3 = k3,
            tau3 = tau3,
            vars = vars
          )
        folder_test$fatigue <-
          fatigue_fn(
            data = folder_test,
            k3 = k3,
            tau3 = tau3,
            tau2 = tau2,
            vars = vars
          )
        folder_test$predicted <-
          perf_model(
            data = folder_test,
            P0 = P0,
            k1 = k1,
            k3 = k3,
            tau1 = tau1,
            tau2 = tau2,
            tau3 = tau3,
            vars = vars,
            target = target
          )
        folder_test$folder <- rep(k, nrow(folder_test))

        dfs <- rbind(dfs, folder_test)
        theta_df <- rbind(theta_df, theta)

        # Compute RMSE, MAE and Rsquared on test data
        rmse_vec <-
          c(rmse_vec,
            caret::RMSE(pred = folder_test[which(folder_test[, "base"] == "test"), "predicted"],
                        obs = folder_test[which(folder_test[, "base"] == "test"), "perf"]))
        MAE_vec <-
          c(MAE_vec,
            caret::MAE(pred = folder_test[which(folder_test[, "base"] == "test"), "predicted"],
                       obs = folder_test[which(folder_test[, "base"] == "test"), "perf"]))
        Rsq_vec <-
          c(Rsq_vec,
            caret::R2(pred = folder_test[which(folder_test[, "base"] == "test"), "predicted"],
                      obs = folder_test[which(folder_test[, "base"] == "test"), "perf"]))


      } # Close kfold
      return(
        list(
          "dfs" = dfs,
          "theta" = theta_df,
          "rmse_vec" = rmse_vec,
          "MAE_vec" = MAE_vec,
          "Rsq_vec" = Rsq_vec
        )
      )
    } else {
      # Close validation

      # Model training
      if (is.null(specify) == FALSE) {
        res_optim <-
          optimx::optimx(
            par = theta_init,
            fn = RSS,
            data = df,
            method = optim.method,
            lower = specify[["lower"]],
            upper = specify[["upper"]]
          )
      } else {
        res_optim <-
          optimx::optimx(
            par = theta_init,
            fn = RSS,
            data = df,
            target = target,
            vars = vars,
            method = optim.method,
            lower = c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1),
            upper = c(P0_init, 1, 1, 80, 40, 10)
          )
      }

      P0 <- res_optim$p1
      k1 <- res_optim$p2
      k3 <- res_optim$p3
      tau1 <- res_optim$p4
      tau2 <- res_optim$p5
      tau3 <- res_optim$p6

      theta <- c(
        P0 = P0,
        k1 = k1,
        k3 = k3,
        tau1 = tau1,
        tau2 = tau2,
        tau3 = tau3
      )

      df$perf <- real_perf(data = df, target = target)
      df$adaptation <-
        adaptation_fn(
          data = df,
          k1 = k1,
          tau1 = tau1,
          vars = vars
        )
      df$k2i <- k2i_fn(
        data = df,
        k3 = k3,
        tau3 = tau3,
        vars = vars
      )
      df$fatigue <-
        fatigue_fn(
          data = df,
          k3 = k3,
          tau3 = tau3,
          tau2 = tau2,
          vars = vars
        )
      df$predicted <-
        perf_model(
          data = df,
          P0 = P0,
          k1 = k1,
          k3 = k3,
          tau1 = tau1,
          tau2 = tau2,
          tau3 = tau3,
          vars = vars,
          target = target
        )

      # Compute RMSE, MAE and Rsquared on test data
      rmse_vec <- caret::RMSE(pred = df[, "predicted"],
                              obs = df[, "perf"])
      MAE_vec <- caret::MAE(pred = df[, "predicted"],
                            obs = df[, "perf"])
      Rsq_vec <- caret::R2(pred = df[, "predicted"],
                           obs = df[, "perf"])

      return(
        list(
          "data" = df,
          "theta" = theta,
          "rmse_vec" = rmse_vec,
          "MAE_vec" = MAE_vec,
          "Rsq_vec" = Rsq_vec
        )
      )
    }
  }



