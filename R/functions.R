# rm(list = ls())
# library(usethis)
# library(devtools)
# library(roxygen2)

# Initiate an example data frame
# load("./data/example_data.rda")

#' @title  Adaptations response modeling
#' @description To model adaptation induced by each training session.
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param k1 denotes the gain term for adaptation occurrences.
#' @param tau1 denotes the time constant of the adaptations exponential decrease.
#' @param vars is a list that contains \emph{input} (i.e. session training loads) and \emph{time} (i.e. the time between two consecutive inputs) numerical vectors.
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


#' @title Variation in the fatiguing effects.
#' @description To make the fatigue impulse response varying over time in accordance with system input.
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param k3 denotes the gain term for fatigue occurrences.
#' @param tau3 denotes the time constant of the fatigue remanence exponential decrease.
#' @param vars is a list that contains \emph{input} (i.e. session training loads) and \emph{time} (i.e. the time between two consecutive inputs) numerical vectors.
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


#' @title Fatigue response modeling.
#' @description To calculate the variable fatigue component for each training input.
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param k3 denotes the gain term for fatigue occurrences.
#' @param tau2 denotes the time constant of fatigue exponential decrease.
#' @param tau3 denotes the time constant of the fatigue remanence exponential decrease.
#' @param vars is a list that contains \emph{input} (i.e. session training loads) and \emph{time} (i.e. the time between two consecutive inputs) numerical vectors.
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
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param target is a character that indicates the performances column name.
#'
#' @return a numerical value
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
#' @description This function extracts known performances from a data frame object.
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
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
#' @description The function models the performance according to Busso Variable dose-response model \(Busso, 2003\) and based on previously defined parameters.
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param P0 denotes the basic level of performance.
#' @param k1 denotes the gain term for adaptations response.
#' @param k3 denotes the gain term for fatigue response.
#' @param tau1 denotes the time constant of the adaptations exponential decrease.
#' @param tau2 denotes the time constant of fatigue exponential decrease.
#' @param tau3 denotes the time constant of the fatigue remanence exponential decrease.
#' @param vars is a list that contains \emph{input} (i.e. session training loads) and \emph{time} (i.e. the time between two consecutive inputs) numerical vectors.
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



#' @title Residual sum of squares.
#' @description This function calculates the residual sum of squares from observations and predicted performances.
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param theta is a vector of parameters \emph{P0}, \emph{k1}, \emph{k3}, \emph{tau1}, \emph{tau2}, \emph{tau3}.
#' @param target is a character that indicates the performance column name.
#' @param vars is a list that contains \emph{input} (i.e. session training loads) and \emph{time} (i.e. the time between two consecutive inputs) numerical vectors.
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
#' @description This function trains a dose-response model \insertCite{busso2003variable}{sysmod} with or without validation procedure. The model is defined as following :
#' @param data A data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param vars A list that contains \emph{input} (i.e. session training loads) and \emph{time} (i.e. the time between two consecutive inputs) numerical vectors.
#' @param target A character that indicates the performances column name.
#' @param date_ID A character that indicates the date time object name.
#' @param specify default is \code{"NULL"}. Alternatively, a list of \code{"theta_init"} numeric vector that contains initial values for \emph{P0}, \emph{k1}, \emph{k3}, \emph{tau1}, \emph{tau2}, \emph{tau3} parameters,
#'  a numeric vector for lower bounds named \emph{lower}, a numeric vector for upper bounds named \emph{upper} and a character defining the method for optimisation \code{optim.method} has to be specified.
#'   \emph{PO} denotes the initial level of performance. The first performance can be extracted through the function [init_perf].
#' @param validation.method default is \code{"none"}. Alternatively, data splitting or cross-validation can be specified (see details).
#' @param specs default is \code{"NULL"}. If a validation method is specified, a list that contains splitting arguments has to be specified (see details)
#'
#' @details \emph{validation.method} can be used to learn model and evaluate within a cross-validation procedure. Methods available are \code{"none"}, \code{"simple"} and \code{"TS-CV"}
#'  for no validation, a simple data split into training / testing sets and time-series cross-validation respectively.
#'  For the \code{"simple"} method, specify decimals as the proportion of the data to be used for model training in \code{"initialWindow"}, the proportion of data used for model evaluation in
#'  \code{"horizon"} and logical term for the \code{"fixedWindow"} within a list.
#'  For the \code{"TS-CV"} method, specify numeric values for \code{"initialWindow"}, \code{"horizon"} and logical term for the \code{"fixedWindow"} within a list.
#' Each of the optimization algorithm with constraints and used by the optimx package can be specified in optim.method.
#'
#' @return A list describing the model output and its performances.
#'
#' @note Model performances (RMSE, MAE and R squared) are calculated on test data for validation = c("simple", "TS-CV").
#'
#' @examples
#' #' DO NOT RUN : no validation, default optimisation specs.
#' model_results <- sysmod(data = example_data,
#'       vars = list("input" = example_data$training_load, "time" = example_data$rest),
#'       target = "perf", date_ID = "datetime",
#'       validation.method = "none")
#'
#' P0_init = init_perf(data = example_data, target = "perf")
#' theta_init <- c(P0_init = P0_init, k1_init = 0.5, k3_init = 0.1, tau1_init = 40, tau2_init = 20, tau3_init = 5)
#' lower <- c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1)
#' upper <- c(P0_init, 1, 1, 80, 40, 10)
#' DO NOT RUN : no validation, custom optimisation.
#'model_results <- sysmod(data = example_data,
#'      vars = list("input" = example_data$training_load, "time" = example_data$rest),
#'      target = "perf", date_ID = "datetime",
#'      specify = list("theta_init" = theta_init, "lower" = lower, "upper" = upper, "optim.method" = "nlm"),
#'      validation.method = "none")
#'
#' DO NOT RUN : simple split example, custom optimisation.
#' model_results <- sysmod(data = example_data,
#'     vars = list("input" = example_data$training_load, "time" = example_data$rest),
#'     target = "perf", date_ID = "datetime",
#'     specify = list("theta_init" = theta_init, "lower" = lower, "upper" = upper, "optim.method" = "nlm"),
#'     validation.method = "simple",
#'     specs = list("initialWindow" = 0.8, "horizon" = 0.2, "fixedWindow" = FALSE))
#'
#' DO NOT RUN : TS-CV example, custom optimisation.
#'model_results <- sysmod(data = example_data,
#'      vars = list("input" = example_data$training_load, "time" = example_data$rest),
#'      target = "perf", date_ID = "datetime",
#'      specify = list("theta_init" = theta_init, "lower" = lower, "upper" = upper, "optim.method" = "nlm"),
#'      validation.method = "TS-CV",
#'      specs = list("initialWindow" = 50, "horizon" = 15, "fixedWindow" = FALSE))
#'
#'@author Frank Imbach <frankimbach@gmail.com>
#'@references
#'\insertAllCited{}
#'@export
sysmod <-
  function(data,
           vars,
           target,
           date_ID,
           specify = NULL,
           validation.method = "none",
           specs = NULL) {
    require(tidyverse)
    require(optimx)
    require(caret)

    # Verify arguments
    if (names(vars[1]) != "input" | names(vars[2]) != "time") {
      stop("vars has to be a list that contains `ìnput` and `time` numerical objects ")
    }
    if (is.character(target) == FALSE) {
      stop("`target` has to be a character")
    }
    if (is.character(date_ID) == FALSE) {
      stop("`date_ID` has to be a character")
    }
    if (validation.method != "none" & validation.method != "simple" & validation.method != "TS-CV") {
      stop("validation.method is mis-specified")
    }
    if (names(specify[1]) != "theta_init" | names(specify[2]) != "lower" | names(specify[3]) != "upper" | names(specify[4]) != "optim.method"){
      stop("specify should be a list that contain `theta_init`, `lower`, `ùpper`, `optim.method` named objects")
    }
    if (names(specs[1]) != "initialWindow" | names(specs[2]) != "horizon" | names(specs[3]) != "fixedWindow") {
      stop("specs should be a list that contain `ìnitialWindow`, `horizon` and `fixedWindow` objects")
    }

    if(is.na(which(data[, target] == 0)[1]) == FALSE) {
      df <-
        data %>% dplyr::slice(-c(which(data[, target] == 0)))   # A reduced data frame that allow to model with no performance days. Used for datetime indexed data frames in Time series CV
    } else {
      df <- data
    }

    # Initiate and optimize parameters
    if (is.null(specify) == FALSE) {
      theta_init <- specify[["theta_init"]]


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



    # No validation -----------------------------------------------------------


    if (validation.method == "none") {

      if (is.null(specify) == FALSE) {
        res_optim <-
          optimx::optimx(
            par = specify[["theta_init"]],
            fn = RSS,
            data = data,
            target = target,
            vars = vars,
            method = specify[["optim.method"]],
            lower = specify[["lower"]],
            upper = specify[["upper"]]
          )
      } else {
        res_optim <-
          optimx::optimx(
            par = theta_init,
            fn = RSS,
            data = data,
            target = target,
            vars = vars,
            method = "nlm",
            lower = c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1),
            upper = c(P0_init, 1, 1, 80, 40, 10)
          )
      }

      P0 <- res_optim[[1]]
      k1 <- res_optim[[2]]
      k3 <- res_optim[[3]]
      tau1 <- res_optim[[4]]
      tau2 <- res_optim[[5]]
      tau3 <- res_optim[[6]]

      theta <- c(
        P0 = P0,
        k1 = k1,
        k3 = k3,
        tau1 = tau1,
        tau2 = tau2,
        tau3 = tau3
      )

      data$perf <- real_perf(data = data, target = target)
      data$adaptation <-
        adaptation_fn(
          data = data,
          k1 = k1,
          tau1 = tau1,
          vars = vars
        )
      data$k2i <- k2i_fn(
        data = data,
        k3 = k3,
        tau3 = tau3,
        vars = vars
      )
      data$fatigue <-
        fatigue_fn(
          data = data,
          k3 = k3,
          tau3 = tau3,
          tau2 = tau2,
          vars = vars
        )
      data$predicted <-
        perf_model(
          data = data,
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
      rmse_vec <-
        caret::RMSE(pred = data[which(data[, target] != 0), "predicted"],
                    obs = data[which(data[, target] != 0), "perf"])
      MAE_vec <-
        caret::MAE(pred = data[which(data[, target] != 0), "predicted"],
                   obs = data[which(data[, target] != 0), "perf"])
      Rsq_vec <-
        caret::R2(pred = data[which(data[, target] != 0), "predicted"],
                  obs = data[which(data[, target] != 0), "perf"])



      return(
        list(
          "data" = data,
          "theta" = theta,
          "rmse_vec" = rmse_vec,
          "MAE_vec" = MAE_vec,
          "Rsq_vec" = Rsq_vec
        )
      )
    } # Close "No validation"


    # Simple data split for validation ----------------------------------


    if (validation.method == "simple") {
      time_slice <- caret::createTimeSlices(
        y = df[, target],
        initialWindow = round(specs[["initialWindow"]] * nrow(df)),
        horizon = round(specs[["horizon"]] * nrow(df)),
        fixedWindow = specs[["fixedWindow"]]
        )

      # split performance dataframe
      folder_train <- df[unlist(time_slice$train),]
      folder_test <- df[unlist(time_slice$test),]

      # retrieve date information of split
      datetime_train_min <- min(folder_train[, date_ID])
      datetime_train_max <- max(folder_train[, date_ID])

      datetime_test_min <- min(folder_test[, date_ID])
      datetime_test_max <- max(folder_test[, date_ID])

      # split datetime data frame
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
            par = specify[["theta_init"]],
            fn = RSS,
            data = folder_train,
            target = target,
            vars = vars,
            method = specify[["optim.method"]],
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
            method = "nlm",
            lower = c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1),
            upper = c(P0_init, 1, 1, 80, 40, 10)
          )
      }

      P0 <- res_optim[[1]]
      k1 <- res_optim[[2]]
      k3 <- res_optim[[3]]
      tau1 <- res_optim[[4]]
      tau2 <- res_optim[[5]]
      tau3 <- res_optim[[6]]

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

      # Compute RMSE, MAE and Rsquared on test data
      rmse_vec <- caret::RMSE(pred = folder_test[, "predicted"],
                              obs = folder_test[, "perf"])
      MAE_vec <- caret::MAE(pred = folder_test[, "predicted"],
                            obs = folder_test[, "perf"])
      Rsq_vec <- caret::R2(pred = folder_test[, "predicted"],
                           obs = folder_test[, "perf"])

      return(
        list(
          "data" = folder_test,
          "theta" = theta,
          "rmse_vec" = rmse_vec,
          "MAE_vec" = MAE_vec,
          "Rsq_vec" = Rsq_vec
        )
      )
    } # Close "simple"



    # TS-CV -------------------------------------------------------------------


    # Times Series Cross Validation
    if (validation.method == "TS-CV") {
      time_slice <- caret::createTimeSlices(
        df[, target],
        initialWindow = specs[["initialWindow"]],
        horizon = specs[["horizon"]],
        fixedWindow = specs[["fixedWindow"]]
      )

      for (k in 1:length(time_slice$train)) {
        # Compute the model for each folder. The last folder will return the full model

        # split performance dataframe
        folder_train <- df[unlist(time_slice$train[k]),]
        folder_test <- df[unlist(time_slice$test[k]),]

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
              par = specify[["theta_init"]],
              fn = RSS,
              data = folder_train,
              target = target,
              vars = vars,
              method = specify[["optim.method"]],
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
              method = "nlm",
              lower = c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1),
              upper = c(P0_init, 1, 1, 80, 40, 10)
            )
        }

        P0 <- res_optim[[1]]
        k1 <- res_optim[[2]]
        k3 <- res_optim[[3]]
        tau1 <- res_optim[[4]]
        tau2 <- res_optim[[5]]
        tau3 <- res_optim[[6]]

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


      } # Close TS-CV


      return(
        list(
          "dfs" = dfs,
          "theta" = theta_df,
          "rmse_vec" = rmse_vec,
          "MAE_vec" = MAE_vec,
          "Rsq_vec" = Rsq_vec
        )
      )
    }
  }



# test --------------------------------------------------------------------

# load("./data/example_data.rda")
# P0_init = init_perf(data = example_data, target = all_of("perf"))
# theta_init <- c(P0_init = P0_init, k1_init = 0.5, k3_init = 0.1, tau1_init = 40, tau2_init = 20, tau3_init = 5)
# lower <- c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1)
# upper <- c(P0_init, 1, 1, 80, 40, 10)
# model_results <- sysmod(data = example_data,
#                        vars = list("input" = example_data$training_load, "time" = example_data$rest),
#                               target = "perf", date_ID = "datetime",
#                               specify = list("theta_init" = theta_init, "lower" = lower, "upper" = upper, "optim.method" = "nlm"),
#                               validation.method = "simple",
#                               specs = list("initialWindow" = 0.8, "horizon" = 0.2, "fixedWindow" = FALSE))




