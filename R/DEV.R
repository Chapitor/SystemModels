
# data --------------------------------------------------------------------

example_data <- data.frame("training_load" = c(rnorm(100, mean = 1000, sd = 150),
                                               rep(0,50)),
                           "rest" = c(rnorm(100, mean= 2, sd=1),
                                      rep(1,50)),
                           "perf" = c(seq(from = 10, to = 30, length.out = 100),
                                      rep(0,50)),
                           "datetime" = seq(ISOdate(2020, 1, 1), by = "day", length.out = 150))


#' @title Dose-response modeling
#' @description This function trains a dose-response model with or without validation procedure. The model is defined as following :
#' @param data is a data frame object that contains at least training loads, performances and time between two consecutive sessions.
#' @param vars is a list that contains \emph{input} (i.e. session training loads) and \emph{time} (i.e. the time between two consecutive inputs) numerical vectors.
#' @param target is a character that indicates the performances column name.
#' @param date_ID is a character that indicates the date time object name.
#' @param specify default is \code{"NULL"}. Alternatively, a list of  \code{"theta_init} numeric vector that contains initial values for \emph{P0}, \emph{k1}, \emph{k3}, \emph{tau1}, \emph{tau2}, \emph{tau3} parameters,
#'  a numeric vector for lower bounds named \emph{lower}, a numeric vector for upper bounds named \emph{upper} and a character defining the method for optimisation \code{optim.method} has to be specified.
#'   \emph{PO} is the initial level of performance and can be extracted through the function [init_perf].
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
           validation.method = "none",
           specs = NULL) {
    require(tidyverse)
    require(optimx)
    require(caret)

    df <-
      data %>% dplyr::slice(-c(which(data[, target] == 0)))   # A reduced data frame that allow to model with no performance days. Used for datetime indexed data frames in Time series CV

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
      # Optimization is initial values and boundaries are provided
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
      rmse_vec <- caret::RMSE(pred = data[which(data[,target] != 0), "predicted"],
                              obs = data[which(data[,target] != 0), "perf"])
      MAE_vec <- caret::MAE(pred = data[which(data[,target] != 0), "predicted"],
                            obs = data[which(data[,target] != 0), "perf"])
      Rsq_vec <- caret::R2(pred = data[which(data[,target] != 0), "predicted"],
                           obs = data[which(data[,target] != 0), "perf"])



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
        initialWindow = specs[["initialWindow"]] * nrow(df),
        horizon = specs[["horizon"]] * nrow(df),
        fixedWindow = specs[["fixedWindow"]]
      )

      # split performance dataframe
      folder_train <- df[unlist(time_slice$train), ]
      folder_test <- df[unlist(time_slice$test), ]

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
        folder_train <- df[unlist(time_slice$train[k]), ]
        folder_test <- df[unlist(time_slice$test[k]), ]

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




# TEST TS-CV --------------------------------------------------------------------
data <- example_data
target <- "perf"
vars <- list("input" = example_data$training_load, "time" = example_data$rest)
date_ID <- "datetime"


P0_init <- init_perf(data, target)
k1_init <- 0.5
k3_init <- 0.1
tau1_init <- 40
tau2_init <- 20
tau3_init <- 5

theta_init <- c(P0_init, k1_init, k3_init, tau1_init, tau2_init, tau3_init)
lower <- c(P0_init - 0.10 * P0_init, 0, 0, 10, 1, 1)
upper <- c(P0_init, 1, 1, 80, 40, 10)

specs <- list("initialWindow" = 50, "horizon" = 15, "fixedWindow" = FALSE)
validation.method = "TS-CV"
specify = list("theta_init" = theta_init, "lower" = lower, "upper" = upper, "optim.method" = "nlm")


a <- sysmod(data = example_data,
            vars = vars,
            target = target,
            date_ID = date_ID,
            specify = specify, validation.method = validation.method, specs = specs)


# TEST simple CV ----------------------------------------------------------
specs <- list("initialWindow" = 0.8, "horizon" = 0.2, "fixedWindow" = FALSE)
validation.method = "simple"

a <- sysmod(data = example_data,
            vars = vars,
            target = target,
            date_ID = date_ID,
            specify = specify, validation.method = validation.method, specs = specs)


# TEST no CV --------------------------------------------------------------
specs <- NULL
validation.method = "none"

a <- sysmod(data = example_data,
            vars = vars,
            target = target,
            date_ID = date_ID,
            specify = specify,
            validation.method = validation.method,
            specs = specs)








# # simple CV
# validation.method = "simple"
# specs <- list("p" = list("initialWindow" = 0.8, "horizon" = 0.2, "fixedWindow" = FALSE))
# date_ID <- "datetime"



