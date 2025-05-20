library(bsts)
library(lubridate)
library(dplyr)
library(ggplot2)
library(forecast)
library(zoo)
library(car) 

data <- read.csv("bitcoin_processed_with_features.csv", header = TRUE)
data$date <- as.Date(data$date, format = "%m/%d/%Y")

data$close_log <- log(data$close)

split_date <- as.Date("2023-07-01")
train_data <- data[data$date < split_date, ]
test_data <- data[data$date >= split_date, ]

y_ts <- ts(train_data$close, frequency = 365)
decomp <- stl(y_ts, s.window = 'periodic')
plot(decomp)

xreg_vars <- names(train_data %>% select(-c(close, date, close_log)))

xreg_train <- train_data[, xreg_vars]
xreg_train <- xreg_train[, colSums(is.na(xreg_train)) == 0]

remaining_vars <- xreg_vars
max_retries <- length(xreg_vars)

for (i in 1:max_retries) {
  formula_str <- paste("close_log ~", paste(remaining_vars, collapse = "+"))
  lm_model <- lm(formula_str, data = train_data)
  
  alias_check <- alias(lm_model)
  if (!is.null(alias_check$Complete)) {
    aliased_vars <- rownames(alias_check$Complete)
    remaining_vars <- setdiff(remaining_vars, aliased_vars)
    message("Removed aliased variables: ", paste(aliased_vars, collapse=", "))
    next
    library(bsts)
    library(lubridate)
    library(dplyr)
    library(ggplot2)
    library(forecast)
    library(zoo)
    library(car)  
    
    data <- read.csv("bitcoin_processed_with_features - Copy.csv", header = TRUE)
    data$date <- as.Date(data$date, format = "%m/%d/%Y")
    
    data$close_log <- log(data$close)
    
    split_date <- as.Date("2023-07-01")
    train_data <- data[data$date < split_date, ]
    test_data <- data[data$date >= split_date, ]
    

    y_ts <- ts(train_data$close, frequency = 365)
    decomp <- stl(y_ts, s.window = 'periodic')
    plot(decomp)
    
    xreg_vars <- names(train_data %>% select(-c(close, date, close_log)))
    
    xreg_train <- train_data[, xreg_vars]
    xreg_train <- xreg_train[, colSums(is.na(xreg_train)) == 0]
    
    remaining_vars <- xreg_vars
    max_retries <- length(xreg_vars)
    
    for (i in 1:max_retries) {
      formula_str <- paste("close_log ~", paste(remaining_vars, collapse = "+"))
      lm_model <- lm(formula_str, data = train_data)
      
      alias_check <- alias(lm_model)
      if (!is.null(alias_check$Complete)) {
        aliased_vars <- rownames(alias_check$Complete)
        remaining_vars <- setdiff(remaining_vars, aliased_vars)
        message("Removed aliased variables: ", paste(aliased_vars, collapse=", "))
        next
      }
      
      vif_values <- tryCatch(
        vif(lm_model),
        error = function(e) {
          message("VIF error: ", e$message)
          return(NULL)
        }
      )
      
      if (is.null(vif_values)) break
      
      if (all(vif_values <= 5)) break
      
      worst_offender <- names(which.max(vif_values))
      remaining_vars <- remaining_vars[remaining_vars != worst_offender]
      message("Removed high VIF variable: ", worst_offender)
    }
    
    if (length(remaining_vars) == 0) stop("No variables remaining after feature selection!")
    
    xreg_train <- as.matrix(train_data[, remaining_vars])
    xreg_test <- as.matrix(test_data[, remaining_vars])
    
    xreg_means <- colMeans(xreg_train, na.rm = TRUE)
    xreg_sds <- apply(xreg_train, 2, sd, na.rm = TRUE)
    xreg_train_scaled <- scale(xreg_train, center = xreg_means, scale = xreg_sds)
    xreg_test_scaled <- scale(xreg_test, center = xreg_means, scale = xreg_sds)
    
    pacf(train_data$close_log, main = "PACF of Training Data", lag.max = 5)
    ar_lags <- 1  
    
    ss <- AddStudentLocalLinearTrend(list(), train_data$close_log)
    
    if (sd(decomp$time.series[, "seasonal"]) > 0.1*sd(train_data$close)) {
      ss <- AddSeasonal(ss, train_data$close_log, season.duration = 7 ,nseasons = 52)
    }
    
    ss <- AddAr(ss, train_data$close_log, lags = ar_lags)
    
    model <- bsts(
      train_data$close_log,
      state.specification = ss,
      niter = 5000,
      ping = 0,
      seed = 54321,
      expected.model.size = 5,
      xreg = xreg_train_scaled
    )
    
    burn <- SuggestBurn(0.1, model)
    
    par(mfrow = c(3,1))
    for(i in 1:3) {
      plot(model$coefficients[,i], type = 'l', 
           main = paste("Trace Plot:", colnames(xreg_train_scaled)[i]))
    }
    
    residuals <- colMeans(model$one.step.prediction.errors[-(1:burn),])
    par(mfrow = c(1,2))
    acf(residuals, main = "ACF of Residuals")
    qqnorm(residuals, main = "Q-Q Plot of Residuals")
    
    n_test <- nrow(test_data)
    forecasts <- numeric(n_test)
    forecast_dates <- numeric(n_test)
    
    for(i in 1:n_test) {
      current_date <- split_date + days(i-1)
      cat("\rProcessing:", as.character(current_date), " ", i, "/", n_test)
      
      current_train <- data[data$date <= current_date, ]
      y_current <- current_train$close_log
      
      xreg_current <- as.matrix(current_train[, remaining_vars])
      xreg_current_scaled <- scale(xreg_current, 
                                   center = xreg_means, 
                                   scale = xreg_sds)
      
      current_model <- bsts(
        y_current,
        state.specification = ss,
        niter = 2000,
        ping = 0,
        xreg = xreg_current_scaled
      )
      
      p <- predict(current_model, horizon = 1, 
                   newdata = xreg_test_scaled[i,,drop=FALSE], 
                   burn = 500)
      forecasts[i] <- exp(p$mean)  
      forecast_dates[i] <- as.numeric(current_date)
    }
    
    actual_values <- test_data$close
    forecast_values <- forecasts
    
    metrics <- list(
      MAPE = mean(abs(actual_values - forecast_values)/actual_values),
      RMSE = sqrt(mean((actual_values - forecast_values)^2)),
      MAE = mean(abs(actual_values - forecast_values)),
      Directional = mean(
        (actual_values > test_data$close[1]) == (forecast_values > test_data$close[1])
      )
      
      cat("\n==== Final Model Performance ====\n")
      sprintf("MAPE: %.2f%%\n", metrics$MAPE*100))
    sprintf("RMSE: %.2f\n", metrics$RMSE))
cat(sprintf("MAE: %.2f\n", metrics$MAE))
cat(sprintf("Directional Accuracy: %.2f%%\n", metrics$Directional*100))

forecast_df <- data.frame(
  Date = as.Date(forecast_dates, origin = "2023-07-01"),
  Forecast = forecast_values,
  Actual = actual_values
)

ggplot(forecast_df, aes(x = Date)) +
  geom_line(aes(y = Actual, color = "Actual"), size = 1) +
  geom_line(aes(y = Forecast, color = "Forecast"), size = 1, linetype = 2) +
  geom_ribbon(aes(ymin = Forecast*0.95, ymax = Forecast*1.05), 
              fill = "blue", alpha = 0.2) +
  scale_color_manual(values = c("Actual" = "red", "Forecast" = "blue")) +
  labs(title = "Bitcoin Price Forecast with Enhanced BSTS Model",
       subtitle = paste("MAPE:", round(metrics$MAPE*100, 1), "%",
                        "Directional Accuracy:", round(metrics$Directional*100, 1), "%"),
       y = "Price (USD)") +
  theme_minimal() +
  theme(legend.position = "bottom")theme(legend.position = "bottom")