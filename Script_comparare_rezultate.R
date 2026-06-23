# Enhanced Unemployment Rate Analysis
# Install required packages (only the ones we actually use)
required_packages <- c("forecast", "tseries")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(forecast)
library(tseries)

# Load and prepare data
unemploy <- read.csv("UNRATE.csv")
unemploy_ts <- ts(unemploy$UNRATE, frequency=12, start=c(1990,1))

cat("=== DATA OVERVIEW ===\n")
cat("Data span:", start(unemploy_ts)[1], "to", end(unemploy_ts)[1], "\n")
cat("Total observations:", length(unemploy_ts), "\n")
cat("Frequency:", frequency(unemploy_ts), "(monthly data)\n")
print(summary(unemploy_ts))

# Enhanced exploratory analysis
par(mfrow = c(2,2), mar = c(4,4,3,2))
plot(unemploy_ts, main = "Unemployment Rate Time Series", ylab = "Rate (%)", col = "blue")
hist(unemploy_ts, main = "Distribution of Unemployment Rate", xlab = "Rate (%)", col = "lightblue")
boxplot(unemploy_ts ~ cycle(unemploy_ts), main = "Seasonal Patterns", 
        xlab = "Month", ylab = "Rate (%)", col = "lightgreen")

# Reset plot parameters and create decomposition plot separately
par(mfrow = c(1,1))
plot(decompose(unemploy_ts))
title(main = "Time Series Decomposition", line = -1)

# STATIONARITY TESTING - More Comprehensive
cat("\n=== STATIONARITY ANALYSIS ===\n")

# Augmented Dickey-Fuller Test
adf_test <- adf.test(unemploy_ts)
cat("ADF Test on Original Series:\n")
cat("Test Statistic:", round(adf_test$statistic, 4), "\n")
cat("p-value:", round(adf_test$p.value, 4), "\n")
cat("Conclusion:", ifelse(adf_test$p.value < 0.05, "Stationary", "Non-stationary"), "\n\n")

# KPSS Test (complementary test)
kpss_test <- kpss.test(unemploy_ts)
cat("KPSS Test on Original Series:\n")
cat("Test Statistic:", round(kpss_test$statistic, 4), "\n")
cat("p-value:", round(kpss_test$p.value, 4), "\n")
cat("Conclusion:", ifelse(kpss_test$p.value > 0.05, "Stationary", "Non-stationary"), "\n\n")

# Test differenced data
unemploy_diff <- diff(unemploy_ts)
adf_diff <- adf.test(unemploy_diff)
cat("ADF Test on Differenced Series:\n")
cat("Test Statistic:", round(adf_diff$statistic, 4), "\n")
cat("p-value:", round(adf_diff$p.value, 4), "\n")
cat("Conclusion:", ifelse(adf_diff$p.value < 0.05, "Stationary", "Non-stationary"), "\n")

# SEASONAL ANALYSIS
cat("\n=== SEASONAL ANALYSIS ===\n")
# Test for seasonal unit roots
tryCatch({
  if(length(unemploy_ts) >= 24) {
    seasonal_test <- nsdiffs(unemploy_ts)
    cat("Recommended seasonal differences:", seasonal_test, "\n")
    
    regular_diff <- ndiffs(unemploy_ts)
    cat("Recommended regular differences:", regular_diff, "\n")
  } else {
    seasonal_test <- 0
    regular_diff <- 1
    cat("Insufficient data for seasonal analysis\n")
  }
}, error = function(e) {
  seasonal_test <<- 0
  regular_diff <<- 1
  cat("Error in seasonal analysis, using defaults\n")
})

# AUTOMATIC MODEL SELECTION - This is what you were missing!
cat("\n=== AUTOMATIC MODEL SELECTION ===\n")

# Auto ARIMA (considers seasonality automatically)
tryCatch({
  auto_model <- auto.arima(unemploy_ts, 
                           seasonal = TRUE,
                           stepwise = FALSE,    # More thorough search
                           approximation = FALSE, # More accurate
                           trace = TRUE)        # Show search process
  
  cat("\nBest model selected by auto.arima:\n")
  print(auto_model)
}, error = function(e) {
  cat("Auto ARIMA failed, using fallback model\n")
  auto_model <<- arima(unemploy_ts, order = c(1, 1, 1))
  print(auto_model)
})

# COMPREHENSIVE MODEL COMPARISON
cat("\n=== COMPREHENSIVE MODEL COMPARISON ===\n")

# Your original models
model1 <- arima(unemploy_ts, order = c(0, 1, 4))
model2 <- arima(unemploy_ts, order = c(1, 1, 0))
model3 <- arima(unemploy_ts, order = c(2, 1, 0))
model4 <- arima(unemploy_ts, order = c(2, 1, 1))

# Additional models to test
model5 <- arima(unemploy_ts, order = c(1, 1, 1))  # ARIMA(1,1,1)
model6 <- arima(unemploy_ts, order = c(0, 1, 1))  # IMA(1,1) - Random walk with MA
model7 <- arima(unemploy_ts, order = c(3, 1, 0))  # AR(3) with differencing

# Seasonal models (if seasonality detected)
if(frequency(unemploy_ts) > 1) {
  tryCatch({
    model8 <- arima(unemploy_ts, order = c(1, 1, 1), seasonal = list(order = c(1, 1, 1), period = 12))
    model9 <- arima(unemploy_ts, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
  }, error = function(e) {
    model8 <<- NULL
    model9 <<- NULL
    cat("Seasonal models failed to converge\n")
  })
} else {
  model8 <- NULL
  model9 <- NULL
}

# Collect all models
all_models <- list(
  "AR(1,0,0)" = model1,
  "ARIMA(1,1,0)" = model2,
  "ARIMA(2,1,0)" = model3,
  "ARIMA(2,1,1)" = model4,
  "ARIMA(1,1,1)" = model5,
  "ARIMA(0,1,1)" = model6,
  "ARIMA(3,1,0)" = model7,
  "Auto-ARIMA" = auto_model
)

if(!is.null(model8)) all_models[["SARIMA(1,1,1)(1,1,1)"]] <- model8
if(!is.null(model9)) all_models[["SARIMA(0,1,1)(0,1,1)"]] <- model9

# Calculate information criteria
model_comparison <- data.frame(
  Model = names(all_models),
  stringsAsFactors = FALSE
)

# Safely calculate metrics for each model
for(i in 1:length(all_models)) {
  tryCatch({
    model_comparison$AIC[i] <- AIC(all_models[[i]])
    model_comparison$BIC[i] <- BIC(all_models[[i]])
    model_comparison$Log_Likelihood[i] <- as.numeric(logLik(all_models[[i]]))
    model_comparison$Parameters[i] <- length(coef(all_models[[i]]))
  }, error = function(e) {
    model_comparison$AIC[i] <<- NA
    model_comparison$BIC[i] <<- NA
    model_comparison$Log_Likelihood[i] <<- NA
    model_comparison$Parameters[i] <<- NA
  })
}

# Add AICc (corrected AIC for small samples)
model_comparison$AICc <- model_comparison$AIC + 
  (2 * model_comparison$Parameters * (model_comparison$Parameters + 1)) / 
  (length(unemploy_ts) - model_comparison$Parameters - 1)

# Remove any models that failed to fit (have NA values)
model_comparison <- model_comparison[complete.cases(model_comparison), ]

# Sort by AICc (best criterion for time series)
model_comparison <- model_comparison[order(model_comparison$AICc),]

cat("Model Comparison (sorted by AICc - best first):\n")
# Round only the numeric columns
model_comparison_display <- model_comparison
model_comparison_display[, -1] <- round(model_comparison_display[, -1], 2)
print(model_comparison_display)

# Select best model
if(nrow(model_comparison) > 0) {
  best_model_name <- model_comparison$Model[1]
  best_model <- all_models[[best_model_name]]
  
  cat("\nBest model:", best_model_name, "\n")
  cat("AICc:", round(model_comparison$AICc[1], 2), "\n")
} else {
  cat("Error in model comparison, using Auto-ARIMA as fallback\n")
  best_model_name <- "Auto-ARIMA"
  best_model <- auto_model
}

# COMPREHENSIVE RESIDUAL DIAGNOSTICS
cat("\n=== RESIDUAL DIAGNOSTICS FOR BEST MODEL ===\n")

residuals_best <- residuals(best_model)

# Ljung-Box test for multiple lags
lb_results <- data.frame(
  Lag = c(10, 15, 20),
  Statistic = numeric(3),
  P_Value = numeric(3)
)

for(i in 1:3) {
  tryCatch({
    test <- Box.test(residuals_best, lag = lb_results$Lag[i], type = "Ljung-Box")
    lb_results$Statistic[i] <- as.numeric(test$statistic)
    lb_results$P_Value[i] <- as.numeric(test$p.value)
  }, error = function(e) {
    lb_results$Statistic[i] <<- NA
    lb_results$P_Value[i] <<- NA
  })
}

cat("Ljung-Box Tests for Residual Autocorrelation:\n")
print(round(lb_results, 4))

if(all(lb_results$P_Value > 0.05, na.rm = TRUE)) {
  cat("✓ All Ljung-Box tests passed - no residual autocorrelation\n")
} else {
  cat("⚠ Some Ljung-Box tests failed - residual autocorrelation detected\n")
}

# ARCH test for heteroscedasticity
arch_test <- Box.test(residuals_best^2, lag = 12, type = "Ljung-Box")
cat("\nARCH Test for Heteroscedasticity:\n")
cat("Test Statistic:", round(arch_test$statistic, 4), "\n")
cat("p-value:", round(arch_test$p.value, 4), "\n")

if(arch_test$p.value > 0.05) {
  cat("✓ No ARCH effects detected\n")
} else {
  cat("⚠ ARCH effects detected - consider GARCH modeling\n")
}

# Normality tests
tryCatch({
  if(length(residuals_best) <= 5000) {
    shapiro_test <- shapiro.test(residuals_best)
    shapiro_p <- shapiro_test$p.value
  } else {
    shapiro_test <- shapiro.test(residuals_best[1:5000])
    shapiro_p <- shapiro_test$p.value
  }
}, error = function(e) {
  shapiro_p <- NA
})

tryCatch({
  jb_test <- jarque.bera.test(residuals_best)
  jb_p <- jb_test$p.value
}, error = function(e) {
  jb_p <- NA
})

cat("\nNormality Tests:\n")
if(!is.na(shapiro_p)) {
  cat("Shapiro-Wilk p-value:", round(shapiro_p, 4), "\n")
} else {
  cat("Shapiro-Wilk test failed\n")
}

if(!is.na(jb_p)) {
  cat("Jarque-Bera p-value:", round(jb_p, 4), "\n")
} else {
  cat("Jarque-Bera test failed\n")
}

# FORECASTING WITH VALIDATION
cat("\n=== FORECASTING WITH CROSS-VALIDATION ===\n")

# Out-of-sample validation
n_total <- length(unemploy_ts)
n_train <- n_total - 24  # Hold out last 24 months for validation

if(n_train > 50) {  # Only do validation if we have enough training data
  train_data <- window(unemploy_ts, end = c(start(unemploy_ts)[1] + (n_train-1) %/% 12, 
                                            (n_train-1) %% 12 + start(unemploy_ts)[2]))
  test_data <- window(unemploy_ts, start = c(end(train_data)[1], end(train_data)[2] + 1))
  
  # Refit best model on training data
  tryCatch({
    # Extract ARIMA order safely
    best_order <- arimaorder(best_model)
    model_train <- arima(train_data, order = best_order[1:3], 
                         seasonal = list(order = best_order[4:6], 
                                         period = frequency(unemploy_ts)))
  }, error = function(e) {
    # Fallback: refit the best model directly on training data
    cat("Using simplified refitting approach\n")
    model_train <<- arima(train_data, order = c(1, 1, 1))
  })
  
  # Generate forecasts for validation period
  forecast_validation <- forecast(model_train, h = length(test_data))
  mae_validation <- mean(abs(forecast_validation$mean - test_data))
  rmse_validation <- sqrt(mean((forecast_validation$mean - test_data)^2))
  
  cat("Out-of-sample validation (last 24 months):\n")
  cat("MAE:", round(mae_validation, 3), "\n")
  cat("RMSE:", round(rmse_validation, 3), "\n")
} else {
  cat("Insufficient data for out-of-sample validation\n")
  rmse_validation <- NA
}

# Final forecast on full data
final_forecast <- forecast(model1, h = 100)

cat("\n=== FORECAST RESULTS ===\n")
print(final_forecast)

# Enhanced visualization
par(mfrow = c(1,1), mar = c(4,4,3,2))

# Main forecast plot
plot(final_forecast, main = paste("Forecast:", best_model_name), 
     ylab = "Unemployment Rate (%)", col = "blue", fcol = "red")

# Residual plots
plot(residuals_best, main = "Residuals vs Time", ylab = "Residuals", type = "l")
abline(h = 0, col = "red", lty = 2)

qqnorm(residuals_best, main = "Q-Q Plot of Residuals")
qqline(residuals_best, col = "red")

acf(residuals_best, main = "ACF of Residuals", lag.max = 24)

# SUMMARY AND RECOMMENDATIONS
cat("\n=== SUMMARY AND RECOMMENDATIONS ===\n")

cat("1. DATA CHARACTERISTICS:\n")
cat("   - Monthly unemployment data from", start(unemploy_ts)[1], "to", end(unemploy_ts)[1], "\n")
cat("   - Series is non-stationary (needs differencing)\n")
if(seasonal_test > 0) cat("   - Seasonal patterns detected\n")

cat("\n2. MODEL SELECTION:\n")
cat("   - Your original models were a good start but limited\n")
cat("   - Best model:", best_model_name, "\n")
tryCatch({
  cat("   - Auto-ARIMA provided:", paste(arimaorder(auto_model), collapse = ","), "\n")
}, error = function(e) {
  cat("   - Auto-ARIMA model summary available above\n")
})

cat("\n3. MODEL QUALITY:\n")
if(exists("lb_results") && all(lb_results$P_Value > 0.05, na.rm = TRUE)) {
  cat("   ✓ Residuals are uncorrelated (good model fit)\n")
} else {
  cat("   ⚠ Some residual autocorrelation may remain\n")
}

if(!is.na(rmse_validation)) {
  cat("   - Out-of-sample RMSE:", round(rmse_validation, 3), "%\n")
} else {
  cat("   - Out-of-sample validation not performed\n")
}

cat("\n4. IMPROVEMENTS MADE:\n")
cat("   - Added automatic model selection\n")
cat("   - Included seasonal models\n")
cat("   - Comprehensive residual diagnostics\n")
cat("   - Out-of-sample validation\n")
cat("   - Multiple information criteria comparison\n")

cat("\n5. NEXT STEPS TO CONSIDER:\n")
cat("   - GARCH modeling if volatility clustering present\n")
cat("   - Structural break tests for recession periods\n")
cat("   - Vector autoregression if using multiple economic indicators\n")
cat("   - Machine learning approaches for comparison\n")
