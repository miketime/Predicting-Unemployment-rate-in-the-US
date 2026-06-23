#install.packages(c("forecast","tseries","ggplot2","gridExtra","tidyquant"))
library(forecast)

somaj <- read.csv("UNRATE.csv")
somaj_ts <- ts(somaj$UNRATE, frequency=12, start=c(1990,1))

plot(somaj_ts, main = 'Unemployment rate in the US between 1990-2025', ylab = 'unemployment rate %', xlab = 'year')
summary(somaj_ts)
length(somaj_ts)
sum(is.na(somaj_ts))

somaj_diff <- diff(somaj_ts)

par(mfrow = c(2,1), mar = c(2,4,4,2))
plot(somaj_ts, main = 'unemployment data')
plot(somaj_diff, main = 'unemployment data')

par(mfrow = c(2,1), mar = c(2,3,3,2))
acf(somaj_diff, lag.max=40, main = "ACF of differenced data")
pacf(somaj_diff, lag.max=40, main = "PACF of differenced data")
par(mfrow = c(1,1))

autocovariance <- function(x, h) {
  n <- length(x)
  if (h >= n) {
    return(0)
  }
  if (h < 0) {
    h <- abs(h)  
  }
  x_mean <- mean(x)
  if (h == 0) {
    gamma_h <- sum((x - x_mean)^2) / n
  } else {
    gamma_h <- sum((x[1:(n-h)] - x_mean) * (x[(h+1):n] - x_mean)) / n
  }
  return(gamma_h)
}

get_autocovariances <- function(x, max_lag) {
  gammas <- numeric(max_lag + 1)
  for (h in 0:max_lag) {
    gammas[h + 1] <- autocovariance(x, h)
  }
  return(gammas)
}

durbin_levinson <- function(gamma_vec, max_order) {
  phi_matrix <- matrix(0, nrow = max_order, ncol = max_order)
  v_vec <- numeric(max_order + 1)
  v_vec[1] <- gamma_vec[1]  
  
  for (n in 1:max_order) {
    numerator <- gamma_vec[n + 1] 
    if (n > 1) {
      for (k in 1:(n-1)) {
        correction <- phi_matrix[n-1, k] * gamma_vec[n + 1 - k]
        numerator <- numerator - correction
      }
    }
    phi_nn <- numerator / v_vec[n]
    phi_matrix[n, n] <- phi_nn
    
    if (n > 1) {
      for (k in 1:(n-1)) {
        phi_matrix[n, k] <- phi_matrix[n-1, k] - phi_nn * phi_matrix[n-1, n-k]
      }
    }
    v_vec[n + 1] <- v_vec[n] * (1 - phi_nn^2)
  }
  return(list(phi_matrix = phi_matrix, v_vec = v_vec))
}

set.seed(123)
ar2_data <- arima.sim(list(order = c(2, 0, 0), ar = c(1, -0.9)), n = 200)
plot(ar2_data, main = "AR(2) test data with parameters (1, -0.9)", type = "l")

ar2_gammas <- get_autocovariances(ar2_data, 10)
dl_ar2_result <- durbin_levinson(ar2_gammas, 10)
pacf_our_ar2 <- diag(dl_ar2_result$phi_matrix)

r_pacf_ar2 <- pacf(ar2_data, plot = FALSE, lag.max = 10)
r_pacf_values_ar2 <- as.numeric(r_pacf_ar2$acf)

comparison_ar2 <- data.frame(
  lag = 1:10,
  our_pacf = pacf_our_ar2,
  r_pacf = r_pacf_values_ar2,
  difference = abs(pacf_our_ar2 - r_pacf_values_ar2)
)
comparison_ar2

plot(1:10, pacf_our_ar2, type = "b", col = "red", main = "PACF Comparison: AR(2) Data", 
     ylab = "PACF", xlab = "Lag", ylim = c(-1, 1))
lines(1:10, r_pacf_values_ar2, type = "b", col = "blue")
legend("topright", c("DL Algorithm", "R PACF"), col = c("red", "blue"), lty = 1)

gamma_diff <- get_autocovariances(somaj_diff, 15)
dl_unemployment <- durbin_levinson(gamma_diff, 15)
pacf_unemployment <- diag(dl_unemployment$phi_matrix)

n_diff <- length(somaj_diff)
threshold <- 1.96 / sqrt(n_diff)

pacf_analysis <- data.frame(
  lag = 1:10,
  pacf_value = pacf_unemployment[1:10],
  threshold = threshold,
  significant = abs(pacf_unemployment[1:10]) > threshold
)
pacf_analysis

plot(1:10, pacf_unemployment[1:10], type = "b", col = "red", 
     main = "PACF Analysis for Unemployment Data (Differenced)", 
     ylab = "PACF", xlab = "Lag")
abline(h = c(-threshold, threshold), lty = 2, col = "blue")
abline(h = 0, col = "gray")

model1 <- arima(somaj_ts, order = c(1, 1, 0))
model2 <- arima(somaj_ts, order = c(2, 1, 0))
model3 <- arima(somaj_ts, order = c(3, 1, 0))
model4 <- arima(somaj_ts, order = c(1, 1, 1))
model5 <- arima(somaj_ts, order = c(2, 1, 1))
auto_model <- auto.arima(somaj_ts)
print(model1)
print(model2)
print(model3)
print(model4)
print(auto_model)


model_perfect <- model2
summary(model_perfect)
tsdiag(model_perfect)

forecasts <- predict(model_perfect, n.ahead = 20)
forecast_values <- forecasts$pred
forecast_se <- forecasts$se
lower_ci <- forecast_values - 1.96 * forecast_se
upper_ci <- forecast_values + 1.96 * forecast_se

plot(somaj_ts, xlim = c(2015, 2027), ylim = c(1, 12),
     main = "Unemployment Rate Forecasts",
     ylab = "Unemployment Rate (%)", xlab = "Year")

forecast_time <- seq(from = end(somaj_ts)[1] + (end(somaj_ts)[2] + 1)/12, 
                     by = 1/12, length.out = 20)
lines(forecast_time, forecast_values, col = "red", lwd = 2)
lines(forecast_time, lower_ci, col = "red", lty = 2)
lines(forecast_time, upper_ci, col = "red", lty = 2)
legend("topleft", c("Historical Data", "Forecast", "95% CI"), 
       col = c("black", "red", "red"), lty = c(1, 1, 2), lwd = c(1, 2, 1))

auto_forecasts <- forecast(auto_model, h = 20)
plot(somaj_ts, xlim = c(2015, 2027), ylim = c(1, 12),
     main = "Forecast Comparison: Manual vs auto.arima",
     ylab = "Unemployment Rate (%)", xlab = "Year")
lines(forecast_time, forecast_values, col = "red", lwd = 2)
lines(forecast_time, as.numeric(auto_forecasts$mean), col = "blue", lwd = 2)
legend("topleft", c("Historical data", "Manual model", "auto.arima"), 
       col = c("black", "red", "blue"), lwd = c(1, 2, 2))

forecast_comparison <- data.frame(
  period = 1:20,
  manual_forecast = as.numeric(forecast_values),
  auto_forecast = as.numeric(auto_forecasts$mean),
  difference = abs(as.numeric(forecast_values) - as.numeric(auto_forecasts$mean))
)
forecast_comparison

