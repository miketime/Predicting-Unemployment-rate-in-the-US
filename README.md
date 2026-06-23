# US Unemployment Rate Forecasting with ARIMA

Time series analysis research paper building and forecasting an ARIMA model for the US unemployment rate (1990–2025), including a from-scratch implementation of the Durbin-Levinson algorithm for PACF computation.

**Author:** Isfan Mihai-George
**Professor:** Cidota Marina Anca

## Overview

This project models and forecasts the monthly US unemployment rate using Box-Jenkins ARIMA methodology. It covers:

- Stationarity analysis (differencing, ACF/PACF inspection)
- A custom implementation of the **Durbin-Levinson algorithm** to compute partial autocorrelations and recursively estimate AR model order, validated against R's built-in `pacf()` on simulated AR(2) data
- Manual ARIMA order selection vs. `auto.arima()`, compared by AIC
- 20-step-ahead forecasting with 95% confidence intervals, benchmarked against the Federal Reserve's own unemployment projections

## Dataset

- **Source:** [FRED – Federal Reserve Bank of St. Louis, UNRATE series](https://fred.stlouisfed.org/series/UNRATE)
- **Range:** April 1990 – April 2025, monthly, seasonally adjusted
- **Observations:** 421, no missing values
- Summary stats: min 3.4%, median 5.4%, mean 5.7%, max 14.8% (COVID-19 spike, 2020)

## Methodology

1. **EDA & stationarity** — original series shows slow-decaying ACF (non-stationary); first-order differencing produces rapid ACF/PACF cutoff after lag 1–2, confirming stationarity.
2. **Durbin-Levinson algorithm** — implemented from scratch (`autocovariance()`, `get_autocovariances()`, `durbin_levinson()`) to recursively solve for PACF coefficients. Validated on simulated AR(2) data (φ₁=1, φ₂=-0.9, n=200): max absolute difference vs. R's `pacf()` < 0.001.
3. **Order selection** — PACF of differenced data significant at lags 1–2 (threshold ±1.96/√n), suggesting AR(2) structure.
4. **Model comparison** (by AIC):

   | Model | AIC |
   |---|---|
   | ARIMA(1,1,0) | 707.32 |
   | **ARIMA(2,1,0)** | **703.70** |
   | ARIMA(3,1,0) | 705.65 |
   | ARIMA(1,1,1) | 704.69 |
   | auto.arima | 705.48 |

   **ARIMA(2,1,0)** selected as best fit (lowest AIC).
5. **Forecasting** — 20-step-ahead forecast from ARIMA(2,1,0), compared against `auto.arima()` forecast (mean absolute difference ≈ 0.011) and against the Fed's own FOMC projection (~4.3% through 2027).

## Results

- Manual ARIMA(2,1,0) and auto.arima forecasts converge to a stable ~4.19% unemployment rate, closely tracking the Fed's projected 4.3% stagnation through 2027.
- Confirms ARIMA(2,1,0) as an adequate, parsimonious model for this series.

## Repository Structure

```
.
├── SDS_Project_Isfan_Mihai_George_DS_411.R   # Main analysis script (EDA, Durbin-Levinson, ARIMA, forecasting)
├── Script_comparare_rezultate.R              # Extended analysis: ADF/KPSS tests, exhaustive auto.arima search,
│                                              # residual diagnostics (Ljung-Box, ARCH, normality), out-of-sample CV
├── UNRATE.csv                                # FRED unemployment rate dataset
├── SDSProject1_DS411_Isfan_Mihai_George.pdf  # Full written report
├── renv.lock / renv/                         # R environment lockfile for reproducibility
└── SDS_Project_Isfan_Mihai_George_DS_411.Rproj
```

## Requirements

R packages: `forecast`, `tseries`, `ggplot2`, `gridExtra`, `tidyquant`

```r
install.packages(c("forecast", "tseries", "ggplot2", "gridExtra", "tidyquant"))
```

Or restore the exact environment via `renv`:

```r
renv::restore()
```

## Usage

```r
# from the project directory
source("SDS_Project_Isfan_Mihai_George_DS_411.R")
```

## References

- Brockwell, P. J., & Davis, R. A. (1987). *Time series: Theory and methods*. Springer-Verlag.
- Shumway, R. H., & Stoffer, D. S. (2017). *Time series analysis and its applications: With R examples* (4th ed.). Springer.
- Federal Reserve Bank of St. Louis. (2025). *Unemployment rate* [Data set]. FRED. https://fred.stlouisfed.org/series/UNRATE

  # THIS IS A FINISHED PROJECT AND NO FURTHER COMMITS WILL BE MADE
