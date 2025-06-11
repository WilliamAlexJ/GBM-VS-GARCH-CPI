# GBM-VS-GARCH-CPI
Comparing GBM to GARCH (1, 1) using CPI data gathered from FRED between 2010-01-01 -> 2018-12-01

This project demonstrates how to simulate and compare forward forecasts of the U.S. Consumer Price Index (CPI) using two approaches in R:

Geometric Brownian Motion (GBM) Monte Carlo simulation on log-returns.

GARCH(1,1) modeling and path simulation via the rugarch package.


RESULTS
After running 200 Monte Carlo simulations over a 10-year horizon:

GBM Forecast
The CPI is projected to drift steadily upward with a mean level around 256.3 by the end of the period. The 90% confidence band spans roughly 253.6 to 259.2, reflecting moderate uncertainty under purely diffusive dynamics.

GARCH(1,1) Forecast
Incorporating time-varying volatility yields a very similar average outcome (mean ≈ 256.3) but slightly tighter dispersion (90% band ~ 254.3 to 258.2), thanks to volatility clustering in the log-returns.
Overall, both models suggest a continuation of the recent slow but persistent rise in CPI, with the GARCH framework offering marginally more confidence in the middle‐range scenarios.

PLOT
![GBMVSGARCHCPIplot](https://github.com/user-attachments/assets/4f447629-115c-49b2-af8a-5d38c5818201)
