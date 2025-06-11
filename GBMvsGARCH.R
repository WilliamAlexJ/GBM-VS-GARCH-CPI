# CPI Level Simulation: Geometric Brownian Motion vs GARCH(1,1)
# Required packages:
install.packages(c("readr", "xts", "rugarch", "ggplot2", "dplyr", "tidyr", "patchwork"))

library(readr)
library(xts)
library(rugarch)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#Load CPI data

df <- read_csv("CPIAUCSL-3.csv")
# Ensure column names are DATE and CPI
names(df) <- c("DATE", "CPI")
# Parse dates using lubridate for robust handling
library(lubridate)
df$DATE <- ymd(df$DATE)
xts_cpi <- xts(df$CPI, order.by = df$DATE)

#Compute log-returns for GBM. Compute log-returns for GBM
logret <- diff(log(xts_cpi))[-1]
mu     <- mean(logret)
sigma  <- sd(logret)

#Simulation settings
h_years <- 10              
freq     <- frequency(xts_cpi)   
h        <- h_years * freq  
dt       <- 1 / freq        
last_val <- as.numeric(last(xts_cpi))
sims     <- 200           

#Dates for forecast
start_date <- index(last(xts_cpi)) + dt
dates_fut  <- seq(from = start_date, by = dt, length.out = h)

#Geometric Brownian Motion simulation
set.seed(42)
# Adjust drift for discrete simulation: mu_adj = mu - 0.5*sigma^2
drift_disc    <- (mu - 0.5 * sigma^2) * dt
vol_disc      <- sigma * sqrt(dt)
increments_gbm <- matrix(
  rnorm(h * sims, mean = drift_disc, sd = vol_disc),
  nrow = h, byrow = TRUE
)
log_paths     <- apply(increments_gbm, 2, cumsum) + log(last_val)
paths_gbm     <- exp(log_paths)

df_gbm <- as.data.frame(paths_gbm) %>%
  setNames(paste0("Sim", seq_len(sims))) %>%
  mutate(Date = dates_fut) %>%
  pivot_longer(-Date, names_to = "Sim", values_to = "Value") %>%
  group_by(Date) %>%
  mutate(
    Mean  = mean(Value),
    Lower = quantile(Value, probs = 0.025),
    Upper = quantile(Value, probs = 0.975)
  ) %>%
  ungroup()

#GARCH(1,1) simulation via ugarchsim
# Use the fitted model to simulate future log-returns
sim <- ugarchsim(fit, n.sim = h, m.sim = sims, rseed = 42)

# Extract simulated log-returns (matrix h × sims)
ret_sim_mat <- sim@simulation$seriesSim

# Build log-price paths and exponentiate
garch_log_paths <- apply(ret_sim_mat, 2, cumsum) + log(last_val)
paths_garch     <- exp(garch_log_paths)

df_garch <- as.data.frame(paths_garch) %>%
  setNames(paste0("Sim", seq_len(sims))) %>%
  mutate(Date = dates_fut) %>%
  pivot_longer(-Date, names_to = "Sim", values_to = "Value") %>%
  group_by(Date) %>%
  mutate(
    Mean  = mean(Value),
    Lower = quantile(Value, probs = 0.025),
    Upper = quantile(Value, probs = 0.975)
  ) %>%
  ungroup()

#Plot simulations side by side
# Subset a few paths for clarity
set.seed(123)
sample_sims <- paste0("Sim", sample(seq_len(sims), 20))
df_gbm_lines  <- df_gbm %>% filter(Sim %in% sample_sims)
df_garch_lines <- df_garch %>% filter(Sim %in% sample_sims)

plot_common <- theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

gbm_plot <- ggplot() +
  # 95% Ribbon in light blue
  geom_ribbon(data = df_gbm, aes(x = Date, ymin = Lower, ymax = Upper),
              fill = "lightblue", alpha = 0.3) +
  # A subset of simulation paths in light blue
  geom_line(data = df_gbm_lines, aes(x = Date, y = Value, group = Sim),
            color = "lightblue", size = 0.5, alpha = 0.5) +
  # Mean path in dark blue
  geom_line(data = df_gbm, aes(x = Date, y = Mean),
            color = "blue", size = 1) +
  labs(title = "GBM: Simulated CPI", x = "", y = "CPI Level") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  plot_common

garch_plot <- ggplot() +
  geom_ribbon(data = df_garch, aes(x = Date, ymin = Lower, ymax = Upper),
              fill = "lightblue", alpha = 0.3) +
  geom_line(data = df_garch_lines, aes(x = Date, y = Value, group = Sim),
            color = "lightblue", size = 0.5, alpha = 0.5) +
  geom_line(data = df_garch, aes(x = Date, y = Mean),
            color = "blue", size = 1) +
  labs(title = "GARCH(1,1): Simulated CPI", x = "", y = "") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  plot_common

# Combine side by side
ggb <- gbm_plot + garch_plot + plot_layout(ncol = 2, guides = "collect")
print(ggb)









#Extract summary metrics from simulations
gbm_metrics <- df_gbm %>%
  group_by(Date) %>%
  summarise(
    Mean       = unique(Mean),
    SD         = sd(Value),
    Median     = quantile(Value, 0.5),
    P05        = quantile(Value, 0.05),
    P95        = quantile(Value, 0.95)
  )

garch_metrics <- df_garch %>%
  group_by(Date) %>%
  summarise(
    Mean       = unique(Mean),
    SD         = sd(Value),
    Median     = quantile(Value, 0.5),
    P05        = quantile(Value, 0.05),
    P95        = quantile(Value, 0.95)
  )

# Display final-horizon summary
final_date <- tail(gbm_metrics$Date, 1)
final_gbm   <- filter(gbm_metrics, Date == final_date)
final_garch <- filter(garch_metrics, Date == final_date)

cat("Final horizon (", final_date, ") metrics:
", sep = "")
cat("GBM -> Mean:", round(final_gbm$Mean, 2),
    ", SD:", round(final_gbm$SD, 2),
    ", 5th pct:", round(final_gbm$P05, 2),
    ", 95th pct:", round(final_gbm$P95, 2), "
")
cat("GARCH -> Mean:", round(final_garch$Mean, 2),
    ", SD:", round(final_garch$SD, 2),
    ", 5th pct:", round(final_garch$P05, 2),
    ", 95th pct:", round(final_garch$P95, 2), "
")
