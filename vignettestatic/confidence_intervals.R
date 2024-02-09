
library(ssdtools)
fit <- ssd_fit_dists(data = ssddata::ccme_silver)
set.seed = 99

t1 <- system.time(hc1 <- ssd_hc(fit, ci = TRUE, multi_est = FALSE, multi_ci = FALSE, weighted = FALSE))
t2 <- system.time(hc2 <- ssd_hc(fit, ci = TRUE, multi_est = TRUE, multi_ci = TRUE, weighted = FALSE))
t3 <- system.time(hc3 <- ssd_hc(fit, ci = TRUE, multi_est = TRUE, multi_ci = TRUE, weighted = TRUE))
t4 <- system.time(hc4 <- ssd_hc(fit, ci = TRUE, multi_est = FALSE, multi_ci = FALSE, weighted = TRUE))

compare_dat <- data.frame(
  method =c("weighted mean", "mixture (fixed weights)", "mixture", "weighted sample"),
  time = c(t1["elapsed"], t2["elapsed"], t3["elapsed"], t4["elapsed"]), 
  ucl = c(hc1$ucl, hc2$ucl, hc3$ucl, hc4$ucl),
  lcl = c(hc1$lcl, hc2$lcl, hc3$lcl, hc4$lcl),
  est = c(hc1$est, hc2$est, hc3$est, hc4$est))

save(hc1, hc2, hc3, hc4, t1, t2, t2, t4, compare_dat, file = "vignettestatic/confidence_intervals.RData")

