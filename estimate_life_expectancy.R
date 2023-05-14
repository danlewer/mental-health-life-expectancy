# -----------------------
# libraries and functions
# -----------------------

library(data.table)
library(MASS) # for MASS::mvrnorm (sample from multivariate normal distribution)

# vectorized poisson confidence interval
vpt <- function(x, t, prefix = '', ...) {
  f <- function(xl, tl) {
    if (is.na(xl) | is.na(tl)) return(c(0, 0, 0))
    y <- poisson.test(xl, tl, ...)
    c(xl/tl, y$conf.int[1:2])
  }
  `colnames<-`(t(mapply(f, xl = x, tl = t)), paste0(prefix, c('rate', 'lower', 'upper')))
}

# life table
life.table <- function(mx, cohort = 100000) { # mx is the mortality rate
  n <- length(mx) + 1
  qx <- 2 * mx / (2 + mx) # probability that a person of age x will die within one year
  qx <- c(qx, 1) # forced method - mortality rate max age + 1 is 100%
  lx <- c(1, cumprod(1 - qx)) * cohort # number surviving
  dx <- -c(diff(lx), lx[n] * qx[n]) # number of deaths
  t <- (lx + c(lx[-1], 0)) / 2
  Tx <- rev(cumsum(rev(t))) # person-years after age x
  ex <- Tx / lx # life expectancy
  data.frame(lx = lx, dx = dx, t = t, Tx = Tx, ex = ex)[1:n,]
}

# add transparency to colours
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

# -----------------
# load example data
# -----------------

d <- read.csv(url('https://raw.githubusercontent.com/danlewer/autism-life-expectancy/main/simulated-data/A_noIDsim.csv'))
setDT(d)
setnames(d, old = c('deaths', 'pys'), new = c('all_cause', 'follow_up'))

# generate random cause-specific data
set.seed(3)
d[, external := floor(runif(.N, 0.25, 0.75) * all_cause)]
d[, external_suicide := rbinom(.N, external, prob = 0.5)]
d[, external_accident := external - external_suicide]
d[, other := all_cause - external]

# -----------------------------------------------------
# fit poisson models and visualise association with age
# -----------------------------------------------------

# make summary data by age groups
# -------------------------------

# add age group to data

age_lims <- 0:18 * 5
d[, age_group := findInterval(age, 0:18 * 5)]
d[, age_group := factor(age_group, seq_along(age_lims), age_lims)]

# summarise data by age groups and add labels

by_age_group <- d[, lapply(.SD, sum), c('exposure', 'age_group'), .SDcols = c('all_cause', 'external_suicide', 'external_accident', 'other', 'follow_up')]
pr <- CJ(age_group = seq(15, 90, 5), exposure = 0:1)
pr[, age := age_group + 2]
pr[, age_group := factor(age_group, age_lims, age_lims)]
labs <- data.table(labs = paste0(seq(15, 90, 5), '-', seq(19, 94, 5)), age_group = seq(15, 90, 5))
labs[, age_group := factor(age_group, age_lims, age_lims)]
labs$labs[1] <- '18-19'
pr <- labs[pr, on = 'age_group']
pr$age[pr$age == 17] <- 18.5
by_age_group <- merge(pr, by_age_group, all = T, by = c('age_group', 'exposure'))

# function to calculate observed rates, and rates estimated from poisson model

group <- function (outcome = 'all_cause', test_poly = 2) {
  actual_rates <- vpt(x = by_age_group[, get(outcome)], t = by_age_group$follow_up) * 1e5
  fs <- lapply(test_poly, function (x) as.formula(paste0(outcome, '~poly(age,', x, ')*exposure+offset(log(follow_up))')))
  ms <- lapply(fs, function (x) glm(x, data = d, family = 'poisson'))
  best_model <- which.min(sapply(ms, AIC))
  m <- ms[[best_model]]
  nd <- by_age_group[, .(exposure = exposure, age = age, follow_up = 1)]
  pred <- predict(m, newdata = nd, se.fit = T)
  f <- m$family$linkinv
  pred_rates <- data.table(pred_rate = f(pred$fit) * 1e5,
                           pred_lower = f(pred$fit - qnorm(0.975) * pred$se.fit) * 1e5,
                           pred_upper = f(pred$fit + qnorm(0.975) * pred$se.fit) * 1e5)
  cbind(by_age_group[, c('age_group', 'exposure', 'labs', 'age')],
        deaths = by_age_group[, get(outcome)],
        actual_rates, pred_rates)
}

# function to plot observed and modelled rates

plot_group <- function (outcome = 'all_cause', off = 0.5, ...) {
  pd <- group(outcome, ...)
  pd$x <- pd$age + ifelse(pd$exposure == 0, off, -off)
  pd$col <- fifelse(pd$exposure == 0, 'blue', 'red')
  pd$col2 <- add.alpha(pd$col, alpha = 0.2)
  pd$rate[pd$deaths == 0] <- NA
  pd$upper[pd$deaths == 0] <- NA
  pd$lower[pd$deaths == 0] <- NA
  pd$upper <- pmin(pd$upper, 1e5)
  pd$pred_upper <- pmin(pd$pred_upper, 1e5)
  pd$lower <- pmax(pd$lower, 1)
  pd$pred_lower <- pmax(pd$pred_lower, 1)
  plot(1, type = 'n', xlim = c(15, 95), ylim = c(0, log(100000)), axes = F, xlab = NA, ylab = NA)
  axis(1, pr$age, pr$labs, pos = 0, las = 2)
  rect(15, 0, 95, log(100000))
  axis(2, log(10^(0:5)), 10^(0:5), las = 2, pos = 15)
  lapply(split(pd, f = x$exposure), function (z) {
    with (z, {
      polygon(x = c(x, rev(x)), y = log(c(pred_lower, rev(pred_upper))), col = col2, border = NA)
      points(x, log(rate), pch = 19, col = col)
      arrows(x, log(lower), x, log(upper), angle = 90, code = 3, length = 0.03, col = col)
      lines(x, log(pred_rate), col = 'white', lwd = 3)
      lines(x, log(pred_rate), col = col)
    })
  })
}

# plots

par(mfrow = c(2, 2))
plot_group('all_cause', test_poly = 2); title(main = 'All cause')
plot_group('external_suicide', test_poly = 2); title(main = 'External (suicides)')
plot_group('external_accident', test_poly = 2); title(main = 'External (accidents)')
plot_group('other', test_poly = 2); title(main = 'Other')

# check cause-specific modelled rates add up

causes <- c('all_cause', 'other', 'external_suicide', 'external_accident')
crude_rates <- sapply(causes, function (x) group(x, test_poly = 2)$rate)
crude_rates <- as.data.frame.matrix(crude_rates)
crude_rates$sum_cause_specific <- rowSums(crude_rates[, c('other', 'external_suicide', 'external_accident')])
modelled_rates <- sapply(causes, function (x) group(x, test_poly = 2)$pred_rate)
modelled_rates <- as.data.frame.matrix(modelled_rates)
modelled_rates$sum_cause_specific <- rowSums(modelled_rates[, c('other', 'external_suicide', 'external_accident')])

# ------------------------------------
# estimate life expectancy and deficit 
# ------------------------------------

# add age squared (this can also be done using `poly(age, 2`), but manual coding makes the model easier to use
d[, age2 := age ^ 2]

# fit poisson models
m0 <- glm(all_cause ~ age + age2 + offset(log(follow_up)), data = d[exposure == 0], family = 'poisson')
m1 <- glm(all_cause ~ age + age2 + offset(log(follow_up)), data = d[exposure == 1], family = 'poisson')
