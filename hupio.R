
# -----------------------
# libraries and functions
# -----------------------

library(data.table)
library(MASS) # for mvrnorm - simulate values from multivariate normal distribution

# vectorized poisson confidence intervals
vpt <- function(x, t, form = F, digs = 2, ...) {
  f <- function(xl, tl) {
    if (is.na(xl) | is.na(tl)) return(c(0, 0, 0))
    y <- poisson.test(xl, tl, ...)
    c(xl/tl, y$conf.int[1:2])
  }
  res <- mapply(f, xl = x, tl = t)
  return(if (form) {
    res <- format(round(res, digs), digits = digs, nsmall = digs)
    res <- apply(res, 2, function(x) paste0(x[1], '(', x[2], '-', x[3], ')'))
    res <- gsub(' ', '', res)
    gsub('\\(', ' (', res)
  } else {
    `colnames<-`(t(res), c('rate', 'lower', 'upper'))
  })
}

# life tables
life.table <- function(mx, cohort = 100000) {
  n <- length(mx) + 1
  qx <- 2 * mx / (2 + mx)
  qx <- c(qx, 1) # forced method - mortality rate max age + 1 is 100%
  lx <- c(1, cumprod(1 - qx)) * cohort
  dx <- -c(diff(lx), lx[n] * qx[n])
  t <- (lx + c(lx[-1], 0)) / 2
  Tx <- rev(cumsum(rev(t)))
  ex <- Tx / lx
  data.frame(lx = lx, dx = dx, t = t, Tx = Tx, ex = ex)[1:n,]
}

# add transparency to colours
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

# --------------------------------------------------
# estimate SMR with monte-carlo confidence intervals
# --------------------------------------------------

smrV <- function(d, B = 10) {
  smrd <- d[, .(deaths = sum(deaths), pys = sum(pys)), c('exposure', 'sex', 'age_grp')]
  full_table <- CJ(exposure = 0:1, sex = 1:2, age_grp = as.factor(seq(15, 90, 5)))
  smrd <- merge(full_table, smrd, by = c('exposure', 'sex', 'age_grp'), all = T)
  smrd[is.na(smrd)] <- 0
  Bs <- rbindlist(rep(list(smrd), B))
  Bs$b <- rep(seq_len(B), each = nrow(smrd))
  Bs$r <- rpois(nrow(Bs), Bs$deaths)
  BsRef <- Bs[exposure == 0]
  BsRef$rt <- BsRef$r / BsRef$pys
  Bs <- Bs[exposure == 1][BsRef, on = c('sex', 'age_grp', 'b')]
  Bs <- Bs[pys != 0 & !is.nan(rt)]
  Bs$e <- Bs$pys * Bs$rt
  Bs <- Bs[, .(o = sum(r), e = sum(e)), b]
  quantile(Bs$o / Bs$e, probs = c(0.5, 0.025, 0.975))
}

# -------------------------------------------
# use model to estimate annual mortality rate
# -------------------------------------------

age_lims <- 0:18 * 5
pr <- CJ(age_grp = seq(15, 90, 5), exposure = 0:1)
pr[, age := age_grp + 2]
pr[, age_grp := factor(age_grp, age_lims, age_lims)]
labs <- data.table(labs = paste0(seq(15, 90, 5), '-', seq(19, 94, 5)), age_grp = seq(15, 90, 5))
labs[, age_grp := factor(age_grp, age_lims, age_lims)]
labs$labs[1] <- '18-19'
pr <- labs[pr, on = 'age_grp']
pr$age[pr$age == 17] <- 18.5

# function

mblt <- function (data, SEX = 1, B = 10, polyn = 3) { # model-based life table
  
  d2 <- data[sex == SEX]
  d2$age1 <- d2$age
  d2$age2 <- d2$age ^ 2
  d2$age3 <- d2$age ^ 3
  d2$age4 <- d2$age ^ 4
  d2$age5 <- d2$age ^ 5
  f <- as.formula(paste0('deaths ~', paste0('age', 1:polyn, collapse = '+'), '+ offset(log(pys))'))
  m <- glm(deaths ~ poly(age, polyn)*exposure + offset(log(pys)), data = d2, family = 'poisson')
  m0 <- glm(f, data = d2[exposure == 0], family = 'poisson')
  m1 <- glm(f, data = d2[exposure == 1], family = 'poisson')
  
  # summary table with age groups
  d2g <- d2[, lapply(.SD, sum), c('exposure', 'age_grp'), .SDcols = c('deaths', 'pys')]
  d2g <- cbind(d2g, with(d2g, vpt(deaths, pys)) * 100000)
  d2g[, upper := pmin(upper, 100000)]
  d2g[, lower := pmax(lower, 1)]
  d2g$rate[d2g$deaths == 0] <- NA
  d2g$lower[d2g$deaths == 0] <- NA
  d2g$upper[d2g$deaths == 0] <- NA
  d2g <- merge(pr, d2g, all = T, by = c('age_grp', 'exposure'))
  d2g[, pys := 1]
  pred <- predict(m, newdata = d2g, se.fit = T)
  f <- m$family$linkinv
  d2g[, pred_rate := f(pred$fit) * 100000]
  d2g[, pred_rate_lower := f(pred$fit - qnorm(0.975) * pred$se.fit) * 100000]
  d2g[, pred_rate_upper := f(pred$fit + qnorm(0.975) * pred$se.fit) * 100000]
  
  # predicted rates
  nd <- CJ(age = 18:100, pys = 1, exposure = 0)
  nd[, age1 := age]
  nd[, age2 := age ^ 2]
  nd[, age3 := age ^ 3]
  nd[, age4 := age ^ 4]
  nd[, age5 := age ^ 5]
  nd[, rt0 := predict(m0, newdata = nd, type = 'response')]
  nd[, rt1 := predict(m1, newdata = nd, type = 'response')]
  nd$rt0[nrow(nd)] <- 1
  nd$rt1[nrow(nd)] <- 1
  
  # life tables
  lt0 <- life.table(mx = nd$rt0)
  lt1 <- life.table(mx = nd$rt1)
  
  # simulate confidence intervals for life expectancy
  Bf <- function (model) {
    Bs <- mvrnorm(B, mu = coef(model), Sigma = vcov(model))
    apply(Bs, 1, function (x) {
      age_matrix <- outer(18:100, 1:10, FUN = '^')[, 1:polyn]
      mx <- exp(rowSums(cbind(x[1], t(x[-1] * t(age_matrix)))))
      lt <- life.table(mx)
      lt$ex[1]
    })
  }
  BsEx <- rbind(Bf(m0), Bf(m1))
  
  # outputs
  list(ex = c(lt0$ex[1], lt1$ex[1]), nd = nd, BsEx = BsEx, d2g = d2g)
  
}

# extract results

sr <- function (x) {
  list (
    life_expectancy = `names<-`(c(x$ex + 18, -diff(x$ex)), c('unexposed', 'exposed', 'diff')),
    conf_ints = `colnames<-`(cbind(apply(x$BsEx, 1, quantile, probs = c(0.5, 0.025, 0.975)) + 18,
                                   quantile(x$BsEx[1,] - x$BsEx[2,], probs = c(0.5, 0.025, 0.975))), c('unexposed', 'exposed', 'diff'))
  )}

sr2 <- function (x, dig = 2) {
  x <- rbind(x$life_expectancy, x$conf_ints[2:3,])
  x <- format(round(x, dig), nsmall = dig, digits = dig)
  r <- apply(x, 2, function(y) paste0(y[1], '(', y[2], ', ', y[3], ')'))
  r <- gsub('\\(', ' (', gsub(' ', '', r))
  gsub(',', ', ', r)
}

# -----------
# do analysis
# -----------

do_analysis <- function (data, polyn = 2) {
  d <- data[age <= 94]
  # add age groups
  age_lims <- 0:18 * 5
  d$age_grp <- findInterval(d$age, age_lims)
  d$age_grp <- factor(d$age_grp, seq_along(age_lims), age_lims)
  # SMR
  smr <- smrV(d, B = 10000)
  # life tables
  male <- mblt(SEX = 1, B = 10000, data = d, polyn = polyn)
  female <- mblt(SEX = 2, B = 10000, data = d, polyn = polyn)
  # results
  list(smr = smrV(d, B = 10000),
       male = sr2(sr(male)),
       female = sr2(sr(female)),
       male_table = male$d2g,
       female_table = female$d2g)
}

d <- fread('https://raw.githubusercontent.com/danlewer/hupio/main/mortality/single-year-of-age-rates/all_cause.csv', col.names = c('exposure', 'age', 'sex', 'deaths', 'pys'))
d[, exposure := fifelse(exposure == T, 1, 0)]
d[, sex := fifelse(sex == 'male', 1, 2)]
d <- d[order(sex, age, exposure)]

set.seed(4668)

res <- do_analysis(d, polyn = 2)

# -------------------
# plot modelled rates
# -------------------

make_plot_data <- function (data, off = 0.5) {
  data$x <- data$age + fifelse(data$exposure == 0, off, -off)
  data$col <- fifelse(data$exposure == 0, 'blue', 'red')
  data$col2 <- add.alpha(data$col, alpha = 0.2)
  data
}

pd_male <- make_plot_data(res$male_table)
pd_female <- make_plot_data(res$female_table)

plot_lines <- function (x) {
  y <- split(x, f = x$exposure)
  lapply(y, function (z) {
    with (z, {
      polygon(x = c(x, rev(x)), y = log(c(pred_rate_lower, rev(pred_rate_upper))), col = col2, border = NA)
      points(x, log(rate), pch = 19, col = col)
      arrows(x, log(lower), x, log(upper), angle = 90, code = 3, length = 0.03, col = col)
      lines(x, log(pred_rate), col = 'white', lwd = 3)
      lines(x, log(pred_rate), col = col)
    })
  })
}

#png('age_group_modelled_rates.png', height = 9, width = 9, units = 'in', res = 300)

par(mfrow = c(1, 2), mar = c(4, 0, 3, 0), oma = c(0, 5, 0, 0))

plot(1, type = 'n', xlim = c(15, 95), ylim = c(0, log(100000)), axes = F, xlab = NA, ylab = NA)
plot_lines(pd_male)
axis(1, pr$age, pr$labs, pos = 0, las = 2)
rect(15, 0, 95, log(100000))
axis(2, log(10^(0:5)), 10^(0:5), las = 2, pos = 15)
title(main = 'Male', line = 0.5)

plot(1, type = 'n', xlim = c(15, 95), ylim = c(0, log(100000)), axes = F, xlab = NA, ylab = NA)
plot_lines(pd_female)
axis(1, pr$age, pr$labs, pos = 0, las = 2)
rect(15, 0, 95, log(100000))
title(main = 'Female', line = 0.5)

mtext('Mortality rate per 100,000 person-years', side = 2, outer = T, line = 3.5)

#dev.off()

