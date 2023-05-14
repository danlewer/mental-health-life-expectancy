library(stringi)

# ---------------
# causes of death
# ---------------

# external - V01:V99, W00:W99, X00:X99, Y00:Y99
external <- c(paste0('V', stri_pad(1:99, width = 2, pad = '0')),
              paste0('W', stri_pad(0:99, width = 2, pad = '0')),
              paste0('X', stri_pad(0:99, width = 2, pad = '0')),
              paste0('Y', stri_pad(0:99, width = 2, pad = '0')))

# suicide - X60:X84 and Y10:Y34
suicide <- c(paste0('X', stri_pad(60:84, width = 2, pad = '0')),
             paste0('Y', stri_pad(10:34, width = 2, pad = '0')))

n <- 2000
poss_entry_dates <- seq(as.Date('2000-01-01'), as.Date('2015-01-01'), by = 'day')

set.seed(4)
d <- data.frame(exposure = rep(c(T, F), each = n/2),
           age_at_entry = sample(18:45, n, T) + sample(0:365, n, T)/365,
           sex = sample(c('m', 'f'), n, replace = T),
           entry_date = sample(poss_entry_dates, n, T))
d$exit_date <- d$entry_date + runif(n, 0, 1000)
d$died <- rbinom(n, 1, 0.05)
d$icd10 <- sample(c(external, rep(paste0('I', stri_pad(0:99, width = 2, pad = '0')), 5)), n, T)
d$icd10[d$died == 0] <- NA_character_
