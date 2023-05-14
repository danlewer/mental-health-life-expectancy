library(data.table)

# -----------------
# load example data
# -----------------

d <- fread('https://raw.githubusercontent.com/danlewer/life-expectancy-mh/main/lexis-expansion/example_individual_data.csv')
# format dates
date_cols <- c('entry_date', 'exit_date')
d[, (date_cols) := lapply(.SD, as.Date), .SDcols = date_cols]
# if your data has birth date, no need to calculate it
d[, date_of_birth := entry_date - age_at_entry * 365.25]
d[, id := .I]

# -------------------------
# categorise cause of death
# -------------------------

# if you do not have cause of death, just run
# d[, cause := NA]

# external - V01:V99, W00:W99, X00:X99, Y00:Y99
external <- c(paste0('V', stri_pad(1:99, width = 2, pad = '0')),
              paste0('W', stri_pad(0:99, width = 2, pad = '0')),
              paste0('X', stri_pad(0:99, width = 2, pad = '0')),
              paste0('Y', stri_pad(0:99, width = 2, pad = '0')))

# suicide - X60:X84 and Y10:Y34
suicide <- c(paste0('X', stri_pad(60:84, width = 2, pad = '0')),
             paste0('Y', stri_pad(10:34, width = 2, pad = '0')))

d[, cause := 'other']
d$cause[d$icd10 %in% external] <- 'external_accident'
d$cause[d$icd10 %in% suicide] <- 'external_suicide'
d$cause[d$died == 0] <- NA_character_
d[, cause := factor(cause, levels = c('external_suicide', 'external_accident', 'other'))]

# --------------------------------------------
# create aggregated data with time-varying age
# --------------------------------------------

# make dates into integers

date_cols <- c('entry_date', 'exit_date', 'date_of_birth')
d[, (date_cols) := lapply(.SD, as.integer), .SDcols = date_cols]

# function for creating lexis-expanded data

expand <- function (data, number_splits = 5, cause_specific = T)  {
  nr <- seq_len(nrow(data))
  ind <- findInterval(nr, quantile(nr, 0:(number_splits-1)/number_splits))
  x <- split(data, ind)
  l <- lapply(seq_along(x), function (z) {
    days <- x[[z]][, .(day = seq.int(from = entry_date, to = exit_date)), c('id', 'date_of_birth', 'sex', 'exposure', 'died', 'exit_date', 'cause')]
    days[, age := (day - date_of_birth) / 365.25]
    days[, age := floor(age)]
    days$died[days$day != days$exit_date] <- NA_integer_
    days$cause[days$day != days$exit_date] <- NA_character_
    follow_up <- days[, .(follow_up = .N), c('exposure', 'sex', 'age')]
    all_cause_deaths <- days[, .(all_cause = sum(died, na.rm = T)), c('exposure', 'sex', 'age')]
    output <- all_cause_deaths[follow_up, on = c('exposure', 'sex', 'age')]
    if (cause_specific == T) {
      cause_specific_deaths <- dcast(days, exposure + sex + age ~ cause, value.var = 'id', fun.aggregate = length, drop = F)
      output <- cause_specific_deaths[output, on = c('exposure', 'sex', 'age')]
    }
    output[is.na(output)] <- 0L
    print (paste0(z, '/', number_splits))
    return (output)
  })
  agg <- rbindlist(l)
  agg[, lapply(.SD, sum), by = c('exposure', 'sex', 'age')]
}

# use function to create aggregated lexis data
# a higher value of 'number_splits' will use less RAM but may take longer

lexis_data <- expand(data = d, number_splits = 5, cause_specific = T)
lexis_data <- lexis_data[order(exposure, sex, age)]
fwrite(lexis_data, 'lexis_results.csv')
