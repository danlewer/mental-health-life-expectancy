library(stringi)

# https://academic.oup.com/ije/article/41/6/1585/741862

# external - V01:V99, W00:W99, X00:X99, Y00:Y99
external <- c(paste0('V', stri_pad(1:99, width = 2, pad = '0')),
              paste0('W', stri_pad(0:99, width = 2, pad = '0')),
              paste0('X', stri_pad(0:99, width = 2, pad = '0')),
              paste0('Y', stri_pad(0:99, width = 2, pad = '0')))

# suicide - X60:X84 and Y10:Y34
suicide <- c(paste0('X', stri_pad(60:84, width = 2, pad = '0')),
             paste0('Y', stri_pad(10:34, width = 2, pad = '0')))
