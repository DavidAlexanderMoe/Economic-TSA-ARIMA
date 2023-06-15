TSA ARIMA
================
David Alexander Moe
2023-06-15

``` r
rm(list = ls())

# libraries and functions
library(forecast)
```

    ## Registered S3 method overwritten by 'quantmod':
    ##   method            from
    ##   as.zoo.data.frame zoo

``` r
library(tseries)
library(tsoutliers)
library(urca)
library(FinTS)
```

    ## Caricamento del pacchetto richiesto: zoo

    ## 
    ## Caricamento pacchetto: 'zoo'

    ## I seguenti oggetti sono mascherati da 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

    ## 
    ## Caricamento pacchetto: 'FinTS'

    ## Il seguente oggetto è mascherato da 'package:forecast':
    ## 
    ##     Acf

``` r
source("G:/Il mio Drive/TRIENNALE/TS Analysis/Economic-TSA-ARIMA/Functions/TSA-Useful-Functions.R")
source("G:/Il mio Drive/TRIENNALE/TS Analysis/Economic-TSA-ARIMA/Functions/CalendarEffects-Functions.R")
source("G:/Il mio Drive/TRIENNALE/TS Analysis/Economic-TSA-ARIMA/Functions/TSA-Predict-Functions.R")
```

``` r
filedati <- "G:\\Il mio Drive\\TRIENNALE\\TS Analysis\\Economic-TSA-ARIMA\\Production in industry - monthly data (Italy).csv"

#### Read data
dati <- read.table(file = filedati, header = TRUE, sep = ",", quote = "", 
                   na.strings = ".", check.names = FALSE, comment.char = "")
name <- "Produzione industriale in Italia"

dati$TIME_PERIOD <- paste(dati$TIME_PERIOD,"-15",sep="")
dati$TIME_PERIOD <- gsub("-","",dati$TIME_PERIOD)
dati$TIME_PERIOD <- as.Date(x=dati$TIME_PERIOD, format= "%Y %m %d")

start <- as.numeric( c( format(dati$TIME_PERIOD[1], "%Y"), format(dati$TIME_PERIOD[1], "%m") ) )

y <- ts(data=dati$OBS_VALUE, start= start, frequency=12)
```

``` r
# some external variables

#### calendar effects
cal <- .calendarEffects(time = dati$TIME_PERIOD, country = "it")
```

    ## Caricamento del pacchetto richiesto: timeDate

    ## Caricamento del pacchetto richiesto: car

    ## Caricamento del pacchetto richiesto: carData

``` r
#cal <- cal[, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat","Sun","sh", "lh", "eh"), drop = FALSE]
cal <- cal[, c("wd", "lh", "eh"), drop = FALSE]

cal <- as.matrix(cal)

#### Drift if needed
drift <- cbind(drift = 1 : NROW(y))
```
