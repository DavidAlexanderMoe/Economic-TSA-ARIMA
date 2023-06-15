################################################################################
##
## File:    TSA-Useful-Functions.R
## 
## Purpose: Some useful functions.
## 
################################################################################

.circle <- function(win = 1.3)
{
  ##############################################################################
  ## Arguments
  ##  win:     (numeric[1]) Window size.
  ## Value
  ##  Plot a unit radius circle.  
  ##############################################################################
  
  #### Initialize plot
  plot(x = c(-win, win), c(-win, win), type = "n", xlab = "", ylab = "")
  
  #### Prepare "circle data"
  radius <- 1
  theta <- seq(from = 0, to = 2 * pi, length = 200)
  
  #### Draw the circle
  lines(x = radius * cos(theta), y = radius * sin(theta))
  
  #### Add axes
  abline(h = 0)
  abline(v = 0)
  
  #### Answer
  invisible(NULL)
}
# ------------------------------------------------------------------------------


.acf.theoretical <- function(model, lag.max, 
                             ylim = c(-1, 1), pacf = FALSE, plot = TRUE)
{
  #### Compute
  S <- ifelse( NROW(model$S) > 0, model$S[1], 12) 
  model <- .sarma2larma(ar = model$ar, ma = model$ma, 
                        sar = model$sar, sma = model$sma, S = S)
  acf <- ARMAacf(ar = model$ar, ma = model$ma, lag.max = lag.max, pacf = pacf)
  #### Plot
  if (plot)
  {
    ind <- (1 + !pacf) : NROW(acf)
    x1 <- acf[ind]  
    main <- paste0("Theoretical ", ifelse(pacf, "PACF", "ACF"))
    ylab <- paste0(ifelse(pacf, "Partial ", ""), "ACF")
    plot(x1, type = "h", main = main, xlab = "Lag", ylab = ylab, 
         ylim = ylim)
    abline(a = 0, b = 0)
  }
  #### Answer
  acf
}
# ------------------------------------------------------------------------------


################################################################################
## Functions to simulate
################################################################################

#### Merge AR(p) and SAR(P); the same for MA(q) and SMA(Q) 
.sarma2larma <- function(ar = NULL, ma = NULL, sar = NULL, sma = NULL, S = 12)
{
  #### Adjust
  if ( NROW(ar) > 0 ) { ar <- -ar }
  if ( NROW(sar) > 0 ) { sar <- -sar }
  #### model
  list(
    ar = -.long(p = ar, ps = sar, s = S),
    ma =  .long(p = ma, ps = sma, s = S) )
}
# ------------------------------------------------------------------------------


#### Merge short and seasonal components 
.long <- function(p, ps, s)
{
  #### Settings
  np  <- NROW( p )
  nps <- NROW( ps )
  
  #### 
  cp  <- if ( np > 0 ) { c(1, p) } else { 1 }
  cps <- if ( nps > 0 )
  { 
    ind <- seq(from = s, by = s, length.out = nps)
    x1 <- numeric(s * nps)
    x1[ind] <- ps
    c(1, x1)
  }
  else
  {
    1
  }
  
  #### Answer
  convolve(cp, rev(cps), type = "open")[-1]
}
# ------------------------------------------------------------------------------


#### Simulate ARIMA(p,d,q)x(P,D,Q)S
.arima.sim <- function(
    model = list(d = 0, D = 0, S = 12, 
                 ar = NULL, ma = NULL, sar = NULL, sma = NULL), 
    innov, nburn)
{
  #### Orders
  p  <- NROW(model$ar)  
  d <- ifelse( NROW(model$d) > 0, model$d[1], 0) 
  q <- NROW(model$ma)
  ps <- NROW(model$sar)
  ds <- ifelse( NROW(model$D) > 0, model$D[1], 0) 
  qs <- NROW(model$sma)
  S <- ifelse( NROW(model$S) > 0, model$S[1], 12) 
  
  #### Long form of parameters
  model <- .sarma2larma(ar = model$ar, ma = model$ma, 
                        sar = model$sar, sma = model$sma, S = S)
  
  #### Set innov as a ts() object
  innov <- ts(data = innov, start = 1, frequency = S)
  
  #### Set Arima object
  require(forecast)
  obj <- Arima(y = innov, order = c(p, d, q), seasonal = c(ps, ds, qs),
               include.mean = TRUE)
  # print(obj)
  obj$model$phi <- model$ar
  obj$model$theta <- model$ma
  
  #### Simulate
  simulate(object = obj, innov = innov, future = FALSE)[(nburn + 1) : NROW(innov)]
}
# ------------------------------------------------------------------------------


################################################################################
## Functions to plot AR, MA roots
################################################################################

.arma.roots <- function(fit)
{
  #### Coefficients
  coef <- fit$coef
  ar <- coef[substr(names(coef), 1, 2) == "ar"]
  ma <- coef[substr(names(coef), 1, 2) == "ma"]
  sar <- coef[substr(names(coef), 1, 3) == "sar"]
  sma <- coef[substr(names(coef), 1, 3) == "sma"]
  S <- fit$arma[5]
  
  #### Long form of parameters
  long <- .sarma2larma(ar = ar, ma = ma, sar = sar, sma = sma, S = S)
  #### Roots
  ar.r <- polyroot(c(1, -long$ar))
  ma.r <- polyroot(c(1, long$ma))
  
  #### Answer
  list(coef = list(ar = ar, ma = ma, sar = sar, sma = sma), 
       coef.long = long, 
       root = list(ar = ar.r, ma = ma.r))
}
# ------------------------------------------------------------------------------


################################################################################
## Extact Arima settings
################################################################################

.constant.type <- function(fit)
{
  if (any( names(fit$coef) == "drift")) {"drift"} 
  else if (any( names(fit$coef) == "mean") ) {"mean"} 
  else if (any( names(fit$coef) == "intercept") ) {"intercept"} 
  else {NULL}
}
# -----------------------------------------------------------------------------

.Arima.settings <- function(fit)
{
  constant.type <- .constant.type(fit = fit)
  list(
    constant.type = constant.type,
    include.constant = NROW(constant.type) > 0,
    include.drift = NROW(constant.type) > 0 && constant.type == "drift", 
    order = fit$arma[c(1,6,2)],
    seasonal = list(order = fit$arma[c(3,7,4)]) )
}
# -----------------------------------------------------------------------------


################################################################################
## Transformation diagnostic
################################################################################

.trsf.test <- function(fit, msg = "")
{
  #### Load
  require(sandwich)
  #### Extract
  mod <- fit
  fit <- fitted(mod)
  res <- residuals(mod)
  #### Fit
  lm1 <- lm( log(res^2) ~ fit )
  
  #### Coef
  coef   <- coef(lm1)
  #### vcov
  vcov    <- vcovHC(x = lm1, type = "const")
  vcovHC  <- vcovHC(x = lm1)
  vcovHAC <- vcovHAC(x = lm1)
  #### Test
  coef  <- coef["fit"]
  se    <- sqrt(vcov["fit", "fit"])   
  seHC  <- sqrt(vcovHC["fit", "fit"]) 
  seHAC <- sqrt(vcovHAC["fit", "fit"])
  tab <- data.frame(
    estimate     = rep.int(coef, 3),  
    se           = rep.int(se, 3),
    H0           = c(          "fit = 0",          "fit = 1",        "fit = 2"   ),  
    tstat        = c(    (coef - 0) / se,    (coef - 1) / se,   (coef - 2) / se  ), 
    "se(HC)"     = rep.int(seHC, 3),
    "tstat(HC)"  = c(  (coef - 0) / seHC,  (coef - 1) / seHC,  (coef - 2) / seHC ), 
    "se(HAC)"    = rep.int(seHAC, 3),
    "tstat(HAC)" = c( (coef - 0) / seHAC, (coef - 1) / seHAC, (coef - 2) / seHAC ), 
    check.names = FALSE )   
  
  #### Print
  if ( msg != "" )
  {
    cat(msg)
    print( tab )
  }
  #### Answer
  tab 
}
# ------------------------------------------------------------------------------


################################################################################
## UR test with RW + drift under H0
################################################################################

.ur.drift <- function(y, 
                      lags = 1, selectlags = c("Fixed", "AIC", "BIC")) 
{
  selectlags<-match.arg(selectlags)
  if (ncol(as.matrix(y)) > 1) 
    stop("\ny is not a vector or univariate time series.\n")
  if (any(is.na(y))) 
    stop("\nNAs in y.\n")
  y <- as.vector(y)
  lag <- as.integer(lags)
  if (lag < 0) 
    stop("\nLags must be set to an non negative integer value.\n")
  CALL <- match.call()
  DNAME <- deparse(substitute(y))
  x.name <- deparse(substitute(y))
  lags <- lags + 1
  z <- diff(y)
  n <- length(z)
  x <- embed(z, lags)
  z.diff <- x[, 1]
  z.lag.1 <- y[lags:n]
  tt <- lags:n
  if (lags > 1) {
    if(selectlags!="Fixed"){
      critRes<-rep(NA, lags)
      for(i in 2:(lags)){
        z.diff.lag = x[, 2:i]
        result <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)  
        critRes[i]<-AIC(result, k = switch(selectlags, "AIC" = 2, "BIC" = log(length(z.diff))))
      }
      lags<-which.min(critRes)
    }
    z.diff.lag = x[, 2:lags]
    result <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
    tau <- coef(summary(result))[2, 3]
    teststat <- as.matrix(tau)
    colnames(teststat) <- c('t')
  }
  else {
    result <- lm(z.diff ~ z.lag.1 + 1)
    tau <- coef(summary(result))[2, 3]
    teststat <- as.matrix(tau)
    colnames(teststat) <- c('t')
  }
  rownames(teststat) <- 'statistic'
  testreg <- summary(result)
  res <- residuals(testreg)
  cval.t2 <- qt(df = result$df.residual, p = c(0.01, 0.05, 0.10))
  cvals <- rbind(t2 = cval.t2)
  testnames <- 't'
  colnames(cvals) <- c("1pct", "5pct", "10pct")
  rownames(cvals) <- testnames
  
  new("ur.df", y = y, cval=cvals, lags=lag, teststat = teststat, testreg=testreg, res=res, test.name="UR(drift) Test")
}