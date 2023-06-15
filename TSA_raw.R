################################################################################
##
## File:  TS-Economica.csv
## 
## Scopo: Produzione industriale mensile in Italia
##      
################################################################################

################################################################################
## Pulisco
################################################################################

rm(list = ls())

################################################################################
## Librerie e funzioni
################################################################################

library(forecast)
library(tseries)
library(tsoutliers)
library(urca)
library(FinTS)

source("G:/Il mio Drive/TRIENNALE/TS Analysis/TSA-Useful-Functions.R")
source("G:/Il mio Drive/TRIENNALE/TS Analysis/CalendarEffects-Student-Functions.R")
source("G:/Il mio Drive/TRIENNALE/TS Analysis/TSA-Predict-Student-Functions.R")

################################################################################
## Leggo i dati
################################################################################

filedati <- "G:\\Il mio Drive\\TRIENNALE\\TS Analysis\\Production in industry - monthly data (Italy).csv"

#### Read data
dati <- read.table(file = filedati, header = TRUE, sep = ",", quote = "", 
                   na.strings = ".", check.names = FALSE, comment.char = "")
name <- "Produzione industriale in Italia"

### ts() e formattazione date
dati$TIME_PERIOD <- paste(dati$TIME_PERIOD,"-15",sep="")
dati$TIME_PERIOD <- gsub("-","",dati$TIME_PERIOD)
dati$TIME_PERIOD <- as.Date(x=dati$TIME_PERIOD, format= "%Y %m %d")

start <- as.numeric( c( format(dati$TIME_PERIOD[1], "%Y"), format(dati$TIME_PERIOD[1], "%m") ) )

y <- ts(data=dati$OBS_VALUE, start= start, frequency=12)


################################################################################
## Variabili esterne
################################################################################

#### Effetti di calendario
cal <- .calendarEffects(time = dati$TIME_PERIOD, country = "it")

#cal <- cal[, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat","Sun","sh", "lh", "eh"), drop = FALSE]
cal <- cal[, c("wd", "lh", "eh"), drop = FALSE]

# uso wd al posto dei gg della settimana perche vedo che piu o meno sono tutti e 5 significativi allo stesso modo
# wd ? centrato sullo 0 quindi ? normale venga basso il coeff
# confronto tramite BIC (prendo il piu basso) il modello con i gg della settimana e quello con wd ->prendo con wd

cal <- as.matrix(cal)

#### Drift in caso servisse
drift <- cbind(drift = 1 : NROW(y))


################################################################################
## Preliminary analysis
################################################################################

#### Ts plot, acf, pacf of the original series
par(mfrow = c(3,1))
x <- y
plot(x, type = "l", main = name, ylab = "N.I.",xlab="Tempo")
Acf(x = x, type = "correlation", na.action = na.pass, lag.max = 60, main = name)
Acf(x = x, type = "partial",     na.action = na.pass, lag.max = 60, main = name)
#lag in anni
#andamento stagionale visibile dall' acf (decadimenti a lag multipli di 12 (ho gli anni))
#e picco al fi a lag 12 (sarebbe 1 sull pacf) -> spia di non staz
#pi? di questo non posso dire

#adesso guardiamo i residui del modello ARIMA(0,0,0)x(0,1,0)[12]
#se il proc ? un ma puro sull acf va giu in maniera esponenziale con "l'onda"
#nella parziale se ? un ar puro ci dovrebbero essere n lag signif quanto ? l ordine del processo
#per capire l ordine dell ar guardo la parziale e per l ma invece la totale
#in questo caso alzo il p


#### UNIT ROOT TESTS:
#### DF/ADF tests with trend
#### The following rows show how to make DF test in practice.
##   DF is however not recommended because gives severely biased results 
##   in case the DGP is more complex than the equation estimated by DF.
##   This time series prove this (compare DF with ADF results)
df.1 <- ur.df(y = y, type = "trend", lags = 0, selectlags = "Fixed") 
df.2 <- ur.df(y = y, type = "drift", lags = 0, selectlags = "Fixed") 
df.3 <- ur.df(y = y, type = "none",  lags = 0, selectlags = "Fixed")


#### (DGP:   RW + drift; 
##    Model: AR(1) + trend (+ other possible stationary terms))

adf.1 <- ur.df(y = y, type = "trend", lags = 24, selectlags = "AIC")
print( adf.1@teststat )
print( adf.1@cval )
##   1) H0(1) -> tau3 A and H0(2) -> phi3 A => no trend, go to 2   
##   2) H0(3) -> phi2 A                     => no drift, go to 3   

#### (DGP:   RW; 
##    Model: AR(1) + drift (+ other possible stationary terms))

adf.2 <- ur.df(y = y, type = "drift", lags = 24, selectlags = "AIC")
print( adf.2@teststat )
print( adf.2@cval )

#accetto tutto -> radici unitarie (rw senza drift)
#non staz --> d+D>0

#quante radici ho e di chi ? la colpa?
#rifiuto 
#la unit root che ho sembrerebbe essere colpa della stagionalit? quindi
#faccio le differenze dodicesime e noto che non ho ur dalla parte non stagionale
dod <- diff(y,lag=12)
adf.3 <- ur.df(y = dod, type = "trend", lags = 24, selectlags = "AIC")
print( adf.3@teststat )
print( adf.3@cval )



################################################################################
## ARIMA modeling ARIMA(p,d,q)x(P,D,Q)_S
################################################################################

#### Modeling
# AIC = -2 * loglikelihood +      2 * #par
# BIC = -2 * loglikelihood + log(T) * #par

xreg <- NULL
fit <- Arima(y = y, 
             order = c(0, 0, 0), seasonal = list(order = c(0, 1, 0)),
             xreg = xreg, include.constant = FALSE)
print(summary(fit))
base <- fit

#ARIMA(0,0,0)(0,1,0)[12]  AIC=2624.3   AICc=2624.31   BIC=2628.21
#ARIMA(1,0,0)(0,1,0)[12]  AIC=2428.83  AICc=2428.87   BIC=2436.66
#ARIMA(2,0,0)(0,1,0)[12]  AIC=2410.39  AICc=2410.45   BIC=2422.13
#ARIMA(3,0,0)(0,1,0)[12]  AIC=2394.99  AICc=2395.09   BIC=2410.64
#ARIMA(3,0,0)(1,1,0)[12]  AIC=2371.9   AICc=2372.06   BIC=2391.46
#ARIMA(3,0,0)(2,1,0)[12]  AIC=2347.6   AICc=2347.83   BIC=2371.08
#ARIMA(3,0,0)(2,1,1)[12]  AIC=2291.87  AICc=2292.18   BIC=2319.27
#ARIMA(3,0,0)(2,1,2)[12]  AIC=2247.26  AICc=2247.65   BIC=2278.56 ? il modello meno peggio

#commentare i residui dell ultimo modello -> ormai dell' andamento non importa pi? di tanto, 
#ci interessa fondamentalmente che iniziano ad assomigliare ad un WN
#test LB: rifiuto quasi sempre tranne a lag 22 e 27 (provenienti da H-#parametri (npa1+lag1-npar1 mi rende il vettore che mi serve))
#sembrerebbe che meglio di cosi non posso fare

#ho provato anche a modellare con d=1 ma ho ottenuto risultati migliori con d=0 soprattutto nell' acf


#### ARIMA (no external regressors)
xreg <- NULL
fit <- Arima(y = y, 
             order = c(3, 0, 0), seasonal = list(order = c(2, 1, 2)),
             xreg = xreg, include.constant = FALSE)
print(summary(fit))
fit1 <- fit

#### ARIMA + external regressors (var di calendario)
xreg <- cal
fit <- Arima(y = y,
             order = c(3, 0, 0), seasonal = list(order = c(0, 1, 1)),
             xreg = xreg, include.constant = FALSE)
print(summary(fit))
fit2 <- fit

#### Root analysis
#parte ar(p) con A(L)=0 e ma(q) con B(L)=0
par(mfrow = c(1,2))

#fit1
root <- .arma.roots(fit = fit1)
.circle(win = 2)
points(root$root$ar, col = "red")
points(root$root$ma, col = "blue")

#fit2
root <- .arma.roots(fit = fit2)
.circle(win = 2)
points(root$root$ar, col = "red")
points(root$root$ma, col = "blue")

## Comment: 
### fit1 shows MA roots (stemming from the seasonal part) critically close to 1
### in module.
### parte MA invertibile (radici B(L)B*(L^s)=0) -> arima invertibile
### alcune radici ar molto vicine a quelle ma -> abbasso Q da 2 a 1 nel modello+cal e rifaccio
### e ottengo risultati decisamenti migliori avendo chiaramente radici pi? separate




################################################################################
## ARIMA modeling with anomalies
################################################################################

#####################################
## ARIMA
#####################################

#### Copy model
fit <- fit1
#### Extract settings
settings <- .Arima.settings(fit = fit)
#### Prepare xreg
xreg <- NULL

#### Fit
fit <- tso(y = y, xreg = xreg,
           types = c("AO", "LS", "TC"), delta = 0.7, cval = 5,
           maxit = 10, maxit.iloop = 100, maxit.oloop = 10,
           tsmethod = "arima",
           args.tsmethod = list( order = settings$order, seasonal = settings$seasonal) )
fit1.o <- fit

# delta=0.7 di default
# cval=n lo alziamo se troviamo troppe anomalie oppure se l algoritmo non riesce a convergere
# gli altri sono un po alzati rispetto al default 


#### Reporting
print(fit)
plot(fit)

#ind -> osservazione
#time -> anno/mese
#stat test -> DI SOLITO MA NON SEMPRE sono maggiori del cval messo in tso()
#i punti rossi sono le oss quando si manifestano le anomalie
#la riga blu ? la ts aggiustata depurata dalle anomalie
#la grigia ? la ts vera per fare un confronto
#da scrivere nel report anche questo

#### Extract outlier effects dal modello stimato sopra
oeff <- outliers.effects(mo = fit$outliers, n = NROW(y), pars = coef(fit$fit),
                         weights = FALSE)

# weights=FALSE per avere i valori degli outliers tra 0 e 1, senno sarebbero moltiplicati per il coeff X

## Plot weighted effects
par(mfrow=c(1,1))
plot(x = dati$TIME_PERIOD, main="effetti pesati delle anomalie", xlab="anno", ylab="effetto", y = rowSums(oeff), type = "l")
#### Estimate again
xreg <- as.matrix( oeff )
fit <- Arima(y = y,
             order = settings$order, seasonal = settings$seasonal,
             include.constant = settings$include.constant,
             xreg = xreg)
fit3 <- fit

# fit3 deve essere identico al modello stimato con tso() perch? ? quello che user? per le previsioni (tso non lo legge)



#####################################
## ARIMA + calendar effects
#####################################

#### Copy model
fit <- fit2
#### Extract settings
settings <- .Arima.settings(fit = fit)
#### Prepare xreg
xreg <- cal
#xreg <- if ( settings$include.drift ) { cbind(drift, xreg) } else { xreg }
#### Fit
fit <- tso(y = y, xreg = xreg,
           types = c("AO", "LS", "TC"), delta = 0.7, cval = 5,
           maxit = 10, maxit.iloop = 100, maxit.oloop = 10,
           # tsmethod = "auto.arima",
           # args.tsmethod = list(allowdrift = false, ic = "bic", trace = true) )
           tsmethod = "arima",
           args.tsmethod = list( order = settings$order, seasonal = settings$seasonal) )
fit2.o <- fit

#caso in cui il modello con out ed eff di calendario sia il migliore (non accade spesso)
#non ho variazioni nelle anomalie -> rimangono sempre uguali

### Reporting
print(fit)
plot(fit) #aggiungere legenda la blu ? la ts aggiustata
### Extract elements
## Outlier effects
oeff <- outliers.effects(mo = fit$outliers, n = NROW(y), pars = coef(fit$fit),
                         weights = FALSE)
## Plot weighted effects
par(mfrow=c(1,1))
plot(x = dati$TIME_PERIOD, main="effetti pesati delle anomalie", xlab="anno", ylab="effetto", y = rowSums(oeff), type = "l")

#### Estimate again
xreg <- cbind(cal, oeff)
fit <- Arima(y = y,
             order = settings$order, seasonal = settings$seasonal,
             include.constant = settings$include.constant,
             xreg = xreg)
fit4 <- fit

#trovo 3 anomalie in piu per? nel fit4 ho cval=4 e in fit3 ho cval=5
#porto in fondo fit4 per migliori diagnostiche e bic minore


#### Root analysis dei modelli con anomalie
#parte ar(p) con A(L)=0 e ma(q) con B(L)=0
par(mfrow = c(1,2))

#fit3
root <- .arma.roots(fit = fit3)
.circle(win = 2)
points(root$root$ar, col = "red")
points(root$root$ma, col = "blue")
legend("bottomright", inset= 0, title="ARIMA + Outliers", 
       legend=c("Radici della parte AR", "Radici della parte MA"),
       col=c("red", "blue"), pch=1, cex=.70)

#fit4
root <- .arma.roots(fit = fit4)
.circle(win = 2)
points(root$root$ar, col = "red")
points(root$root$ma, col = "blue")
legend("bottomright", inset= 0, title="ARIMA + Calendario + Outliers", 
       legend=c("Radici della parte AR", "Radici della parte MA"),
       col=c("red", "blue"), pch=1, cex=.70)


################################################################################
## Diagnostiche
################################################################################

#### Select the model fit4: modello con eff di calendario e anomalie
fit <- fit4

#### Useful quantities
npar1  <- NROW(fit$coef)                            ## Number of parameters
lag1   <- npar1 + c(1, 2, 5, 10, 15, 20)
res1   <- residuals(fit)                            ## Residuals
resst1 <- ( res1 - mean(res1) ) / sqrt(fit$sigma2)  ## Standardized residuals

#### Ts plot, acf, pacf, Ljung-Box of residuals
par(mfrow = c(2,1), cex.main=1.15)
main <- "Residui"
x1 <- res1
plot(x1, type = "l", main = main, ylab = "Residui")
Acf(x = x1, type = "correlation", lag.max = 60, na.action = na.pass, main = main)


# lb <- Box.test(x = res1, lag = 10, type = "Ljung-Box", fitdf = NROW(fit$coef))
cat("\n", paste("Ljung-Box dei", main, "a diversi lag\n") )
lb <- mapply(FUN = Box.test, lag = lag1, 
             MoreArgs = list(x = x1, type = "Ljung-Box", fitdf = npar1))[1:3, , drop = FALSE]
print(rbind(lag1,lb))

#### Ts plot, acf of residuals^2
par(mfrow = c(2,1))
main <- "Residui^2"
x1 <- res1^2
plot(x1, type = "l", main = main, ylab = "")
Acf(x = x1, type = "correlation", lag.max = 60, na.action = na.pass, main = main)

#Pacf inutile se i residui sono sempre piu vicini al WN perche acf simile pacf
#un WN ? formato da vc IID, quando facciamo l'acf e controlliamo che siano tutte nelle bande non stiamo controllando
#l 'indip ma l' incorrelazione seriale che sono cose diverse (indip->incorr ma non il contrario)
#se sono indip allora anche trasf dei residui sono incorrelate
#provo valori assoluti e res^2 che probabilmente saranno simili
#se vengono 1-2 lag fuori dalle bande puo andare bene lo stesso
#nel caso avessi lag fuori dalle bande (molti e\o alla fine dell acf) o mi va giu piano piano oppure 
#va giu veloce ma ho gobbe ai multipli del lag -> questo signifca che res non sono 
#indip e/o poi che sono eteroschedastici che pero si puo correggere questa eteroschedasticita 
#si fa arch test e se ce bisogno tsrf test per decidere quale trasformazione fare

#inutile fare LB test perche si vede gia dall acf che le cose non sono perfette ma posso usarlo per conferma

#### Ts plot, acf of |residuals|
par(mfrow = c(2,1), cex.main=1.15)
main <- "|Residui|"
x1 <- abs(res1)
plot(x1, type = "l", main = main, ylab = "")
Acf(x = x1, type = "correlation", lag.max = 60, na.action = na.pass, main = main)

#### Another diagnostic: the ARCH test
cat("\n-----------------------------------------------------------------
  ARCH based preliminary analyses\n")
cat("ARCH test on demeaned log-returns\n")
lag <- c(1, 2, 3, 6, 12, 24)
at <- mapply(FUN = ArchTest, lags = lag, 
             MoreArgs = list(x = x1 , demean = TRUE))
print(at[1:3,])
#demean vuol dire che levo la media ai residui


#Trasformazione
trt <- .trsf.test(fit = fit4, msg = "Transformation check on 'fit4'\n")



#### Distribuzione non condizionata dei residui
## Plot
par(mfrow = c(1,2))
hist(x = resst1, breaks = 25, freq = FALSE, main = "Residuals", xlab = "")
x1 <- seq(from = min(resst1), to = max(resst1)+1, length.out = 100) 
lines(x = x1, y = dnorm(x = x1, mean = 0, sd = 1), col = "red")
qqnorm(y = resst1, main = "Normal Q-Q Plot",
       xlab = "Theoretical quantiles", ylab = "Sample quantiles",
       plot.it = TRUE)
abline(a = 0, b = 1, col = "red")

## Test di normalità
print( shapiro.test(x = res1 ) )

#W = 0.97978, p-value = 3.464e-05 no normalità


################################################################################
## Forecasts
################################################################################

###########################################
## Ex-post forecasts: all made 1-step ahead
###########################################

#il controllo ex post non ? un controllo da prendere per un po pi? con fiducia le prev ex-ante
#? comunque un controllo su sole J oss pero meglio di nulla

#### Settings
J  <- 12                                              ## How many ex-post forecast to compute
## quelle per calcolare le misure di errore
H  <- 1                                               ## Forecasting horizon
## H=1 ? obbligatorio
#avendo 252 oss, nelle ex-ante parto a fare previsioni da 240 (H=12 -> 252-12=240)
#adesso da 239 (252-12-H=239) pero +1 =240
t1 <- .predict.t1(nobs = NROW(y), J = J, n.ahead = H) ## 1st obs in the ex-post period (needed below)
## calcolo su T oss, l'ultima oss che mi sta dentro It -> Y(T-j)

#nelle previsioni ex post gli S.E. non crescono ma sono tutte piu o meno uguali, sono tutte previsioni un passo avanti
#quindi l incertezza ? la stessa per tutte

#### No external regressors
pred1.1 <- .predict(object = fit1, n.ahead = H, t = t1, y = y,
                    fixed.n.ahead = TRUE)

#fixed.n.ahead = TRUE voglio solo la previsione all'orizzonte H fisso

#### If we have external regressors
newxreg <- cal
pred2.1 <- .predict(object = fit2, n.ahead = H, t = t1, y = y, xreg = newxreg,
                    fixed.n.ahead = TRUE)

#### If we have outliers
newxreg <- .oeff.4.predict(object = fit1.o, n.ahead = 0)
pred3.1 <- .predict(object = fit3, n.ahead = H, t = t1, y = y, xreg = newxreg,
                    fixed.n.ahead = TRUE)

# n.ahead=0 perch? non voglio aggiungere alcuna oss
# prima avevo =H perche volevo 12 oss in pi?

#### If we have external regressors and outliers
x2 <- .oeff.4.predict(object = fit2.o, n.ahead = 0)
newxreg <- as.matrix(cbind(cal, x2))
pred4.1 <- .predict(object = fit4, n.ahead = H, t = t1, y = y, xreg = newxreg,
                    fixed.n.ahead = TRUE)

#### Naive
predn.1 <- .predict.naive(fit = fit4, J = J, n.ahead = H)

#### Bands
x1 <- .pred.bands(pred = pred1.1, alpha = 0.05)
x2 <- .pred.bands(pred = pred2.1, alpha = 0.05)
x3 <- .pred.bands(pred = pred3.1, alpha = 0.05)
x4 <- .pred.bands(pred = pred4.1, alpha = 0.05)

#### Error Measures
em1.1  <- .ErrorMeasures(y = y, fit = x1$mean, naive = predn.1)
em2.1  <- .ErrorMeasures(y = y, fit = x2$mean, naive = predn.1)
em3.1  <- .ErrorMeasures(y = y, fit = x3$mean, naive = predn.1)
em4.1  <- .ErrorMeasures(y = y, fit = x4$mean, naive = predn.1)
emn.1  <- .ErrorMeasures(y = y, fit = predn.1, naive = predn.1)
## Print
ErrorMeas <- data.frame(
  model = c("Arima", "Arima + Calendar", "Arima + Outliers", "Arima + Calendar + Outliers", "Naive"),
  h = H,
  rbind( em1.1, em2.1, em3.1, em4.1, emn.1, deparse.level = 0 ) )
print( ErrorMeas )

#quale funziona meglio? guardo le misure scalate e prendo il modello che le ha piu basse
#quanto guadagnano i modelli?
#ScMEA e ScRMSE sono =1 ?

#### Plot
ind  <- (NROW(y) - J + 1) : NROW(y)
ind1 <- 1 : NROW(ind)
par(mfrow = c(1,1))
ylim <- range(
  x1$lower[ind1], x1$upper[ind1], x2$lower[ind1], x2$upper[ind1],
  x3$lower[ind1], x3$upper[ind1], x4$lower[ind1], x4$upper[ind1] )
time <- dati$TIME_PERIOD[ind]
plot(x = time, y = y[ind], ylim = c(60,125),
     main = "(Ex-post) Previsioni per i 12 mesi passati", xlab = "time", ylab = "Produzione")
lines(x = time, y = x1$mean[ind1],  col = "red")
lines(x = time, y = x2$mean[ind1],  col = "cyan")
lines(x = time, y = x3$mean[ind1],  col = "violet")
lines(x = time, y = x4$mean[ind1],  col = "blue")
lines(x = time, y = predn.1[ind1],  col = "green", lty = "solid")
lines(x = time, y = x4$lower[ind1], col = "blue", lty = "dotted")
lines(x = time, y = x4$upper[ind1], col = "blue", lty = "dotted")
legend("bottomleft", inset=0, title="Arima: ", 
       legend=c("Semplice","+ calendario", "+ anomalie", "+ calendario + outliers", "naive"), 
       col=c("red","cyan", "violet", "blue" ,"green"), lty=1, cex=0.75)

#il grafico ottenuto ? spesso simile al grafico ottenuto nelle ex-ante pero ho 2 varianti importanti:
#adesso ho dei dati (i punti) per giudicare se le previsioni vanno bene (su cui poi sovrapporre le previsione)
#nelle ex ante ho origine fissa e orizzonte mobile , mentre nel ex post ho orizzonte =1 sempre uguale pero 
#cambia l origine di modo che la distanza sia sempre 1 ed avere previsioni confrontabili


#######################################################
## Ex-ante (genuine) forecasts: from 1 to H steps ahead
#######################################################

#### Settings
H  <- 12      ## Maximum forecasting horizon
t1 <- NROW(y) ## Last obs in the info set (T)

#### Genero Time
## Useful for extending calendar and for plotting
time1 <- .extend.time(x = dati$TIME_PERIOD, n.ahead = H, by = "month")

#### Senza regressori esterni
newxreg <- NULL
pred1 <- .predict(object = fit1, n.ahead = H, t = t1, y = y, xreg = newxreg,
                  fixed.n.ahead = FALSE)

#fit1 ? il modello arima base senza nulla
#ha le radici tutte sul bordo di raggio 1 e quindi non va tanto bene
#n.ahead n di passi avanti
#t ? T
#fixed... se considero H mobile ovvero voglio FINO ad H
#$pred sono le previsioni con s.e. condizionati
#all'aumentare di H, le var di prev aumentano (ho gli se qua)
#l'se della prima previsione ? = sqrt(fit1$sigma2) simile 


#### Se ho regressori esterni
x1 <- .calendarEffects(time = time1)[, colnames(cal), drop = FALSE]
#mi fabbrico le v ar esterne che mi servono per la previsione
#time1 quello che sta fuori il periodo di previsione
newxreg <- as.matrix(rbind(cal, x1))
pred2 <- .predict(object = fit2, n.ahead = H, t = t1, y = y, xreg = newxreg,
                  fixed.n.ahead = FALSE)

#### Se ho outliers
x2 <- .oeff.4.predict(object = fit1.o, n.ahead = H)
newxreg <- x2
pred3 <- .predict(object = fit3, n.ahead = H, t = t1, y = y, xreg = newxreg,
                  fixed.n.ahead = FALSE)

#### Se ho regressori esterni e outliers
x2 <- .oeff.4.predict(object = fit2.o, n.ahead = H)
newxreg <- as.matrix( cbind( rbind(cal, x1), x2) )
pred4 <- .predict(object = fit4, n.ahead = H, t = t1, y = y, xreg = newxreg,
                  fixed.n.ahead = FALSE)

#### Naive
predn <- .predict.naive(fit = fit4, J = H, n.ahead = H) 

#### Bande

#in caso di trasf di var allora media e mediana potrebbero avere valori un po diversi fra loro
#differenza dovuta a:
#median ? trasformazione inversa (uso invarianza)
#mean valore atteso aggiunto del contributo della varianza
#per il grafico usaiamo la mean

x1 <- .pred.bands(pred = pred1, alpha = 0.05)
x2 <- .pred.bands(pred = pred2, alpha = 0.05)
x3 <- .pred.bands(pred = pred3, alpha = 0.05)
x4 <- .pred.bands(pred = pred4, alpha = 0.05)

#t time
#medie condiz e se condiz
#upper e lower bande

print( cbind(t = x1$t, pred1 = x1$mean, pred2 = x2$mean, pred3 = x3$mean, pred4 = x4$mean, 
             naive = predn) )
# non interessante da mettere nel report

#### Plot
#interessante da mettere nel report
par(mfrow = c(1,1))
ylim <- range(
  x1$lower, x1$upper, x2$lower, x2$upper,
  x3$lower, x3$upper, x4$lower, x4$upper, predn )
time <- time1
plot(x = time, y = x1$mean, type = "l", col = "red",
     main = "(Ex-ante) Previsioni per i prossimi 12 mesi", xlab = "time", ylab = "Produzione",
     ylim = ylim)
lines(x = time, y = x2$mean,  col = "cyan")
lines(x = time, y = x3$mean,  col = "violet")
lines(x = time, y = x4$mean,  col = "blue")
lines(x = time, y = x4$lower, col = "blue", lty = "dotted")
lines(x = time, y = x4$upper, col = "blue", lty = "dotted")
lines(x = time, y = predn, col = "green")
legend("bottomleft", inset=0, title="Arima: ", 
       legend=c("Semplice","+ calendario", "+ anomalie", "+ calendario + outliers", "naive"), 
       col=c("red","cyan", "violet", "blue" ,"green"), lty=1, cex=0.75)

#contando che il modello migliore sembrerebbe quello con calendario e outliers, guardando le previsioni
#ex ante e mettendole a confronto con le ex post si puo dire che le previsioni per i prossimi mesi sono buone

