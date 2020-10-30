remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# 
# SVARS_s_pcom <- function(name, v1, v2){
#   source("PS2_SVAR_Bootstrap_sin_precios_int.R")
#   assign(paste(deparse(substitute(name)), ".SIRF.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
#   assign(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE), envir = .GlobalEnv)
#   assign(paste(deparse(substitute(name)), ".FEVD.boot", sep=""), SVAR.fevd.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
#   assign(paste(deparse(substitute(name)), ".ERPT.boot", sep=""), SVAR.erpt.boot(name, Amat, Bmat, Yb, pmax, H_ERPT, v1, v2, a, R, cumulative = TRUE), envir = .GlobalEnv)
# }
# 
# SVARS<- function(name, v1, v2){
#   source("PS2_SVAR_Bootstrap.R")
#   assign(paste(deparse(substitute(name)), ".SIRF.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
#   assign(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE), envir = .GlobalEnv)
#   assign(paste(deparse(substitute(name)), ".FEVD.boot", sep=""), SVAR.fevd.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
#   assign(paste(deparse(substitute(name)), ".ERPT.boot", sep=""), SVAR.erpt.boot(name, Amat, Bmat, Yb, pmax, H_ERPT, v1, v2, a, R, cumulative = TRUE), envir = .GlobalEnv)
# }
source("PS2_SVAR_Plots.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("Prueba_vars and ggplot.R")

SVARS<- function(name, v1, v2){
  source("PS2_SVAR_Analysis.R")
  source("PS2_SVAR_Bootstrap.R")
  assign(paste(deparse(substitute(name)), ".SIRF.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".FEVD.boot", sep=""), SVAR.fevd.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".ERPT.boot", sep=""), SVAR.erpt.boot(name, Amat, Bmat, Yb, pmax, H_ERPT, v1, v2, a, R, cumulative = TRUE), envir = .GlobalEnv)
  
}

SVARS_s_pcom<- function(name, v1, v2){
  source("PS2_SVAR_Analysis.R")
  source("PS2_SVAR_Bootstrap_sin_precios_int.R")
  assign(paste(deparse(substitute(name)), ".SIRF.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".FEVD.boot", sep=""), SVAR.fevd.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".ERPT.boot", sep=""), SVAR.erpt.boot(name, Amat, Bmat, Yb, pmax, H_ERPT, v1, v2, a, R, cumulative = TRUE), envir = .GlobalEnv)
}

graficos <- function(name, width, height){
  
  png(file=paste(deparse(substitute(name)), ".SIRF.boot", ".png", sep=""),  width= width, height= height)
  plot1 <- plot.sirf.boot(eval(as.name(paste(deparse(substitute(name)), ".SIRF.boot", sep=""))), m, H)
  print(plot1)
  dev.off()
  
  
  png(file=paste(deparse(substitute(name)), ".SIRF.c.boot", ".png", sep=""), width= width, height= height)
  plot2 <- plot.sirf.boot(eval(as.name(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""))), m, H)
  print(plot2)
  dev.off()
  
  png(file=paste(deparse(substitute(name)), ".FEVD.boot", ".png", sep=""), width= width, height= height)
  plot3 <- plot.fevd.boot(eval(as.name(paste(deparse(substitute(name)), ".FEVD.boot", sep=""))), m, H)
  print(plot3)
  dev.off()
  
  png(file=paste(deparse(substitute(name)), ".ERPT.boot", ".png", sep=""),  width=800, height=600)
  plot4 <- grafico(eval(as.name(paste(deparse(substitute(name)), ".ERPT.boot", sep=""))), H_ERPT)
  print(plot4)
  dev.off()
  return (list(plot1, plot2, plot3, plot4))
}



# Data ####

source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

#Seleccionamos la ventana temporal. Enero de 2005 a Diciembre de 2019.
Yl <- window(Yl.f, start = c(2005, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2005, 01), end = c(2019, 12))


# Inciso 1 ####

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
#MARTIN: No entiendo muy bien por qué, pero solo eliminamos 10 valores, no 12. 
#MARTIN: ¿No deberíamos comernos únicamente los primeros 3 valores?
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation
VAR <- VAR(Ydt, p = p, type = "const")

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs # No. of effective sample observations, excluding "p" starting values

# Ad hoc function
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}

# Re-estimate VAR (no feedback from local vars. to pcom)
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
VAR


# Model checking
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation

# A Matrix
#esta representa a la matriz A0^(-1)
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# SVAR estimation (AB model configuration)
#aca simplemente estimamos el SVAR imponiendo ciertas assumptions embebidas en A y en B
#los NA en A y en B se van a estimar como en la slide 24 de la lecture 4
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR


H <- 18 # Horizon
H_ERPT <- 120 # Horizon for ERPT
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications

#Reporte de resultados inciso 1
#(quiza podriamos hacer graficos mas lindos con ggplot)
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

SVARS(SVAR, 3, 2)
graficos(SVAR, 800, 800)
# 
# # IRF
# #SVAR.SIRF <- SVAR.sirf(SVAR, H)
# plot.sirf.boot(SVAR.SIRF.boot, m, H)
# 
# # Cumulative IRF
# #SVAR.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
# plot.sirf.boot(SVAR.SIRF.boot.c, m, H)
# 
# # FEVD
# #SVAR.FEVD <- SVAR.fevd(SVAR, H)
# plot.fevd.boot(SVAR.FEVD.boot, m, H)
# 
# # HD
# #SVAR.HD <- SVAR.hd(SVAR)
# #plot.hd.boot(Yd, SVAR.HD.boot, m, pmax)
# 
# # ERPT
# #SVAR.ERPT.boot <- SVAR.erpt(SVAR, H_ERPT, 3, 2, cumulative = TRUE)
# plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)
# 
# #PENDIENTE: Revisar cuestiones de estacionariedad de las series.
# 

# Inciso 2 ####
#Estimamos el modelo sin el índice de precios externos.

Yd_2 <- Yd[, 2:3]
popt_2 <- VARselect(Yd_2, lag.max = pmax, type = "const")
popt_2
p_2 <- popt_2$selection[1] # AIC
p_2 <- 1

Yd0_2 <- Yd_2[1:pmax, ] # Initial values
Ydt_2 <- Yd_2[(pmax - p_2 + 1):nrow(Yd_2), ] 

# Estimation
VAR_2 <- VAR(Ydt_2, p = p_2, type = "const")

m_2 <- VAR_2$K # No. of variables in the VAR
N_2 <- VAR_2$obs

roots(VAR_2, modulus = TRUE)
serial.test(VAR_2, lags.bg = 12, type = "ES")

# SVAR estimation

# A Matrix
Amat_2 <- diag(m_2)
for (i in 2:m_2) {
  for (j in 1:(i - 1)) {
    Amat_2[i, j] <- NA
  }
}

# B Matrix
Bmat_2 <- matrix(0, m_2, m_2)
for (i in 1:m_2) {
  Bmat_2[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_2 <- SVAR(VAR_2, Amat = Amat_2, Bmat = Bmat_2, lrtest = FALSE)
SVAR_2

#Reportamos resultados modelo inciso 2.

#hay que volver a llamarla m porque las funciones usadas a continuación están definidas usando m, y la cantidad de variables no se puede especificar de otra manera
#después de haber calculado las cosas que queremos volvemos a darle el valor inicial.
m <- m_2
Yb <- boot.rb.replicate(VAR_2, Yd0_2,  pmax, R)

#SVARS_s_pcom(SVAR_2, 2, 1)

 SVAR_2.SIRF.boot <- SVAR.sirf.boot(SVAR_2, Amat_2, Bmat_2, Yb, pmax, H, a, R)
 png(file="SVAR_2.sirf.boot.png",  width=800, height=800)
 plot.sirf.boot(SVAR_2.SIRF.boot, m, H)
 dev.off()
 
# 
# # Cumulative IRF (bootstrap)
 SVAR_2.SIRF.c.boot <- SVAR.sirf.boot(SVAR_2, Amat_2, Bmat_2, Yb, pmax, H, a, R, cumulative = TRUE)
 png(file="SVAR_2.sirf.c.boot.png",   width=800, height=800)
 plot.sirf.boot(SVAR.SIRF.c.boot, m, H)
 dev.off()
 
# 
# # FEVD (bootstrap)
 SVAR_2.FEVD.boot <- SVAR.fevd.boot(SVAR_2, Amat_2, Bmat_2, Yb, pmax, H, a, R)
 png(file="SVAR_2.FEVD.boot.png",   width=800, height=800)
 plot.fevd.boot(SVAR_2.FEVD.boot, m, H)
 dev.off()
# 
# # ERPT (bootstrap)
SVAR_2.ERPT.boot <- SVAR.erpt.boot(SVAR_2, Amat_2, Bmat_2, Yb, pmax, H_ERPT, 2, 1, a, R, cumulative = TRUE)
png(file="SVAR_2.ERPT.boot.png",   width=800, height=600)
 grafico(SVAR_2.ERPT.boot, H_ERPT)
 dev.off()
# 
# 
# # IRF
# SVAR_2.SIRF.boot <- SVAR.sirf(SVAR_2, H)
# plot.sirf.boot(SVAR_2.SIRF.boot, m, H)
# 
# # Cumulative IRF
# #SVAR_2.SIRF.c <- SVAR.sirf(SVAR_2, H, cumulative = TRUE)
# plot.sirf.boot(SVAR_2.SIRF.c.boot, m, H)
# 
# # FEVD
# #SVAR_2.FEVD <- SVAR.fevd(SVAR_2, H)
# plot.fevd.boot(SVAR_2.FEVD.boot, m, H)

# HD
#SVAR_2.HD <- SVAR.hd(SVAR_2)
#plot.hd(Yd_2, SVAR_2.HD.boot, m, pmax)

# ERPT
#SVAR_2.ERPT <- SVAR.erpt(SVAR_2, H_ERPT, 2, 1, cumulative = TRUE)
plot.erpt.boot(SVAR_2.ERPT.boot, H_ERPT)

#Pendiente: comparar con resultados modelo inciso 1.

# Punto 3 #### 
#run each subsample and repeat using 3 lags

#sub-sample 1 ####: 
#from January 2005 to July 2011


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("PS1_Data.R")
library(vars)
Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) 
Yd.f <- 100 * diff(Yl.f)

Yl <- window(Yl.f, start = c(2003, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2003, 01), end = c(2019, 12))

Yone <- window(Yd.f, start = c(2004, 01), end = c(2011, 07))

pmax <- 12

popt <- VARselect(Yone, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]

Yone0 <- Yone[1:pmax, ] # Initial values, 12M
Yonet <- Yone[(pmax - p + 1):nrow(Yone), ]

VAR <- VAR(Yonet, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Yone0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

#sub-sample 2 ####
#from August 2011 to September 2015

remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) 
Yd.f <- 100 * diff(Yl.f)

Yl <- window(Yl.f, start = c(2003, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2003, 01), end = c(2019, 12))

pmax <- 12

Ytwo <- window(Yd.f, start = c(2010, 08), end = c(2015, 09))

popt <- VARselect(Ytwo, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]


Ytwo0 <- Ytwo[1:pmax, ] # Initial values, 12M
Ytwot <- Ytwo[(pmax - p + 1):nrow(Ytwo), ]

VAR <- VAR(Ytwot, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)


P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Ytwo0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# sub-sample 3
####from October 2015 to August 2019

remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) 
Yd.f <- 100 * diff(Yl.f)

Yl <- window(Yl.f, start = c(2003, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2003, 01), end = c(2019, 12))

pmax <- 12

Ythree <- window(Yd.f, start = c(2014, 10), end = c(2019, 08))

popt <- VARselect(Ythree, lag.max = pmax, type = "const")
popt
p <- popt$selection[2] # [check and change for each sample} /// [change again if autocorrelated errors]

Ythree0 <- Ythree[1:pmax, ] # Initial values, 12M
Ythreet <- Ythree[(pmax - p + 1):nrow(Ythree), ]

VAR <- VAR(Ythreet, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Ythree0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)


# Inciso 4 #### 
#Funcion para ejecutar todos las funciones sobre el SVAR automaticamente.
#El 1° argumento tiene que ser el nombre con el que uno quiere que terminen los outputs.


# Parte A: calcular desde 2013 hasta 2019, usando IPC general.
# Parte A.1: sin precios internacionales
#Lo seteo en 11 porque para el IPC CABA no tenemos la variación de julio de 2012, con lo que no tenemos 12 meses como condición inicial, sino 11.
#Por la misma razón tomo desde agosto de 2012 en lugar de desde julio.


source("Prueba_vars and ggplot.R")
source("PS2_SVAR_Plots.R")


SVARS<- function(name, v1, v2){
  source("PS2_SVAR_Analysis.R")
  source("PS2_SVAR_Bootstrap.R")
  assign(paste(deparse(substitute(name)), ".SIRF.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".FEVD.boot", sep=""), SVAR.fevd.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".ERPT.boot", sep=""), SVAR.erpt.boot(name, Amat, Bmat, Yb, pmax, H_ERPT, v1, v2, a, R, cumulative = TRUE), envir = .GlobalEnv)
  
}

SVARS_s_pcom<- function(name, v1, v2){
  source("PS2_SVAR_Analysis.R")
  source("PS2_SVAR_Bootstrap_sin_precios_int.R")
  assign(paste(deparse(substitute(name)), ".SIRF.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""), SVAR.sirf.boot(name, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".FEVD.boot", sep=""), SVAR.fevd.boot(name, Amat, Bmat, Yb, pmax, H, a, R), envir = .GlobalEnv)
  assign(paste(deparse(substitute(name)), ".ERPT.boot", sep=""), SVAR.erpt.boot(name, Amat, Bmat, Yb, pmax, H_ERPT, v1, v2, a, R, cumulative = TRUE), envir = .GlobalEnv)
}

graficos <- function(name, width, height){
  
  png(file=paste(deparse(substitute(name)), ".SIRF.boot", ".png", sep=""),  width= width, height= height)
  plot1 <- plot.sirf.boot(eval(as.name(paste(deparse(substitute(name)), ".SIRF.boot", sep=""))), m, H)
  print(plot1)
  dev.off()
  
  
  png(file=paste(deparse(substitute(name)), ".SIRF.c.boot", ".png", sep=""), width= width, height= height)
  plot2 <- plot.sirf.boot(eval(as.name(paste(deparse(substitute(name)), ".SIRF.c.boot", sep=""))), m, H)
  print(plot2)
  dev.off()
  
  png(file=paste(deparse(substitute(name)), ".FEVD.boot", ".png", sep=""),  width= width, height= height)
  plot3 <- plot.fevd.boot(eval(as.name(paste(deparse(substitute(name)), ".FEVD.boot", sep=""))), m, H)
  print(plot3)
  dev.off()
  
  png(file=paste(deparse(substitute(name)), ".ERPT.boot", ".png", sep=""),  width=800, height=600)
  plot4 <- grafico(eval(as.name(paste(deparse(substitute(name)), ".ERPT.boot", sep=""))), H_ERPT)
  print(plot4)
  dev.off()
  return (list(plot1, plot2, plot3, plot4))
}



pmax <- 11
Yd_4 <- window(Yd.f, start = c(2012, 08), end = c(2019, 12))
Yd_4_s_pcom <- Yd_4[, 2:3]
popt_4_s_pcom <- VARselect(Yd_4_s_pcom, lag.max = pmax, type = "const")
popt_4_s_pcom
p_4_s_pcom <- popt_4_s_pcom$selection[1] # AIC

Yd0_4_s_pcom <- Yd_4_s_pcom[1:pmax, ] # Initial values
Ydt_4_s_pcom <- Yd_4_s_pcom[(pmax - p_4_s_pcom + 1):nrow(Yd_4_s_pcom), ]

# Estimation
VAR_4_s_pcom <- VAR(Ydt_4_s_pcom, p = p_4_s_pcom, type = "const")

m_4_s_pcom <- VAR_4_s_pcom$K # No. of variables in the VAR
N_4_s_pcom <- VAR_4_s_pcom$obs

N <- N_4_s_pcom
m <- m_4_s_pcom

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

roots(VAR_4_s_pcom, modulus = TRUE)
serial.test(VAR_4_s_pcom, lags.bg = 12, type = "ES")

# SVAR estimation (AB model configuration)
SVAR_4_s_pcom <- SVAR(VAR_4_s_pcom, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR_4_s_pcom

#Resultados parte A.1.
m <- m_4_s_pcom
N <- N_4_s_pcom 


H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR_4_s_pcom, Yd0_4_s_pcom,  pmax, R)

SVARS_s_pcom(SVAR_4_s_pcom, 2, 1)
graficos(SVAR_4_s_pcom)
# 
# # IRF (bootstrap)
# #SVAR_4_s_pcom.SIRF.boot <- SVAR.sirf.boot(SVAR_4_s_pcom, Amat, Bmat, Yb, pmax, H, a, R)
# plot.sirf.boot(SVAR_4_s_pcom.SIRF.boot, m, H)
# 
# # Cumulative IRF (bootstrap)
# #SVAR_4_s_pcom.SIRF.c.boot <- SVAR.sirf.boot(SVAR_4_s_pcom, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
# plot.sirf.boot(SVAR_4_s_pcom.SIRF.c.boot, m, H)
# 
# # FEVD (bootstrap)
# #SVAR_4_s_pcom.FEVD.boot <- SVAR.fevd.boot(SVAR_4_s_pcom, Amat, Bmat, Yb, pmax, H, a, R)
# plot.fevd.boot(SVAR_4_s_pcom.FEVD.boot, m, H)
# 
# #MARTIN: Problema con subscript, revisar.
# # ERPT (bootstrap)
# #SVAR_4_s_pcom.ERPT.boot <- SVAR.erpt.boot(SVAR_4_s_pcom, Amat, Bmat, Yb, pmax, H_ERPT, 2, 1, a, R, cumulative = TRUE)
# plot.erpt.boot(SVAR_4_s_pcom.ERPT.boot, H_ERPT)


#Parte A.2: con precios internacionales
#Es 11 por la misma razón que en la parte anterior.
pmax <- 11

popt_4_c_pcom <- VARselect(Yd_4, lag.max = pmax, type = "const")
popt_4_c_pcom
p_4_c_pcom <- popt_4_c_pcom$selection[1] # AIC

Yd0_4_c_pcom <- Yd_4[1:pmax, ] # Initial values
#MARTIN: No entiendo muy bien por qué, pero solo eliminamos 10 valores, no 12. 
#MARTIN: ¿No deberíamos comernos únicamente los primeros 3 valores?
Ydt_4_c_pcom <- Yd_4[(pmax - p_4_c_pcom + 1):nrow(Yd_4), ] # Starting in Jan-04

# Estimation
VAR_4_c_pcom <- VAR(Ydt_4_c_pcom, p = p_4_c_pcom, type = "const")

m_4_c_pcom <- VAR_4_c_pcom$K # No. of variables in the VAR
N_4_c_pcom <- VAR_4_c_pcom$obs # No. of effective sample observations, excluding "p" starting values
m <- m_4_c_pcom
N <- N_4_c_pcom 


Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# Ad hoc function
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}

# Re-estimate VAR (no feedback from local vars. to pcom)
VAR_4_c_pcom <- restrict(VAR_4_c_pcom, method = "man", resmat = matC(m_4_c_pcom, p_4_c_pcom, 1))
VAR_4_c_pcom


# Model checking
roots(VAR_4_c_pcom, modulus = TRUE)
serial.test(VAR_4_c_pcom, lags.bg = 12, type = "ES")

# SVAR estimation (AB model configuration)
#aca simplemente estimamos el SVAR imponiendo ciertas assumptions embebidas en A y en B
#los NA en A y en B se van a estimar como en la slide 24 de la lecture 4
SVAR_4_c_pcom <- SVAR(VAR_4_c_pcom, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

#Reporte de resultados inciso 4 parte a2
#(quiza podriamos hacer graficos mas lindos con ggplot)

Yb <- boot.rb.replicate(VAR_4_c_pcom, Yd0_4_c_pcom,  pmax, R)
source("PS2_SVAR_Bootstrap.R")
SVARS(SVAR_4_c_pcom, 3, 2)
graficos(SVAR_4_c_pcom)

# 
# # IRF (bootstrap)
# #SVAR_4_c_pcom.SIRF.boot <- SVAR.sirf.boot(SVAR_4_c_pcom, Amat, Bmat, Yb, pmax, H, a, R)
# plot.sirf.boot(SVAR_4_c_pcom.SIRF.boot, m, H)
# 
# # Cumulative IRF (bootstrap)
# #SVAR_4_c_pcom.SIRF.c.boot <- SVAR.sirf.boot(SVAR_4_c_pcom, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
# plot.sirf.boot(SVAR_4_c_pcom.SIRF.c.boot, m, H)
# 
# # FEVD (bootstrap)
# #SVAR_4_c_pcom.FEVD.boot <- SVAR.fevd.boot(SVAR_4_c_pcom, Amat, Bmat, Yb, pmax, H, a, R)
# plot.fevd.boot(SVAR_4_c_pcom.FEVD.boot, m, H)
# 
# # ERPT (bootstrap)
# #SVAR_4_c_pcom.ERPT.boot <- SVAR.erpt.boot(SVAR_4_c_pcom, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
# plot.erpt.boot(SVAR_4_c_pcom.ERPT.boot, H_ERPT)

# 
# # IRF
# SVAR_4_c_pcom.SIRF <- SVAR.sirf(SVAR_4_c_pcom, H)
# plot.sirf(SVAR_4_c_pcom.SIRF, m, H)
# 
# # Cumulative IRF
# SVAR_4_c_pcom.SIRF.c <- SVAR.sirf(SVAR_4_c_pcom, H, cumulative = TRUE)
# plot.sirf(SVAR_4_c_pcom.SIRF.c, m, H)
# 
# # FEVD
# SVAR_4_c_pcom.FEVD <- SVAR.fevd(SVAR_4_c_pcom, H)
# plot.fevd(SVAR_4_c_pcom.FEVD, m, H)
# 
# # HD
# SVAR_4_c_pcom.HD <- SVAR.hd(SVAR_4_c_pcom)
# plot.hd(Yd, SVAR_4_c_pcom.HD, m, pmax)

# ERPT
# #SVAR_4_c_pcom.ERPT <- SVAR.erpt(SVAR_4_c_pcom, H_ERPT, 3, 2, cumulative = TRUE)
# plot.erpt(SVAR_4_c_pcom.ERPT, H_ERPT)


#Parte B: con índice de precios de CABA.
#B.1: Sin precios internacionales.
# CPI, City of Buenos Aires
ipc.ba <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?ids=193.1_NIVEL_GENERAL_JULI_0_13&limit=5000&format=csv"))
ipc.ba <- ts(ipc.ba$nivel_general, start = c(2012, 07), frequency = 12)
ipc.ba <- diff(log(ipc.ba))
ipc.ba <- window(ipc.ba, end = c(2019, 12))
ipc.ba <- window(ipc.ba, start = c(2012, 08))

Yd_caba <- Yd_4[, 1:2]
Yd_caba <- cbind(Yd_caba, ipc.ba)
Yd_caba_s_com <- Yd_caba[, 2:3]

popt_4_s_pcom <- VARselect(Yd_caba_s_com, lag.max = pmax, type = "const")
popt_4_s_pcom
p_4_s_pcom <- popt_4_s_pcom$selection[1] # AIC

Yd0_4_s_pcom <- Yd_caba_s_com [1:pmax, ] # Initial values
Ydt_4_s_pcom <- Yd_caba_s_com [(pmax - p_4_s_pcom + 1):nrow(Yd_4_s_pcom), ]

# Estimation
VAR_4_s_pcom <- VAR(Ydt_4_s_pcom, p = p_4_s_pcom, type = "const")

m <- VAR_4_s_pcom$K # No. of variables in the VAR
N <- VAR_4_s_pcom$obs

m_4_s_pcom_caba <-m
N_4_s_pcom_caba <- N

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

roots(VAR_4_s_pcom, modulus = TRUE)
serial.test(VAR_4_s_pcom, lags.bg = 12, type = "ES")

# SVAR estimation (AB model configuration)
SVAR_4_s_pcom_caba <- SVAR(VAR_4_s_pcom, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR <- SVAR_4_s_pcom_caba

#Resultados parte B.1.

Yb <- boot.rb.replicate(VAR_4_s_pcom, Yd0_4_s_pcom,  pmax, R)
source("PS2_SVAR_Bootstrap_sin_precios_int.R")
SVARS_s_pcom(SVAR_4_s_pcom_caba, 2, 1)
graficos(SVAR_4_s_pcom_caba)

# IRF (bootstrap)
# 
# plot.sirf.boot(SVAR_4_s_pcom_caba.SIRF.boot, m, H)
# 
# # Cumulative IRF (bootstrap)
# plot.sirf.boot(SVAR_4_s_pcom_caba.SIRF.c.boot, m, H)
# 
# # FEVD (bootstrap)
# plot.fevd.boot(SVAR_4_s_pcom_caba.FEVD.boot, m, H)
# 
# #MARTIN: Problema con subscript, revisar.
# # ERPT (bootstrap)
# plot.erpt.boot(SVAR_4_s_pcom_caba.ERPT.boot, H_ERPT)

# 
# 
# # IRF
# SVAR_4_s_pcom_caba.SIRF <- SVAR.sirf(SVAR_4_s_pcom_caba, H)
# plot.sirf(SVAR_4_s_pcom_caba.SIRF, m, H)
# 
# # Cumulative IRF
# SVAR_4_s_pcom_caba.SIRF.c <- SVAR.sirf(SVAR_4_s_pcom_caba, H, cumulative = TRUE)
# plot.sirf(SVAR_4_s_pcom_caba.SIRF.c, m, H)
# 
# # FEVD
# SVAR_4_s_pcom_caba.FEVD <- SVAR.fevd(SVAR_4_s_pcom_caba, H)
# plot.fevd(SVAR_4_s_pcom_caba.FEVD, m, H)
# 
# # ERPT
# SVAR_4_s_pcom_caba.ERPT <- SVAR.erpt(SVAR_4_s_pcom_caba, H_ERPT, 2, 1, cumulative = TRUE)
# plot.erpt(SVAR_4_s_pcom_caba.ERPT, H_ERPT)


#B.2: Con precios internacionales.
Yd_caba_c_com <- Yd_caba[, 2:3]
popt_4_c_pcom <- VARselect(Yd_caba_c_com, lag.max = pmax, type = "const")
popt_4_c_pcom

p_4_c_pcom <- popt_4_c_pcom$selection[1] # AIC

Yd0_4_c_pcom <- Yd_caba_c_com[1:pmax, ] # Initial values

Ydt_4_c_pcom <- Yd_caba_c_com[(pmax - p_4_c_pcom + 1):nrow(Yd_caba_c_com), ] # Starting in Jan-04

# Estimation
VAR_4_c_pcom <- VAR(Ydt_4_c_pcom, p = p_4_c_pcom, type = "const")

m_4_c_pcom_caba <- VAR_4_c_pcom$K # No. of variables in the VAR
N_4_c_pcom_caba <- VAR_4_c_pcom$obs # No. of effective sample observations, excluding "p" starting values
m <- m_4_c_pcom_caba
N <- N_4_c_pcom_caba 
# Ad hoc function

matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}

# Re-estimate VAR (no feedback from local vars. to pcom)
VAR_4_c_pcom <- restrict(VAR_4_c_pcom, method = "man", resmat = matC(m, p, 1))
VAR_4_c_pcom


# Model checking
roots(VAR_4_c_pcom, modulus = TRUE)
serial.test(VAR_4_c_pcom, lags.bg = 12, type = "ES")

# SVAR estimation


Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}


# SVAR estimation (AB model configuration)
#aca simplemente estimamos el SVAR imponiendo ciertas assumptions embebidas en A y en B
#los NA en A y en B se van a estimar como en la slide 24 de la lecture 4
SVAR_4_c_pcom_caba <- SVAR(VAR_4_c_pcom, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
Yb <- boot.rb.replicate(VAR_4_c_pcom, Yd0_4_c_pcom,  pmax, R)

SVARS(SVAR_4_c_pcom_caba, 3, 2)
graficos(SVAR_4_c_pcom_caba)
#Reporte de resultados inciso 4 parte a2
#(quiza podriamos hacer graficos mas lindos con ggplot)

# IRF
plot.sirf.boot(SVAR_4_c_pcom_caba.SIRF.boot, m, H)

# Cumulative IRF
plot.sirf.boot(SVAR_4_c_pcom_caba.SIRF.c.boot, m, H)

# FEVD
plot.fevd.boot(SVAR_4_c_pcom_caba.FEVD.boot, m, H)

# HD
#plot.hd(Yd, SVAR_4_c_pcom_caba.HD.boot, m, pmax)

# ERPT
plot.erpt.boot(SVAR_4_c_pcom_caba.ERPT.boot, H_ERPT)


# Inciso 5 ####

wages.file <- paste(tempfile(), ".csv", sep = "")
download.file("https://infra.datos.gob.ar/catalog/sspm/dataset/157/distribution/157.2/download/remuneracion-bruta-promedio-asalariados-registrados-sector-privado-segun-actividad-sin-estacionalidad.csv
", wages.file, mode = "wb")

wages <- read.csv(wages.file)
wages <- wages[c(1,16)]
wages_num <- as.numeric(wages$rem_bruta_asal_sin_est_total)

wages_num <- wages_num[complete.cases(wages_num)] # delete NAs
wages_num <- ts(wages_num, start = c(1995, 01), end = c(2020, 03), frequency = 12)
wages_num  <- window(wages_num, start = c(2002, 03))

wages_real <- wages_num * tail(pc, n=1)/pc
wages_real <- wages_real / tail(wages_real, n=1)
plot(wages_real)
wages_real  <- window(wages_real, start = c(2004, 12))

var_wages_real <- 100 * diff(log(wages_real))  # log transformation

Yd_5_1 <- cbind(Yd, var_wages_real) # Raw data in log
Yd_5 <- window(Yd_5_1,start = c(2005, 01), end = c(2019, 12))

#Parte A: sin precios internacionales

pmax <- 12
Yd_5_s_pcom <- Yd_5[, 2:4]
popt <- VARselect(Yd_5_s_pcom, lag.max = pmax, type = "const")
p <- popt$selection[1] # AIC
p_5_s_pcom <- p

Yd0 <- Yd_5_s_pcom [1:pmax, ] # Initial values
Ydt <- Yd_5_s_pcom [(pmax - p + 1):nrow(Yd_5_s_pcom), ]

# Estimation
VAR <- VAR(Ydt, p = p, type = "const")

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs

m_5_s_pcom <- m
N_5_s_pcom <- N

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation (AB model configuration)
SVAR_5_s_pcom <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR_5_s_pcom

#Resultados parte A

Yb <- boot.rb.replicate(VAR, Yd0,  pmax, R)
SVARS_s_pcom(SVAR_5_s_pcom, 2, 1)
graficos(SVAR_5_s_pcom)
# 
# # IRF
# #SVAR_5_s_pcom.SIRF <- SVAR.sirf(SVAR, H)
# plot.sirf.boot(SVAR_5_s_pcom.SIRF.boot, m, H)
# 
# # Cumulative IRF
# #SVAR_5_s_pcom.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
# plot.sirf.boot(SVAR_5_s_pcom.SIRF.c.boot, m, H)
# 
# # FEVD
# #SVAR_5_s_pcom.FEVD <- SVAR.fevd(SVAR, H)
# plot.fevd.boot(SVAR_5_s_pcom.FEVD.boot, m, H)
# 
# # ERPT
# #SVAR_5_s_pcom.ERPT <- SVAR.erpt(SVAR, H_ERPT, 2, 1, cumulative = TRUE)
# plot.erpt.boot(SVAR_5_s_pcom.ERPT.boot, H_ERPT)

#Parte B: con precios internacionales

popt <- VARselect(Yd_5, lag.max = pmax, type = "const")

p <- popt$selection[1] # AIC
p_5_c_pcom<- p

Yd0 <- Yd_5[1:pmax, ] # Initial values
#MARTIN: No entiendo muy bien por qué, pero solo eliminamos 10 valores, no 12. 
#MARTIN: ¿No deberíamos comernos únicamente los primeros 3 valores?
Ydt <- Yd_5[(pmax - p + 1):nrow(Yd_5), ] # Starting in Jan-04

# Estimation
VAR <- VAR(Ydt, p = p, type = "const")

m_5_c_pcom <- VAR$K # No. of variables in the VAR
N_5_c_pcom <- VAR$obs # No. of effective sample observations, excluding "p" starting values
m <- m_5_c_pcom
N <- N_5_c_pcom 


Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# Ad hoc function
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}

# Re-estimate VAR (no feedback from local vars. to pcom)
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
VAR


# Model checking
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation (AB model configuration)
#aca simplemente estimamos el SVAR imponiendo ciertas assumptions embebidas en A y en B
#los NA en A y en B se van a estimar como en la slide 24 de la lecture 4
SVAR_5_c_pcom <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

#Reporte de resultados inciso 5

Yb <- boot.rb.replicate(VAR, Yd0,  pmax, R)
SVARS(SVAR_5_c_pcom, 3, 2)
graficos(SVAR_5_c_pcom)
# 
# # IRF
# #SVAR_5_c_pcom.SIRF <- SVAR.sirf(SVAR, H)
# plot.sirf.boot(SVAR_5_c_pcom.SIRF.boot, m, H)
# 
# # Cumulative IRF
# #SVAR_5_c_pcom.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
# plot.sirf.boot(SVAR_5_c_pcom.SIRF.c.boot, m, H)
# 
# # FEVD
# #SVAR_5_c_pcom.FEVD <- SVAR.fevd(SVAR, H)
# plot.fevd.boot(SVAR_5_c_pcom.FEVD.boot, m, H)
# 
# # HD
# #SVAR_5_c_pcom.HD <- SVAR.hd(SVAR)
# #plot.hd.boot(Yd, SVAR_5_c_pcom.HD.boot, m, pmax)
# 
# # ERPT
# #SVAR_5_c_pcom.ERPT <- SVAR.erpt(SVAR, H_ERPT, 3, 2, cumulative = TRUE)
# plot.erpt.boot(SVAR_5_c_pcom.ERPT.boot, H_ERPT)



#Inciso 6 ####

#Necesario instalar el siguiente paquete:
#util para crear dummies temporales
#install.packages("tstools")
library(tstools)

emae.file <- paste(tempfile(), ".csv", sep = "")
download.file("https://infra.datos.gob.ar/catalog/sspm/dataset/143/distribution/143.3/download/emae-valores-anuales-indice-base-2004-mensual.csv", emae.file, mode = "wb")
emae <- read.csv(emae.file)
emae <- emae[c(1,6)]
emae_num <- as.numeric(emae$emae_desestacionalizada_vm)
emae_num <- ts(emae_num, start = c(2004, 01), end = c(2020, 07), frequency = 12)
emae_num <- window(emae_num, start = c(2005, 01), end = c(2019, 12))*100
plot(emae_num)
emae_num <- emae_num[complete.cases(emae_num)] # delete NAs
Yd_6 <- cbind(Yd_5, emae_num) 

anio <-2009

d2009_1<- create_dummy_ts(end_basic = c(anio,06), dummy_start = c(anio,03), dummy_end = c(anio,06), sp = TRUE, start_basic = c(2005,1),  dummy_value = -1, frequency = 12)
d2009_2<- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(anio,08), dummy_end = c(anio,8), sp = TRUE, start_basic = c(anio,7),  dummy_value = 1, frequency = 12)
d2009_2[1] = 1
d2009 <- ts(c(d2009_1,d2009_2),  start = start(d2009_2),frequency = frequency(121))
d2009 <- ts(d2009, start = c(2005, 01), end = c(2019, 12), frequency = 12)

anio <-2012

d2012_1<- create_dummy_ts(end_basic = c(anio,06), dummy_start = c(anio,03), dummy_end = c(anio,06), sp = TRUE, start_basic = c(2005,1),  dummy_value = -1, frequency = 12)
d2012_2<- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(anio,08), dummy_end = c(anio,8), sp = TRUE, start_basic = c(anio,7),  dummy_value = 1, frequency = 12)
d2012_2[1] = 1
d2012 <- ts(c(d2012_1,d2012_2),  start = start(d2012_2),frequency = frequency(121))
d2012 <- ts(c(d2012_1,d2012_2),  start = start(d2012_2),frequency = frequency(121))
d2012 <- ts(d2012, start = c(2005, 01), end = c(2019, 12), frequency = 12)

anio <-2018

d2018_1<- create_dummy_ts(end_basic = c(anio,06), dummy_start = c(anio,03), dummy_end = c(anio,06), sp = TRUE, start_basic = c(2005,1),  dummy_value = -1, frequency = 12)
d2018_2<- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(anio,08), dummy_end = c(anio,8), sp = TRUE, start_basic = c(anio,7),  dummy_value = 1, frequency = 12)
d2018_2[1] = 1
d2018 <- ts(c(d2018_1,d2018_2),  start = start(d2018_2),frequency = frequency(121))
d2018 <- ts(d2018, start = c(2005, 01), end = c(2019, 12), frequency = 12)


dlog.emae.reg <- lm(emae_num ~ -1 + d2009 + d2012 + d2018)
dlog.emae.adj <- dlog.emae.reg$residuals

dlog_residuos <- as.numeric(dlog.emae.adj)
dlog_residuos <- ts(dlog_residuos, start = c(2005, 01), end = c(2019, 12), frequency = 12)
dlog_residuos

Yd_6 <- cbind(Yd_6[,1:4], dlog_residuos) 
Yd_6_s_pcom <- Yd_6[, 2:5]


#Parte A: sin precios internacionales

popt <- VARselect(Yd_6_s_pcom, lag.max = pmax, type = "const")
p <- popt$selection[1] # AIC
p_6_s_pcom <- p

Yd0 <- Yd_6_s_pcom [1:pmax, ] # Initial values
Ydt <- Yd_6_s_pcom [(pmax - p + 1):nrow(Yd_6_s_pcom), ]

# Estimation
VAR <- VAR(Ydt, p = p, type = "const")

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs

m_6_s_pcom <- m
N_6_s_pcom <- N

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation (AB model configuration)
SVAR_6_s_pcom<- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

#Resultados parte A

Yb <- boot.rb.replicate(VAR, Yd0,  pmax, R)
SVARS(SVAR_6_s_pcom, 2, 1)
graficos(SVAR_6_s_pcom, 800, 800)

# IRF
#SVAR_6_s_pcom.SIRF <- SVAR.sirf(SVAR, H)
plot.sirf.boot(SVAR_6_s_pcom.SIRF.boot, m, H)

# Cumulative IRF
#SVAR_6_s_pcom.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
plot.sirf.boot(SVAR_6_s_pcom.SIRF.c.boot, m, H)

# FEVD
#SVAR_6_s_pcom.FEVD <- SVAR.fevd(SVAR, H)
plot.fevd.boot(SVAR_6_s_pcom.FEVD.boot, m, H)

# ERPT
#SVAR_6_s_pcom.ERPT <- SVAR.erpt(SVAR, H_ERPT, 2, 1, cumulative = TRUE)
plot.erpt.boot(SVAR_6_s_pcom.ERPT.boot, H_ERPT)

#Parte B: con precios internacionales

popt <- VARselect(Yd_6, lag.max = pmax, type = "const")

p <- popt$selection[1] # AIC
p_6_c_pcom<- p

Yd0 <- Yd_6[1:pmax, ] # Initial values
#MARTIN: No entiendo muy bien por qué, pero solo eliminamos 10 valores, no 12. 
#MARTIN: ¿No deberíamos comernos únicamente los primeros 3 valores?
Ydt <- Yd_6[(pmax - p + 1):nrow(Yd_6), ] # Starting in Jan-04

# Estimation
VAR <- VAR(Ydt, p = p, type = "const")

m_6_c_pcom <- VAR$K # No. of variables in the VAR
N_6_c_pcom <- VAR$obs # No. of effective sample observations, excluding "p" starting values
m <- m_6_c_pcom
N <- N_6_c_pcom 


Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# Ad hoc function
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}

# Re-estimate VAR (no feedback from local vars. to pcom)
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
VAR


# Model checking
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation (AB model configuration)
#aca simplemente estimamos el SVAR imponiendo ciertas assumptions embebidas en A y en B
#los NA en A y en B se van a estimar como en la slide 24 de la lecture 4
SVAR_6_c_pcom <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR


Yb <- boot.rb.replicate(VAR, Yd0,  pmax, R)
SVARS(SVAR_6_c_pcom, 3, 2)
graficos(SVAR_6_c_pcom, 800, 800)

#Reporte de resultados inciso 4 parte a2
#(quiza podriamos hacer graficos mas lindos con ggplot)
# 
# # IRF
# #SVAR_6_c_pcom.SIRF <- SVAR.sirf(SVAR, H)
# plot.sirf.boot(SVAR_6_c_pcom.SIRF.boot, m, H)
# 
# 
# 
# # Cumulative IRF
# #SVAR_6_c_pcom.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
# plot.sirf.boot(SVAR_6_c_pcom.SIRF.c.boot, m, H)
# png("my_plot.png")
# 
# # Code
# plot.sirf.boot(SVAR_6_c_pcom.SIRF.c.boot, m, H)
# 
# # Close device
# dev.off()
# # FEVD
# #SVAR_6_c_pcom.FEVD <- SVAR.fevd(SVAR, H)
# plot.fevd.boot(SVAR_6_c_pcom.FEVD.boot, m, H)
# 
# # HD
# #SVAR_6_c_pcom.HD <- SVAR.hd(SVAR)
# #plot.hd.boot(Yd, SVAR_6_c_pcom.HD.boot, m, pmax)
# 
# # ERPT
# #SVAR_6_c_pcom.ERPT <- SVAR.erpt(SVAR, H_ERPT, 3, 2, cumulative = TRUE)
# plot.erpt.boot(SVAR_6_c_pcom.ERPT.boot, H_ERPT)

# Punto 7 ####

remove(list = ls(all.names = TRUE))
gc()
library(vars)
library(seasonal)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

M2private <- read_excel("m2.xlsx")
M2private <- ts(M2private, start = c(2003, 12), frequency = 12)
seas.adj <- seas(M2private)
M2adj <- seas.adj$series$s11
source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc, M2adj)
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences
Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

pmax <- 12
popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]
p

Yd0 <- Yd[1:pmax, ] # Initial values, 12M
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # [check as it depends upon the lags]

VAR <- VAR(Ydt, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM
pars.S

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# Punto 8 ####
remove(list = ls(all.names = TRUE))
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pf <- read_excel("pftna.xlsx")
pf <- ts(pf, start = c(2003, 01), frequency = 12)
source("PS1_Data.R")


Yl.f <- cbind(pcom, er, pc)
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl.f <- cbind(Yl.f, pf)
pf.d <- diff(pf)
Yd.f <- cbind(Yd.f, pf.d)

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

pmax <- 12
popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # [check and change for each sample} /// [change again if autocorrelated errors]
p

Yd0 <- Yd[1:pmax, ] # Initial values, 12M
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # [check as it depends upon the lags]


VAR <- VAR(Ydt, p = p, type = "const")
m <- VAR$K
N <- VAR$obs
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))

roots(VAR, modulus = TRUE) #model check
serial.test(VAR, lags.bg = 12, type = "ES") #model check

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)

P <- solve(SVAR$A, SVAR$B)
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) #SVAR #check no feedback PCOM
pars.S

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon IR
H_ERPT <- 120 # Horizon ERPT - variance decomp.
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# Punto 9 ####

#TBC

# Punto 10 ####
#Forma 1: fuente oficial

#Descargamos los valores de la fuente, tomando el promedio mensual:
dolar_ccl <- read.csv(url("https://apis.datos.gob.ar/series/api/series/?collapse=month&collapse_aggregation=avg&ids=168.1_T_CAMBIDRS_D_0_0_29&limit=5000&format=csv"))

#Construimos la serie y la restringimos
dolar_ccl <- ts(dolar_ccl$tipo_cambio_implicito_en_adrs, start = c(2002, 04), frequency = 12)
dolar_ccl <- window(dolar_ccl, start = c(2004, 01), end = c(2019, 12))

#Graficamos para verificar
plot(dolar_ccl, type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="CCL implícito en ADRs", ylim=c(0,80))


#Forma 2: tomando valores de la acción del Banco Galicia
#Para obtener tipo de cambio implícito en los ADR (Contado Contra Liquidación)
#es necesario hacer la siguiente cuenta:
#Dólar CCL = (Precio Local de Acción/Precio del ADR)* Factor de conversión
#En el caso del Grupo Financiero Galicia S.A. (Banco Galicia), el factor correspondiente es 10.


#Utilizamos una descarga automática basada en Yahoo Finance de los valores en NY y en BA, tomando precio de cierre mensual ajustado:
library(tseries)
P_adr <- get.hist.quote(instrument = "GGAL", start = "2004-01-01", end = "2019-12-31", 
                        quote = "AdjClose", provider = "yahoo", compression = "m", 
                        retclass = "zoo")
plot(P_adr,type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="GGAL Adjusted Price", ylim=c(0,70))

P_local <- get.hist.quote(instrument = "GGAL.BA", start = "2004-01-01", end = "2019-12-31", 
                          quote = "AdjClose", provider = "yahoo", compression = "m", 
                          retclass = "zoo")
plot(P_local,type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="GGAL.BA Adjusted Price", ylim=c(0,170))

#Construimos la serie
dolar_implicito <- (P_local/P_adr)*10
dolar_ccl2 = ts(dolar_implicito, start = c(2004, 01), frequency = 12)

#Graficamos para verficar
plot(dolar_ccl2, type = "l",lwd=2, col="black", xlab="",ylab="",bty="n", main="CCL a partir de $GGAL", ylim=c(0,140))


#Analizamos y graficamos las dos alternativas para comparar:
discrepancia = (dolar_ccl-dolar_ccl2)/dolar_ccl
mean(discrepancia)
sd(discrepancia)
ts.plot(discrepancia)
#Conclusión: aunque en promedio no son muy distintas, las medidas difieren bastante en algunos meses. Parece ser mejor usar el oficial, ante la duda.



# Comenzamos el análisis VAR
# Definimos las nuevas variables, con ventana temporal enero2005-diciembre2019
Yl.f_10 <- cbind(pcom, dolar_ccl, pc) # Variables en log
Yl.f_10 <- log(Yl.f_10) # log transformation
Yd.f_10 <- 100 * diff(Yl.f_10) # Variables en log-differences
Yl_10 <- window(Yl.f_10, start = c(2005, 01), end = c(2019, 12))
Yd_10 <- window(Yd.f_10, start = c(2005, 01), end = c(2019, 12))


# Volvemos a llamar al paquete, para evitar <<enmascaramiento>> del paquete ts:
library(vars)

# Lag order selection

SVAR_10popt_10 <- VARselect(Yd_10, lag.max = pmax, type = "const")
popt_10
p_10 <- popt_10$selection[1] # AIC

# Valores iniciales
Yd0_10 <- Yd_10[1:pmax, ] # Initial values
Ydt_10 <- Yd_10[(pmax - p_10 + 1):nrow(Yd_10), ] 

# Estimation
VAR_10 <- VAR(Ydt_10, p = p_10, type = "const")

# Control
m_10 <- VAR_10$K # No. of variables in the VAR
N_10 <- VAR_10$obs
roots(VAR_10, modulus = TRUE)
serial.test(VAR_10, lags.bg = 2, type = "ES") #DUDA: me lo estima con un solo lag, ¿les parece OK hacerlo así?

# Re-estimación con restricciones:
# Re-estimate VAR (no feedback from local vars. to pcom)
VAR_10 <- restrict(VAR_10, method = "man", resmat = matC(m_10, p_10, 1))
VAR_10

# SVAR estimation

# A Matrix
Amat_10 <- diag(m_10)
for (i in 2:m_10) {
  for (j in 1:(i - 1)) {
    Amat_10[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat_10 <- matrix(0, m_10, m_10)
for (i in 1:m_10) {
  Bmat_10[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_10 <- SVAR(VAR_10, Amat = Amat_10, Bmat = Bmat_10, lrtest = FALSE)

# Reportamos los resultados del modelo
m <- m_10
# IRF
SVAR_10.SIRF <- SVAR.sirf(SVAR_10, H)
plot.sirf(SVAR_10.SIRF, m=m_10, H)

# Cumulative IRF
SVAR_10.SIRF.c <- SVAR.sirf(SVAR_10, H, cumulative = TRUE)
plot.sirf(SVAR_10.SIRF.c, m = m_10, H)

# FEVD
SVAR_10.FEVD <- SVAR.fevd(SVAR_10, H)
plot.fevd(SVAR_10.FEVD, m = m_10, H)

# HD 
N <- N_10 #Hay que redifinir para que funcione la función SVAR.hd
SVAR_10.HD <- SVAR.hd(SVAR_10)
plot.hd(Yd_10, SVAR_10.HD, m = m_10, pmax)

# ERPT
SVAR_10.ERPT <- SVAR.erpt(SVAR_10, H_ERPT, 3, 2, cumulative = TRUE)
plot.erpt(SVAR_10.ERPT, H_ERPT)

# Replicación con bootstrap y gráfico final con bandas de confianza
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications
Yb_10 <- boot.rb.replicate(VAR_10, Yd0_10, pmax, R)
SVAR.SIRF.boot_10 <- SVAR.sirf.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_10, m = m_10, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_10 <- SVAR.sirf.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_10, m = m_10 , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_10 <- SVAR.fevd.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_10, m = m_10, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_10 <- SVAR.erpt.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot_10, H_ERPT)

# Punto 11: primer ordenamiento ####
# Definimos las nuevas variables, con ventana temporal enero2005-diciembre2019
pcom_log <- log(pcom)
er_log <- log(er)
brecha_log <- log(dolar_ccl)-log(er)
pc_log <- log(pc)

Yl.f_11 <- cbind(er_log, brecha_log, pc_log) 
Yd.f_11 <- 100 * diff(Yl.f_11) # Variables en log-differences
Yl_11 <- window(Yl.f_11, start = c(2005, 01), end = c(2019, 12))
Yd_11 <- window(Yd.f_11, start = c(2005, 01), end = c(2019, 12))

# Volvemos a llamar al paquete, para evitar <<enmascaramiento>> del paquete ts:
library(vars)

# Lag order selection
popt_11 <- VARselect(Yd_11, lag.max = pmax, type = "const")
popt_11
p_11 <- popt_11$selection[1] # AIC

# Valores iniciales
Yd0_11 <- Yd_11[1:pmax, ] # Initial values
Ydt_11 <- Yd_11[(pmax - p_11 + 1):nrow(Yd_11), ] 

# Estimation
VAR_11 <- VAR(Ydt_11, p = p_11, type = "const")

# Control
m_11 <- VAR_11$K # No. of variables in the VAR
N_11 <- VAR_11$obs
roots(VAR_11, modulus = TRUE)
serial.test(VAR_11, lags.bg = 3, type = "ES")

# SVAR estimation

# A Matrix
Amat_11 <- diag(m_11)
for (i in 2:m_11) {
  for (j in 1:(i - 1)) {
    Amat_11[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat_11 <- matrix(0, m_11, m_11)
for (i in 1:m_11) {
  Bmat_11[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_11 <- SVAR(VAR_11, Amat = Amat_11, Bmat = Bmat_11, lrtest = FALSE)

# Reportamos los resultados del modelo
m <- m_11
N <- N_11 #Hay que redifinir para que funcione la función SVAR.hd
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications

Yb_11 <- boot.rb.replicate(VAR_11, Yd0_11, pmax, R)
SVAR.SIRF.boot_11 <- SVAR.sirf.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_11, m = m_11, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_11 <- SVAR.sirf.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_11, m = m_11 , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_11 <- SVAR.fevd.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_11, m = m_11, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_11 <- SVAR.erpt.boot(SVAR_11, Amat_11, Bmat_11, Yb_11, pmax, H_ERPT, 3, 1, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot_11, H_ERPT)


#Punto 11: segundo ordenamiento ####

#Nuevo orden:
Yl.f_11b <- cbind(pcom_log, brecha_log, er_log, pc_log) 
Yd.f_11b <- 100 * diff(Yl.f_11b) # Variables en log-differences
Yl_11b <- window(Yl.f_11b, start = c(2005, 01), end = c(2019, 12))
Yd_11b <- window(Yd.f_11b, start = c(2005, 01), end = c(2019, 12))

# Comenzamos análisis VAR
popt_11b <- VARselect(Yd_11b, lag.max = pmax, type = "const")
popt_11b
p_11b <- popt_11b$selection[1] # AIC

# Valores iniciales
Yd0_11b <- Yd_11[1:pmax, ] # Initial values
Ydt_11b <- Yd_11b[(pmax - p_11b + 1):nrow(Yd_11b), ] 

# Estimation
VAR_11b <- VAR(Ydt_11b, p = p_11b, type = "const")

# Control
m_11b <- VAR_11b$K # No. of variables in the VAR
N_11b <- VAR_11b$obs
roots(VAR_11b, modulus = TRUE)
serial.test(VAR_11b, lags.bg = 1, type = "ES") #DUDA: me lo estima con un solo lag, ¿les parece OK hacerlo así?

# Re-estimación con restricciones:
# Re-estimate VAR (no feedback from local vars. to pcom)
VAR_11b <- restrict(VAR_11b, method = "man", resmat = matC(m_11b, p_11b, 1))
VAR_11b

# SVAR estimation

# A Matrix
Amat_11b <- diag(m_11b)
for (i in 2:m_11b) {
  for (j in 1:(i - 1)) {
    Amat_11b[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat_11b <- matrix(0, m_11b, m_11b)
for (i in 1:m_11b) {
  Bmat_11b[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_11b <- SVAR(VAR_11b, Amat = Amat_11b, Bmat = Bmat_11b, lrtest = FALSE)


# Replicación con bootstrap y gráfico final con bandas de confianza
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications

Yb_11b <- boot.rb.replicate(VAR_11b, Yd0_11b, pmax, R)
N <- N_11b
m <- m_11b
SVAR.SIRF.boot_11b <- SVAR.sirf.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_11b, m = m_11b, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_11b <- SVAR.sirf.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_11b, m = m_11b , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_11b <- SVAR.fevd.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_11b, m = m_11b, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_11b <- SVAR.erpt.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H_ERPT, 4, 3, a, R, cumulative = TRUE) 
plot.erpt.boot(SVAR.ERPT.boot_11b, H_ERPT)



#Punto 12: primer ordenamiento ####
#Tomamos como períodos con controles de capitales al período octubre 2011-diciembre2015 y desde septiembre 2019.
#Fuentes: AFIP con la Resolución General 3210 y la Resolución General 3819 , Poder Ejecutivo Nacional con DNU 19/609
@@ -1673,3 +1749,77 @@ plot.fevd.boot(SVAR.FEVD.boot_12, m = m_12, H)
SVAR.ERPT.boot_12 <- SVAR.erpt.boot(SVAR_12, Amat_12, Bmat_12, Yb_12, pmax, H_ERPT, 4, 2, a, R, cumulative = TRUE) # DUDA: Acá puse el 4 en vez del 3 porque sino había un problema con las dimensiones y no estimaba, pero NO ESTOY muy seguro. No me termina de quedar claro qué representan estos argumentos (los números) en la función.
plot.erpt.boot(SVAR.ERPT.boot_12, H_ERPT)

# Punto 12: segundo ordenamiento ####

#Armamos vectores de variables. Primero definimos:
library(tstools) 
dummy_cepo1 <- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(2011,10), dummy_end =c(2015,12), sp= NULL, start_basic = c(2004, 01), frequency = 12)
dummy_cepo2 <- create_dummy_ts(end_basic = c(2019,12), dummy_start = c(2019,09), dummy_end =c(2019,12), sp= NULL, start_basic = c(2004, 01), frequency = 12)

#Creamos la variable de la brecha que toma valor 0 cuando no hay controles de capitales:
brecha_con_cepo1 <- brecha_log*dummy_cepo1
brecha_con_cepo1 <- window(brecha_con_cepo1, end = c(2019, 08))
brecha_con_cepo2 <- brecha_log*dummy_cepo2
brecha_con_cepo2 <- window(brecha_con_cepo2, start = c(2019, 09))
brecha_con_cepo_log <- concat_ts(brecha_con_cepo1, brecha_con_cepo2)
#Eliminamos variables intermedias:
remove(brecha_con_cepo1, brecha_con_cepo2)
#Volvemos a llamar al paquete, para evitar <<enmascaramientos>> con paquete ts
library(vars)


#Ahora sí, vectores
Yl.f_12b <- cbind(pcom_log, brecha_con_cepo_log, er_log, pc_log) 
Yd.f_12b <- 100 * diff(Yl.f_12b) # Variables en log-differences
Yl_12b <- window(Yl.f_12b, start = c(2005, 01), end = c(2019, 12))
Yd_12b <- window(Yd.f_12b, start = c(2005, 01), end = c(2019, 12))

# Comenzamos análisis VAR
popt_12b <- VARselect(Yd_12b, lag.max = pmax, type = "const")
popt_12b
p_12b <- popt_12b$selection[1] # AIC

# Valores iniciales
Yd0_12b <- Yd_12b[1:pmax, ] # Initial values
Ydt_12b <- Yd_12b[(pmax - p_12b + 1):nrow(Yd_12b), ] 

# Estimation
VAR_12b <- VAR(Ydt_12b, p = p_12b, type = "const")

# Control
m_12b <- VAR_12b$K # No. of variables in the VAR
N_12b <- VAR_12b$obs
roots(VAR_12b, modulus = TRUE)
serial.test(VAR_12b, lags.bg = 1, type = "ES") #DUDA: me lo estima con un solo lag, ¿les parece OK hacerlo así?

# Re-estimación con restricciones:
# Re-estimate VAR (no feedback from local vars. to pcom)
VAR_12b <- restrict(VAR_12b, method = "man", resmat = matC(m_12b, p_12b, 1))
VAR_12b

# SVAR estimation

# A Matrix
Amat_12b <- diag(m_12b)
for (i in 2:m_12b) {
  for (j in 1:(i - 1)) {
    Amat_12b[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat_12b <- matrix(0, m_12b, m_12b)
for (i in 1:m_12b) {
  Bmat_12b[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_12b <- SVAR(VAR_12b, Amat = Amat_12b, Bmat = Bmat_12b, lrtest = FALSE)


# Replicación con bootstrap y gráfico final con bandas de confianza
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications

Yb_12b <- boot.rb.replicate(VAR_12b, Yd0_12b, pmax, R)
N <- N_12b
m <- m_12b
SVAR.SIRF.boot_12b <- SVAR.sirf.boot(SVAR_12b, Amat_12b, Bmat_12b, Yb_12b, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_12b, m = m_12b, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_12b <- SVAR.sirf.boot(SVAR_12b, Amat_12b, Bmat_12b, Yb_12b, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_12b, m = m_12b , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_12b <- SVAR.fevd.boot(SVAR_12b, Amat_12b, Bmat_12b, Yb_12b, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_12b, m = m_12b, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_12b <- SVAR.erpt.boot(SVAR_12b, Amat_12b, Bmat_12b, Yb_12b, pmax, H_ERPT, 4, 3, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot_12b, H_ERPT)
