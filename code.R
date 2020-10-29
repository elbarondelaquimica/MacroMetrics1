remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")


#Reporte de resultados inciso 1
#(quiza podriamos hacer graficos mas lindos con ggplot)

# IRF
SVAR.SIRF <- SVAR.sirf(SVAR, H)
plot.sirf(SVAR.SIRF, m, H)

# Cumulative IRF
SVAR.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
plot.sirf(SVAR.SIRF.c, m, H)

# FEVD
SVAR.FEVD <- SVAR.fevd(SVAR, H)
plot.fevd(SVAR.FEVD, m, H)

# HD
SVAR.HD <- SVAR.hd(SVAR)
plot.hd(Yd, SVAR.HD, m, pmax)

# ERPT
SVAR.ERPT <- SVAR.erpt(SVAR, H_ERPT, 3, 2, cumulative = TRUE)
plot.erpt(SVAR.ERPT, H_ERPT)

#PENDIENTE: Revisar cuestiones de estacionariedad de las series.


# Inciso 2 ####
#Estimamos el modelo sin el índice de precios internos.

#Duda: ¿hay que hacer esto de nuevo o con los resultados que obtuvimos con las 3 variables alcanza?
Yd_2 <- Yd[, 2:3]
popt_2 <- VARselect(Yd_2, lag.max = pmax, type = "const")
popt_2
p_2 <- popt_2$selection[1] # AIC

Yd0_2 <- Yd_2[1:pmax, ] # Initial values
Ydt_2 <- Yd_2[(pmax - p_2 + 1):nrow(Yd_2), ] # Starting in Jan-04

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

# IRF
SVAR_2.SIRF <- SVAR.sirf(SVAR_2, H)
plot.sirf(SVAR_2.SIRF, m, H)

# Cumulative IRF
SVAR_2.SIRF.c <- SVAR.sirf(SVAR_2, H, cumulative = TRUE)
plot.sirf(SVAR_2.SIRF.c, m, H)

# FEVD
SVAR_2.FEVD <- SVAR.fevd(SVAR_2, H)
plot.fevd(SVAR_2.FEVD, m, H)

# HD
SVAR_2.HD <- SVAR.hd(SVAR_2)
plot.hd(Yd_2, SVAR_2.HD, m, pmax)

# ERPT
SVAR_2.ERPT <- SVAR.erpt(SVAR_2, H_ERPT, 2, 1, cumulative = TRUE)
plot.erpt(SVAR_2.ERPT, H_ERPT)

#Pendiente: comparar con resultados modelo inciso 1.








# Punto 10 ####
# Definimos las nuevas variables, con ventana temporal enero2005-diciembre2019
Yl.f_10 <- cbind(dolar_ccl, pc) # Variables en log
Yl.f_10 <- log(Yl.f_10) # log transformation
Yd.f_10 <- 100 * diff(Yl.f_10) # Variables en log-differences
Yl_10 <- window(Yl.f_10, start = c(2005, 01), end = c(2019, 12))
Yd_10 <- window(Yd.f_10, start = c(2005, 01), end = c(2019, 12))

# Volvemos a llamar al paquete, para evitar <<enmascaramiento>> del paquete ts:
library(vars)

# Lag order selection
popt_10 <- VARselect(Yd_10, lag.max = pmax, type = "const")
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
serial.test(VAR_10, lags.bg = 2, type = "ES")

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
N <- N_10 #Hay que redifinir para que funcione la función SVAR.hd
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
SVAR.ERPT.boot_10 <- SVAR.erpt.boot(SVAR_10, Amat_10, Bmat_10, Yb_10, pmax, H_ERPT, 2, 1, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot_10, H_ERPT)


# Punto 11: primer ordenamiento ####
# Definimos las nuevas variables, con ventana temporal enero2005-diciembre2019
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


# Punto 11: segundo ordenamiento ####
# Definimos las nuevas variables, con ventana temporal enero2005-diciembre2019
Yl.f_11b <- cbind(brecha_log, er_log, pc_log) 
Yd.f_11b <- 100 * diff(Yl.f_11b) # Variables en log-differences
Yl_11b <- window(Yl.f_11b, start = c(2005, 01), end = c(2019, 12))
Yd_11b <- window(Yd.f_11b, start = c(2005, 01), end = c(2019, 12))

# Volvemos a llamar al paquete, para evitar <<enmascaramiento>> del paquete ts:
library(vars)

# Lag order selection
popt_11b <- VARselect(Yd_11b, lag.max = pmax, type = "const")
popt_11b
p_11b <- popt_11b$selection[1] # AIC

# Valores iniciales
Yd0_11b <- Yd_11b[1:pmax, ] # Initial values
Ydt_11b <- Yd_11b[(pmax - p_11b + 1):nrow(Yd_11b), ] 

# Estimation
VAR_11b <- VAR(Ydt_11b, p = p_11b, type = "const")

# Control
m_11b <- VAR_11b$K # No. of variables in the VAR
N_11b <- VAR_11b$obs
roots(VAR_11b, modulus = TRUE)
serial.test(VAR_11b, lags.bg = 3, type = "ES")

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

# Reportamos los resultados del modelo
m <- m_11b
N <- N_11b #Hay que redifinir para que funcione la función SVAR.hd
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications

Yb_11b <- boot.rb.replicate(VAR_11b, Yd0_11b, pmax, R)
SVAR.SIRF.boot_11b <- SVAR.sirf.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_11b, m = m_11b, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_11b <- SVAR.sirf.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_11b, m = m_11b , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_11b <- SVAR.fevd.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_11b, m = m_11b, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_11b <- SVAR.erpt.boot(SVAR_11b, Amat_11b, Bmat_11b, Yb_11b, pmax, H_ERPT, 3, 1, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot_11, H_ERPT)


#Punto 12: primer ordenamiento ####
#Tomamos como períodos con controles de capitales al período octubre 2011-diciembre2015 y desde septiembre 2019.
#Fuentes: AFIP con la Resolución General 3210 y la Resolución General 3819 , Poder Ejecutivo Nacional con DNU 19/609

#Variables
Yl.f_12 <- cbind(er_log, brecha_con_cepo_log, pc_log) 
Yd.f_12 <- 100 * diff(Yl.f_12) # Variables en log-differences
Yl_12 <- window(Yl.f_12, start = c(2005, 01), end = c(2019, 12))
Yd_12 <- window(Yd.f_12, start = c(2005, 01), end = c(2019, 12))

# Comenzamos análisis VAR
popt_12 <- VARselect(Yd_12, lag.max = pmax, type = "const")
popt_12
p_12 <- popt_12$selection[1] # AIC

# Valores iniciales
Yd0_12 <- Yd_12[1:pmax, ] # Initial values
Ydt_12 <- Yd_12[(pmax - p_12 + 1):nrow(Yd_12), ] 

# Estimation
VAR_12 <- VAR(Ydt_12, p = p_12, type = "const")

# Control
m_12 <- VAR_12$K # No. of variables in the VAR
N_12 <- VAR_12$obs
roots(VAR_12, modulus = TRUE)
serial.test(VAR_12, lags.bg = 3, type = "ES")

# SVAR estimation

# A Matrix
Amat_12 <- diag(m_12)
for (i in 2:m_12) {
  for (j in 1:(i - 1)) {
    Amat_12[i, j] <- NA
  }
}

# B Matrix
#esta representa a omega (en el caso en el que no normalizamos el desvio de los errores a 1)
Bmat_12 <- matrix(0, m_12, m_12)
for (i in 1:m_12) {
  Bmat_12[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR_12 <- SVAR(VAR_12, Amat = Amat_12, Bmat = Bmat_12, lrtest = FALSE)


# Replicación con bootstrap y gráfico final con bandas de confianza
a <- 0.95 # Confidence level
R <- 1000 # No. of bootstrap replications

Yb_12 <- boot.rb.replicate(VAR_12, Yd0_12, pmax, R)
N <- N_12
m <- m_12
SVAR.SIRF.boot_12 <- SVAR.sirf.boot(SVAR_12, Amat_12, Bmat_12, Yb_12, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot_12, m = m_12, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot_12 <- SVAR.sirf.boot(SVAR_12, Amat_12, Bmat_12, Yb_12, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_12, m = m_12 , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_12 <- SVAR.fevd.boot(SVAR_12, Amat_12, Bmat_12, Yb_12, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_12, m = m_12, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_12 <- SVAR.erpt.boot(SVAR_12, Amat_12, Bmat_12, Yb_12, pmax, H_ERPT, 3, 1, a, R, cumulative = TRUE) # DUDA: Acá puse el 4 en vez del 3 porque sino había un problema con las dimensiones y no estimaba, pero NO ESTOY muy seguro. No me termina de quedar claro qué representan estos argumentos (los números) en la función.
plot.erpt.boot(SVAR.ERPT.boot_12, H_ERPT)

# Punto 12: segundo ordenamiento ####

#Armamos vectores de variables
#Variables
Yl.f_12b <- cbind(brecha_con_cepo_log, er_log, pc_log) 
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
serial.test(VAR_12b, lags.bg = 3, type = "ES")

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
for (i in 1:m_12) {
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
SVAR.SIRF.c.boot_12b <- SVAR.sirf.boot(SVAR_12b, Amat_12b, Bmat_12, Yb_12, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot_12b, m = m_12b , H)

# FEVD (bootstrap)
SVAR.FEVD.boot_12b <- SVAR.fevd.boot(SVAR_12b, Amat_12b, Bmat_12b, Yb_12b, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot_12b, m = m_12b, H)

# ERPT (bootstrap)
SVAR.ERPT.boot_12b <- SVAR.erpt.boot(SVAR_12b, Amat_12b, Bmat_12b, Yb_12b, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE) # DUDA: Acá puse el 4 en vez del 3 porque sino había un problema con las dimensiones y no estimaba, pero NO ESTOY muy seguro. No me termina de quedar claro qué representan estos argumentos (los números) en la función.
plot.erpt.boot(SVAR.ERPT.boot_12b, H_ERPT)
