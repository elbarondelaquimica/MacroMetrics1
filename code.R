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

