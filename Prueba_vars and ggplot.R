#Prueba VAR and ggplot

#Función para convertir a data frame lo que sale de la función ERPT:
convertir_a_df <- function(SVAR.erpt.boot, H_ERPT){
  t = c(0:H_ERPT)
  data.frame(t, SVAR.erpt.boot$lb, SVAR.erpt.boot$pe, SVAR.erpt.boot$ub)
}

SVAR.erpt.boot$lb <- as.numeric(SVAR.erpt.boot$lb)
#Acá tomo uno estimado para ver si funciona
prueba <- convertir_a_df(SVAR.ERPT.boot_12, H_ERPT = 120)
prueba

#Parece que funciona, así que ahora hay que graficar. Antes, hay que cambiar el formato del dataframe:
# TBC

#Ahora sí, el dichoso gráfico:
library(tidyverse)
grafico1 <- ggplot(data=prueba, aes_string(prueba$t, prueba$lb, prueba$pe, prueba$up)) +
  geom_line(aes(y = prueba$lb), colour = 'lightblue2')