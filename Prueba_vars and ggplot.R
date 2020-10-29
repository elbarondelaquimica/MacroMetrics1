# Conversión para ggplot2

#Paquetes necesarios
library(tidyverse)
library(reshape2)

#Defino la función: hay que especificar el nombre con el que guardamos a la función de ERPT, el horizonte a graficar y el nombre (aunque esto está medio al cohete)
convertir_a_df <- function(SVAR.erpt.boot, H_ERPT){
  t = c(0:H_ERPT)
  df_aux <- data.frame(t, SVAR.erpt.boot$lb, SVAR.erpt.boot$pe, SVAR.erpt.boot$ub)
  melt(df_aux,id="t")
}

#Aplico la función para graficar un caso.
#tomo el punto 12 de ejemplo. El gráfico es el más sencillo, hay que agregarle muchas cosas, a gusto del consumidor (?)
#prueba <- convertir_a_df(SVAR.ERPT.boot_12, H_ERPT = 120)
#grafico <- ggplot(data=prueba,aes(x= t, y=value,colour=variable,group=variable)) +
#  geom_line()
#grafico
