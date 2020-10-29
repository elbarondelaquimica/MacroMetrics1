# Conversión para ggplot2

#Paquetes necesarios
library(tidyverse)
library(reshape2)
library("tibble")
library("gridExtra")
library("dplyr")
library("Lock5Data")
library("ggthemes")
library("fun")
library("zoo")
library("corrplot")
library("maps")
library("mapproj")

#Genero la función

grafico <- function(SVAR.ERPT.boot, H_ERPT){
  t = c(0:120)
  datos <- data.frame(t, SVAR.ERPT.boot$lb, SVAR.ERPT.boot$pe, SVAR.ERPT.boot$ub)
  colnames(datos) <- c("t", "lb", "pe", "ub")
  grafico <- ggplot(data=datos,aes(x= t, y=pe))
  grafico <- grafico + geom_ribbon(aes(x=t, ymax=ub, ymin=lb), fill="olivedrab3", alpha=.25) 
  grafico <-  grafico +geom_line(aes(y = ub), colour = 'darkolivegreen3') 
  grafico <-  grafico + geom_line(aes(y = lb), colour = 'darkolivegreen3')
  grafico <- grafico + geom_line(colour = "darkolivegreen4", size =1.15)
  grafico <- grafico + labs(x = "Mes", y =" ", title = "Exchange rate pass through", subtitle = "En porcentaje")
  grafico <- grafico + theme_minimal()
  grafico <- grafico + theme(axis.ticks = element_blank())
  grafico
}

