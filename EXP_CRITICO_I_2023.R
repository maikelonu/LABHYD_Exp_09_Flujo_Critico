# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: FLUJO NO-UNIFORME/ FLUJO CRITICO

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////

# INFO:
# Analisis grafico avanzado
# ggplot2
# lattice
# Normalizacion y homogenizacion de variables
# Exportaci?n ASCII"
# //////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
rm(list = ls())

# Working directory is selected
setwd("/media/maikel/Trabajo/R_ITC/R_LABHYD/EXP_CRITICO")

# CRAN libraries are loaded
# require(Agreement)
require(DescTools)
require(effects)
require(ggplot2)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(reshape)
require(visreg)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -3)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# ////////////////////////////////////////////////////////
# BLOCK: Declarations
# ////////////////////////////////////////////////////////
base_m <- 0.086 # hydraulic flume base (m)
     q <- 4.00 # water flow (m3/h)
     q <- q/3600 # water flow (m3/s)

     # Triangular Weir dimensions and vectors are defined
block_ID <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9") # control-points IDs
block_X <- c(15.050, 15.505, 15.550, 15.602, 15.670, 15.740, 15.810, 15.870, 16.120) # control-points X distance (m)
block_DZ <- c(0.45, 0.45, 0.45, 0.45, 0.50, 0.55, 0.55, 0.55, 0.45) # hydraulic flume DeltaZ to bottom (cm)
block_ELEV <- c(0.0, 0.0, 2.7, 6.2, 4.6, 3.0, 1.5, 0.0, 0.0) # relative elevation of the Triangular Weir (cm)

# Experimental Data !!!!!!!!!!!!!!
block_Y <- c(10.40, 10.35, 10.25, 9.65, 7.10, 5.10, 3.50, 2.15, 1.90) # experimental water depth (cm)

# Vectors are transformed
    block_X <- (block_X) - 15.05 # normalized to start in 0 (m)
   block_DZ <- block_DZ / 100 # transformed to (m)
 block_ELEV <- block_ELEV / 100 # transformed to (m)
    block_Y <- block_Y / 100 # transformed to (m)
 block_YEFF <- block_Y - block_DZ # effective water depth is calculated for plotting (m) 
block_YEFF2 <- block_YEFF - block_ELEV # effective water depth is calculated for numerica analysis (m)

# A data.frame object is created   
df.weir <- data.frame(block_ID,block_X,block_DZ,block_ELEV,block_Y,block_YEFF,block_YEFF2)

# Static energy is calculated (m)
df.weir$block_energ_EST <- df.weir$block_YEFF2

# Dynamic energy is calculated (m)
df.weir$block_energ_DYM <- (q^2) / (2*9.81*(base_m^2)*((df.weir$block_energ_EST)^2))

# Total energy is calculated (m)
df.weir$block_energ_TOTAL <- df.weir$block_energ_EST + df.weir$block_energ_DYM

# Theoretical total-energy is calculated (m) 
df.weir$block_energ_TOTAL_TEO <- df.weir$block_energ_TOTAL[1]

# Energy losses are calculated (m) 
df.weir$block_energ_LOSS <- df.weir$block_energ_TOTAL_TEO - df.weir$block_energ_TOTAL

# Energy losses are calculated for plotting(m) 
df.weir$block_energ_PLOT <- df.weir$block_energ_TOTAL + df.weir$block_ELEV

# hydraulic area is calculated (m2)
df.weir$area <- (df.weir$block_YEFF2) * base_m

# hydraulic perimeter is calculated (m)
df.weir$perimeter <- ((df.weir$block_YEFF2) * 2) + base_m

# water velocity is calculated (m/s)
df.weir$vel <- (q / df.weir$area)

# Froude number is calculated
df.weir$Froude <- df.weir$vel / ((df.weir$area * 9.81 / base_m) ^ 0.5)

# If-statement & For-loop for Froude number

# A rule variable is created
df.weir$rule <- NA

# A master counter is defined
counter <- length(df.weir$block_X)
for(i in 1:counter) {
  if (df.weir$Froude[i] > 1) {
    df.weir$rule[i] = "SUPER"
  }
  else {
    df.weir$rule[i] = "SUB"
  }
}

# approxfun {stats} is used to interpolate yc (m)
tfun01 <- approxfun(df.weir$Froude, df.weir$block_energ_EST)

# yc is estimated (m)
yc.inter <- round(tfun01(1),4) # Froude number = 1

# approxfun {stats} is used to interpolate yc position (m)
tfun02 <- approxfun(df.weir$Froude, df.weir$block_X)

# yc distance is estimated (m)
yc.inter.dist <- round(tfun02(1),4) # Froude number = 1

# approxfun {stats} is used to interpolate Ec (m)
tfun03 <- approxfun(df.weir$Froude, df.weir$block_energ_TOTAL)

# Ec is estimated (m)
Ec <- round(tfun03(1),4) # Froude number = 1

# yc is calculated based on ec.03 (m)
yc.ec03 <- round(((q^2) / ((base_m^2)*9.81)) ^ (1/3),4)

# yc is calculated based on ec.04 (m)
yc.ec04 <- (2/3)*Ec


cols <- c("pto_control"="firebrick4",
          "Energ_Est"="deepskyblue1",
          "base_canal" = "black",
          "vertedor" = "darkgreen",
          "Energ_Total_Teo" = "magenta3",
          "Energ_Total_Exp" = "lightsalmon3",
          "yc" = "#009999",
          "y_effec" = "#00cc00")

# A ggplot object is created
fg01 <- ggplot() +
  geom_path(aes(x = block_X,y = block_ELEV, colour = "vertedor"),data=df.weir,size = 1.50) +
  geom_hline(aes(yintercept = 0.0, colour = "base_canal"),data=df.weir,size = 1.50, linetype = 1) +
  geom_hline(aes(yintercept = df.weir$block_energ_TOTAL[1],colour = "Energ_Total_Teo"), data=df.weir,size = 0.75,linetype = 5) +
  geom_vline(aes(xintercept = block_X, color = "pto_control"),data=df.weir,size = 0.75,linetype = 2,alpha = 0.50) +
  geom_line(aes(x = block_X,y = block_YEFF,colour = "Energ_Est"),data=df.weir,size = 1.50) +
  geom_line(aes(x = block_X,y = block_energ_PLOT, colour = "Energ_Total_Exp"),data=df.weir,size = 1.00) +
  geom_vline(aes(xintercept = yc.inter.dist, colour = "yc"), data=df.weir,size = 1.25,linetype = 6) +
  geom_linerange(aes(x = block_X,ymin = block_ELEV,ymax = block_YEFF,colour = "y_effec"),data=df.weir,size = 0.50,linetype = 1) +
  scale_colour_manual(name="Numenclatura",values=cols) +  
  geom_text(aes(x = block_X,y = block_energ_TOTAL_TEO,label = rule,vjust = -1.5),data=df.weir,colour = '#666666',size = 3.5,hjust = 0.25,parse = FALSE) +
  geom_text(aes(x = block_X,y = block_ELEV,label = block_ID,vjust = 2.0),data=df.weir,colour = '#666666',size = 7.5,hjust = 0.25,parse = FALSE) +
  scale_y_continuous(limits = c(-0.025,0.15)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Perfil Energetico. Vertedor Triangular") +
  xlab("Distancia (m)") +
  ylab("Energia (m)") +
  theme_bw(base_size = 18.0) 

# A ggplot object is requested
fg01

# A ggplot object is created
fg02 <- ggplot() +
  geom_point(aes(x = block_energ_TOTAL,y = block_energ_EST),data=df.weir,shape = 17,colour = '#0000cc',size = 3.5) +
  geom_path(aes(x = block_energ_TOTAL,y = block_energ_EST),data=df.weir,size = 0.75) +
  geom_hline(aes(yintercept = yc.ec03),data=df.weir,colour = '#009999',size = 0.75,linetype = 6) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  geom_text(aes(x = block_energ_TOTAL,y = block_energ_EST,label = rule,hjust = 0),data=df.weir,size = 3.5,vjust = 1.0,parse = FALSE) +
  geom_text(aes(x = block_energ_TOTAL,y = block_energ_EST,label = block_ID,vjust = -1),data=df.weir,colour = '#cc0000',size = 5.5,parse = FALSE) +
  ggtitle("Energias especificas. Vertedor Triangular") +
  xlab("Energia Total (m)") +
  ylab("Energia Piezometrica(m)") +
  theme_bw(base_size = 18.0)

# A ggplot object is requested
fg02


# round_df function is applied to relevant data.frames
df.output <- round_df(df=df.weir, digits=3)

# Objects to export:
# yc.ec03, yc.ec04, yc.inter, yc.inter.dist, Ec, df.output, fg01, fg02
write.csv(df.output, file = "df.output.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////


