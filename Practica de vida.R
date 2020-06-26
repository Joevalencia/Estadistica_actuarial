#-----------------------
#       APARTADO 1
#-----------------------
install.packages('lifecontingencies')
library(lifecontingencies)

tabla <- read.csv2('PASEM2010.csv')
hombres <- tabla[-3]
names(hombres) <- c('x', 'qx')
hombresLt <- probs2lifetable(hombres[,"qx"],type="qx",name="PASEM 2010 Hombres")    #es un objeto lifetable
mujeres <- tabla[-2]
names(mujeres) <- c('x', 'qx')
mujeresLt <- probs2lifetable(mujeres[,"qx"],type="qx",name="PASEM 2010 MUJERES")    #es un objeto lifetable
tablaList <- list(hombresLt, mujeresLt)    #tabla de vida conjunta

#mujeresAct<-new("actuarialtable",x=mujeresLt@x,lx=mujeresLt@lx,interest=0.0115,
#                 name="MUJERES ACT TABLE")

#hombresAct<-new("actuarialtable",x=hombresLt@x,lx=hombresLt@lx,interest=0.0115,
#                 name="HOBRES ACT TABLE")

#------------------------------
#           APARTADO 2
#-----------------------------

# a) 

a2 <- qxt(hombresLt, x = 49, t = 2)
a2

# b)
b2 <- pxt(hombresLt, x=49, t=2)*qxt(hombresLt, x=51, t=3)
b2

# c)

c2 <- 50000*qxt(hombresLt,x=49,t=2)*(1.015)^(-1)+
  +65000*pxt(hombresLt,x=49,t=2)*qxt(hombresLt,x=51,t=3)*(1.015)^(-3)
round(c2, 2)

# d) 
d2 <- 1 - (a2 + b2 - a2*b2)
d2

#-------------------------------------
             #APARTADO 3
#-------------------------------------

# a)
a3 <- qxyzt(tablaList, x = c(40, 38), t = 2, status = 'joint')
a3

# b)
m_h <- qxt(hombresLt, 40, 2)    # muerte Hombres
v_m <- 1 + pxt(mujeresLt, 38, 2)   # vida Mujeres

b3 <- (m_h*v_m/2)
b3


# c)
c3 <- qxt(hombresLt, 40, 2) * pxt(mujeresLt, 38, 2)
c3

# d) 
d3 <- 50000*Axyzn(tablaList, x = c(40, 38), n = 2, i = 0.0115, status = 'joint')
round(d3, 2)

# e)
c <- 50000*qxt(hombresLt, x=40, t=1)*pxt(mujeresLt, x=38, t=1)*(1.0115)^(-1)
d <- 70000*pxt(hombresLt,x=40, t=1)*qxt(hombresLt, x=41,t=1)*pxt(mujeresLt, x=38, t=2)*(1.0115)^(-2)
e3 <- c+d
round(e3, 2)

#------------------------
#          APARTADO 4
#------------------------

# a) 
a4 <- 150000*Axyzn(tablaList, x = c(45, 43), n = 16, i = 0.012, status = 'last')
round(a4, 2)

# b) 
b4 <- 0
for (i in 1:16){
 
   b <- (2000*(i-1) + 150000)*(Axyzn(tablaList, x = c(45, 43), n=i, i = 0.012, status = 'last')-
        Axyzn(tablaList, x = c(45,43), n = i-1, i=0.012, status = 'last'))
        b4 <- b4 + b
}

b4
#-----------------------------------