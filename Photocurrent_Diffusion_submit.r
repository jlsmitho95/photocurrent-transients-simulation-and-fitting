########################################
##Kinetics of First order recombination#
##This model includes the diffusion    #
########################################
# load libraries
library(deSolve)     # library for solving differential equations
library(ReacTran)    #Generates the grid for mass transport
library(FME)         #Performs nonlinear fitting
library(plotly)
library(RColorBrewer)
library(gganimate)
library(gifski)
library(transformr)
rm(list = ls())

cols <- brewer.pal(12, "Paired")
#############################
##Parameters
#The values are arbitrally selected
k<-c(Sig=1e5,     #Excitation constant, units = s^-1
     kd=5e5,      #Relaxation constant, units = s^-1
     ket1=1.5e5,  #Heterogenous electron transfer constant I, units = cm^2.s^-1.mol^-1
     kps1=1,      #Products separation constant I, units = s^-1
     kb1=3,       #Recombination constant I, units = s^-1
     ket2=1.5e5,  #Heterogeneous electron transfer constant II, units = cm^2.s^-1.mol^-1
     kb2=3,       #Recombination constant II, units = s^-1
     kps2=1)      #Products separation constant II, units = s^-1

#############################################
##Further constants
D.coeff   <-0.0000014         #Diffusion coefficient, units = cm^2.s^-1
Faraday   <-96485             #Faraday's constant, units = C.mol^-1
R         <-8.314462618       #molar Boltzmann constant R = 8.314462618 J.mol-1.K-1, taken from NIST
Temp      <-293.15            #Thermodynamic temperature in Kelvin 20 °C
f         <-Faraday/(R*Temp)  # f = F/RT in 
Galv.Pot  <- 0.1              #Galvani potential -0.2
Galv.Pot  <- Galv.Pot*f       #Dimensionless Galvani potential
Cb.W      <- 10e-3            #Aqueous phase supporting electrolyte concentration mol.L^-1. Ivan's thesis
E.0       <- 8.8541878128e-12 #Vacuum electric permittivity in F.m^-1.

##Debye-Huckel values for organic component##
E.TFT     <- 9.22             #Dielectric constant of TFT at 20 °C. 9.22 at 298.2 K Source: CRC Handbook 2004, pp 6-159
# E.TFT   <- k.TFT*E.0        #TFT permittivity at 20 °C
Cb.org    <- 5e-3             #Organic phase supporting electrolyte concentration mol.L^-1. Ivan's thesis
kappa.TFT <- ((2*1000*Faraday^2*Cb.org)/(E.TFT*E.0*R*Temp))^(1/2) #Reciprocal Debye-Length of organic phase in m^-1.
C.Org.DH  <- kappa.TFT*E.TFT*E.0/100  #Debye-Huckel value for organic phase.

##Debye-Huckel values for the aqueous component##
E.Water   <- 81.39             #Dielectric constant of water at 290 K, at 295K is 79.55. Source: CRC Handbook 2004, pp 6-13
Cb.W      <- 10e-3             #Aqueous phase supporting electrolyte concentration mol.L^-1. Ivan's thesis
kappa.W   <- ((2*1000*Faraday^2*Cb.W)/(E.Water*E.0*R*Temp))^(1/2) #Reciprocal Debye-Length of aqueous phase. Note: this is an approximation, the lithium citrate permittivity must be calculated or found.
C.W.DH    <- kappa.W*E.Water*E.0/100  #Debye-Huckel value in cm^-1 for aqueous phase.


#C.Donor<-0.0000014
N       <-100        #Lenght of space partition
xgrid   <-setup.grid.1D(x.up=0,x.down=0.01,N=N) #Generating the grid  #x.down = 0.1, discretizacion del espacio, metodo de las lineas generar una particion del espacio dx, x.down= 0.1 cm  
# plot(xgrid$x.mid)
x       <-xgrid$x.mid



############################################
# #***Initial concentrations at t=0***
yini                  <-rep(0.0,9*N)      #Initital concentration is set to zero. vector for 9 species
yini[1:N]             <-1.94e-9           #ZnPor   1. Initial concentration at the interface of porphyrin film. mol.cm3^-1
yini[(N+1):(2*N)]     <-0.0               #ZnPor*  2. Initial concentration of porphyrin fil in the photo-excited state. mol.cm3^-1
yini[(2*N+1):(3*N)]   <-5e-6              #D       3. Initial concentration of the electrodonor. mol.cm3^-1
yini[(3*N+1):(4*N)]   <-0.0               #[ZnP-D] 4. Initial concentration interfacial complex. mol.cm3^-1
yini[(4*N+1):(5*N)]   <-0.0               #ZnPor-  5. Initial concentration photoproduct D. mol.cm3^-1
yini[(5*N+1):(6*N)]   <-0.0               #D+      6. Initial concentration photoproduct D. mol.cm3^-1
yini[(6*N+1):(7*N)]   <-3.5e-7            #A       7. Initial concentration of the electracceptor. mol.cm3^-1
yini[(7*N+1):(8*N)]   <-0.0               #[ZnP-A] 8. Initial concentration interfacial complex. mol.cm3^-1
yini[(8*N+1):(9*N)]   <-0.0               #A-      9. Initial concentration of the photoproduct A. mol.cm3^-1
######################################

Model1<-function (t,y,k){
###**Initial conditions***###
  y1<-y[1:N]             #ZnPor
  y2<-y[(N+1):(2*N)]     #ZnPor*
  y3<-y[(2*N+1):(3*N)]   #D
  y4<-y[(3*N+1):(4*N)]   #ZnPor.D
  y5<-y[(4*N+1):(5*N)]   #ZnPor-
  y6<-y[(5*N+1):(6*N)]   #D+
  y7<-y[(6*N+1):(7*N)]   #A
  y8<-y[(7*N+1):(8*N)]   #ZnP.A
  y9<-y[(8*N+1):(9*N)]   #A+
  
  #Reactions rates
  v1<-k[1]*y1 - k[2]*y2      #σI0*[ZnPor] - kd[ZnPor*]
  v2<-k[3]*y2*y3[1]          #ket1[ZnPor][D]
  v3<-k[4]*y4                #kps1[ZnP.D]
  v4<-k[5]*y4                #kb1[ZnP.D]
  v5<-k[6]*y5*y7[1]-k[7]*y8  #ket2[ZnP-][A]-kb2[ZnP.A]
  v6<-k[8]*y8                #kps2[ZnP.A]
  
  #***Diffusion***
  tran1<-tran.1D(C=y3,flux.up=(k[5]*y4[1]-k[3]*y2[1]*y3[1]),flux.down=0.0,D=D.coeff,dx=xgrid$dx)  #Mass transport for D
  tran2<-tran.1D(C=y6,flux.up=(k[4]*y4[1]),flux.down=0.0,D=D.coeff,dx=xgrid$dx)                   #Mass transport for D+
  tran3<-tran.1D(C=y7,flux.up=(k[7]*y8[1]-k[6]*y5[1]*y7[1]),flux.down=0.0,D=D.coeff,dx=xgrid$dx)  #Mass transport for A
  tran4<-tran.1D(C=y9,flux.up=(k[8]*y8[1]),flux.down=0.0,D=D.coeff,dx=xgrid$dx)                   #Mass transport for A-
  
  #***Differential equations****
  dy1<--v1+v4+v6  #d[ZnPor]
  dy2<- v1-v2     #d[ZnPor*]
  dy3<- tran1$dC  #d[D]
  dy4<- v2-v3-v4  #d[ZnP.D]
  dy5<- v3-v5     #d[ZnP-]
  dy6<- tran2$dC  #d[D+]
  dy7<- tran3$dC  #d[A]
  dy8<- v5-v6     #d[ZnP.A]
  dy9<- tran4$dC  #d[A-]
  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in ydot
  return(list(c(dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9)))
}


##########################################################
##########################################################
####Function for generating the photocurrent####
Photocurrent     <-function(k){
  times          <-seq(0,10.0,by=0.05)    #The first 10 seconds 
  #time          <-seq(0,10.0,by=0.05)
  out1           <-ode.1D(y=yini,         #Initial concentration of species 
                          time=times,     #Time interval
                          func=Model1,    #Function with the differential equations
                          parms=k,        #Specified parameter k, rate constants
                          dimens=N,       #Number of data 
                          nspec=9,        #Number of chemical species 
                          method="lsode", #Method to solve the differential equations
                          names=c("ZnP","Znp*","D","ZnP.D","ZnP-","D+","A","ZnP.A","A+"))
  
  #**The light is off***
  # matplot(t(out1[,202:301]), type ="l")
  
  #The initial concentration
  yini2           <-out1[length(out1[,1]),2:length(out1[1,])]  # Is the concentration of all species at time = t
  k[1]=0.0        #No light
  out2            <-lsode(y=yini2,time=times,func=Model1,parms=k)               #ask Ivan 
  out             <-rbind(out1[1:(length(out1[,1])-1),],out2)                   #combine the data provided by ode.1D and lsode 
  # jph1            <-k[3]*out[,N+2]*out[,2*N+2] + k[6]*out[,4*N+2]*out[,6*N+2] #Rate of electrotransfer I and II
  # jph2            <--k[5]*out[,3*N+2] - k[7]*out[,7*N+2]                      #Rate of recombination I and II, we need to correct this to Jose's 
  # jph1            <-Faraday*jph1; jph2<-Faraday*jph2                          #Total current flux for et and b
  # jph             <-jph1+jph2                                                 #Sum of the current fluxes 
  time            <-seq(0,20.0,by=0.05)                                         #Sequence of time, the first 20 seconds
  sigma           <- -Faraday*(out[,4*N+2]+out[7*N+2])

  ##Function to calculate the potentials of equation h12#####  
  Potential     <- function(x,sigma){    
    (f/2)*sigma - C.Org.DH*sinh(x/2) + C.W.DH*sinh((Galv.Pot-x)/2)
  }
  Phi           <- vector(length = length(time))
  for (i in 1:length(time)){
    Phi[i]      <- uniroot(Potential, c(-3*Galv.Pot, 3*Galv.Pot), sigma[i],tol= 1e-18, check.conv = TRUE, trace=1)$root
  }
  
  Phi.W         <- Phi-Galv.Pot                   #Definition of the electric potential in the aqueous phase
  
  c.W           <- C.W.DH*cosh((Galv.Pot-Phi)/2)  #Differential capacitance of the wEDL
  C.O           <- C.Org.DH*cosh(Phi/2)           #Differential capacitance of the oEDL
  
  Psi.O         <- C.O/(C.O+c.W)                  #Relative contribution of the oEDL to the interfacial capacitance
  
  vpsI          <- k[4]*out[,3*N+2]               #v product separation (I) = kps^i*[ZnP.D]
  vpsII         <- k[8]*out[,7*N+2]               #v product separation (II) = kps^ii*[ZnP.A-]
  
  ##Photocurrents calculation ##
  
  j_ph          <- (1-Psi.O)*Faraday*vpsI + Psi.O*Faraday*vpsII   #Photocurrents as equation A32
  
  #Faradaic photocurrents# 
  
  j_f_aq        <- Faraday*vpsII                                  #Faradaic current of W phase 
  j_f_or        <- Faraday*vpsI                                   #Faradaic current of O phase
  
  #Capacitive photocurrents#            
  
  j_c_aq        <- j_ph - j_f_aq                                  #Capacitive current of W phase
  j_c_or        <- j_ph - j_f_or                                  #Capacitive current of O phase
  J_ph          <- cbind(time, j_ph,j_f_aq,j_f_or,j_c_aq,j_c_or)  #Vector with time and photocurrents
  return(J_ph)
 }
  
##### Apply the function to generate the photocurrent ----

Photocurrents <- Photocurrent(k)


######Generate the data frames to export data ####
df1 <- as.data.frame(Photocurrents)
# df2 <- as.data.frame(results[[2]][,1:6])
# df3 <- as.data.frame(results[[3]][,1:6])
# df4 <- as.data.frame(results[[4]][,1:6])
# df5 <- as.data.frame(results[[5]][,1:6])
# df6 <- as.data.frame(results[[6]][,1:6])
# df7 <- as.data.frame(results[[7]][,1:6])
# df8 <- as.data.frame(results[[8]][,1:6])
plot(df1$time, df1$j_ph)

stop("lala")
######Save the data in .csv file####
# write.csv(df1, "General\\v1.csv", row.names=FALSE)
# lty=# write.csv(df2, "kb-i\\v2.csv", row.names=FALSE)
# write.csv(df3, "kb-i\\v3.csv", row.names=FALSE)
# write.csv(df4, "kb-i\\v4.csv", row.names=FALSE)
# write.csv(df5, "kb-i\\v5.csv", row.names=FALSE)
# write.csv(df6, "kb-i\\v6.csv", row.names=FALSE)
# write.csv(df7, "kb-i\\v7.csv", row.names=FALSE)
# write.csv(df8, "kb-i\\v8.csv", row.names=FALSE)


##### Calculation of percetanges of contribution photocurrents###
# j               <- abs(Photocurrents[,3]) + abs(Photocurrents[,5])
# Cont_F_w        <- 100*abs(Photocurrents[,3])/j
# Cont_C_w        <- 100*abs(Photocurrents[,5])/j

##### Calculation of the contribution of photocurrents ----- 

#### Contribution in the aqueous phase --
Cont_F_w        <- 100*Photocurrents[,3]/Photocurrents[,2]
Cont_C_w        <- 100*Photocurrents[,5]/Photocurrents[,2]

#### Contribution in the organic phase ---
Cont_F_O        <- 100*Photocurrents[,4]/Photocurrents[,2]
Cont_C_O        <- 100*Photocurrents[,6]/Photocurrents[,2]


#### Generate the graphs related to the simulation of photocurrents transients ---- 

#### 1. Contribution of photocurrents at the aqueous phase ----

# jpeg("Figures/contribution-aqueous.jpeg", res=600,width = 5, height = 4, units = 'in')
# par(family = "sans")#, ps =10)
# plot(Photocurrents[,1], Cont_F_w, ylim = c(-400, 400), type = "l", col = cols[2],
#      ylab = "Contribution (%)", xlab = "Time (s)", lwd = 2, las = 1, tck = 0.05,
#      xaxs = "i", yaxs = "i", yaxp  = c(-400, 400, 10), cex =2.3)
# lines(Photocurrents[,1],Cont_C_w, ylim = c(min(Cont_C_w[2:401]), max(Cont_F_w[2:401])), col = cols[4], lwd = 2)
# lines(Photocurrents[,1], (Cont_F_w+Cont_C_w), col = cols[6], lty = 2, lwd=2)
# legend("topleft", legend = c(expression(italic(j)[F]^{w}),expression(italic(j)[F]^{w})), lty = c(1,1),
#        y.intersp=1.4, bty="n", col = c(cols[2], cols[4]), lwd =2)
# dev.off()

#### 2. Contribution of photocurrents at the organic phase ----

# jpeg("Figures/contribution-organic.jpeg", res=300,width = 5, height = 4, units = 'in')
# plot(Photocurrents[,1], Cont_F_O, ylim = c(1.5*min(Cont_C_O[2:401]), 1.1*max(Cont_F_O[2:401])), type = "l", col = cols[2],
#      ylab = "Contribution (%)", xlab = "Time (s)", lwd = 2)
# lines(Photocurrents[,1],Cont_C_O, ylim = c(min(Cont_C_O[2:401]), max(Cont_F_O[2:401])), col = cols[4], lwd = 2)
# lines(Photocurrents[,1], (Cont_F_w+Cont_C_w), col = cols[6], lty = 2, lwd=2)
# legend("right", legend = c(expression(italic(j)[F]^{o}),expression(italic(j)[C]^{o}), expression(italic(j)[F]^{o} +italic(j)[C]^{o})), lty = c(1,1,2),
#        y.intersp=1.4, bty="n", col = c(cols[2], cols[4],cols[6]), lwd =2)
# dev.off()


##### The simulated photocurrent transient graph ----- 
p1 <- ggplot(aes(x=time, y= 1e6*j_ph), data = df1) +
  geom_line(size =1.2, color = cols[2]  )+
  # geom_line(aes(y = 1e6*Experimental, color =Potential, linetype = "Exp."), data = full_data, size =1.2)+
  scale_x_continuous(name = bold("Time (s)"), 
                     breaks=seq(0,20,5), expand = c(0,0.1), 
                     minor_breaks = seq(0, 20, 2.5),
                     limits = c(0,20),sec.axis = dup_axis(),
                     guide = guide_axis(minor.ticks = TRUE))+
  scale_y_continuous(name = expression(bold(italic(j)[photo]~(mu*A%.%cm^{-2}))), 
                     breaks = seq(-1,10,1),
                     minor_breaks = seq(-1,10, 0.5), 
                     limits = c(-0.02,
                                1.1*1e6*max(df1$j_ph)
                     ),
                     sec.axis = dup_axis(),
                     expand = c(0,0), 
                     guide = guide_axis(minor.ticks = TRUE))+
  # scale_color_viridis("Scan", discrete = FALSE, option = "G") +
  # scale_linetype_manual("" , values = c("Fitted" =1,
  #                                       "Exp." = 2)
  # )+
  # scale_colour_manual(expression(bold(Delta[o]^w*italic(ϕ)~"(V)")),values = c(brewer.pal(n=8, name="Dark2")))+
  theme(text=element_text(family="Arial"),
        plot.margin=unit(c(0.35,0.5,0.1,0.1), 'cm'), #t, r, b, l
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size =20, face="bold"),
        axis.title.x = element_text(size =16, face="bold"),
        axis.text=element_text(size=20, face= "bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks.length= unit(-0.15, "cm"), #Axis tick
        axis.text.x.top = element_blank(),       # do not show top / right axis labels
        axis.text.y.right = element_blank(),     # for secondary axis
        axis.title.x.top = element_blank(),      # as above, don't show axis titles for
        axis.title.y.right = element_blank(),   # secondary axis either
        legend.position = "right",#c(.15,.85),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "cm"), #change legend key size
        legend.key.width = unit(2, "line"),
        legend.background = element_rect(fill = NA),
        legend.key= element_blank(),
        legend.text = element_text(size = 15, face= "bold"),
        legend.title = element_text(size = 15, face= "bold"),# element_blank(),
        plot.title = element_text(size=14.8, hjust = 1))

p1

### Animation of the photocurrent transient graph ----
anim_p1 <- p1 + transition_reveal(time)

animate(
  anim_p1,
  nframes = 100,
  fps = 25,
  width = 1800,
  height = 1500,
  res = 300,
  renderer = gifski_renderer("j_ph_time.gif")
  
)

# anim_p1

stop("Stop here")

my_label_positions <- seq(-1, 5,by = 1)
my_labels <- sprintf(my_label_positions,    # Specify axis labels
                     fmt = '%#.1f')


###### 3. Simulated photocurrent transient at the aqueous phase -----

# jpeg("Figures/general-photocurrent-water.jpeg", res=300,width = 5, height = 4, units = 'in')
# par(mar = c(3.5, 5, 1, 2), family = "sans", font.axis =2) #pty="s"mar = c(5.5, 6.5, 1, 2) +0.3
# plot(x = Photocurrents[,1], y= Photocurrents[,2]*1e6 , xlab = "Time (s)", ylab = expression(italic(j)[ph]~(mu*A%.%cm^{-2})),
#      cex.lab = 1.5, cex.axis = 1.0, lwd= 3,
#      type = "l", col = cols[2], mgp = c(2.5, 0.5, 0), tck = 0.02, 
#      main = " ", cex.main = 2.0,xaxs = "i", yaxs = "i",yaxt ="n", ylim = c(-1,5))
# axis(side = 2, at = my_label_positions, labels= my_labels, tck = 0.02, las = 1, cex = 1.5,yaxs = "i", mgp = c(1.5, 0.5, 0))
# axis(side = 3, labels = FALSE, tck = 0.02)
# axis(side = 4, labels = FALSE, tck = 0.02)
# lines(x = Photocurrents[,1], y= Photocurrents[,3]*1e6 , col = cols[4], lty = 2, lwd = 3)
# lines(x = Photocurrents[,1], y= Photocurrents[,5]*1e6 , col = cols[6], lty = 3, lwd = 3)
# legend("topright", legend = c(expression(italic(j)[ph]), expression(italic(j)[F]^w), expression(italic(j)[C]^w)),
#        bty ="n", col = c(cols[2], cols[4], cols[6]), lty = c(1,2,3), lwd = 3, cex = 1.5, y.intersp=1.4)
# dev.off()

###### 4. Simulated photocurrent transient at the organic phase -----
# jpeg("Figures/general-photocurrent-organic.jpeg", res=300,width = 5, height = 4, units = 'in')
# my_label_positions <- seq(-2, 6,by = 1)
# my_labels <- sprintf(my_label_positions,    # Specify axis labels
#                      fmt = '%#.1f')
# par(mar = c(3.5, 5, 1, 2), family = "sans", font.axis =2) #pty="s"mar = c(5.5, 6.5, 1, 2) +0.3
# plot(x = Photocurrents[,1], y= Photocurrents[,2]*1e6 , xlab = "Time (s)", ylab = expression(italic(j)[ph]~(mu*A%.%cm^{-2})),
#      cex.lab = 1.5, cex.axis = 1.0, lwd= 3,
#      type = "l", col = cols[2], mgp = c(2.5, 0.5, 0), tck = 0.02,
#      main = "", cex.main = 2,xaxs = "i", yaxs ="i", yaxt ="n", ylim = c(-2,6) )
# axis(side = 2, at = my_label_positions, labels= my_labels, tck = 0.02, las = 1, cex = 1.5,yaxs = "i", mgp = c(1.5, 0.5, 0))
# axis(side = 3, labels = FALSE, tck = 0.02)
# axis(side = 4, at = my_label_positions, labels=FALSE, tck = 0.02)
# lines(x = Photocurrents[,1], y= Photocurrents[,4]*1e6 , col = cols[4], lty = 2, lwd = 3)
# lines(x = Photocurrents[,1], y= Photocurrents[,6]*1e6 , col = cols[6], lty = 3, lwd = 3)
# legend("topright", legend = c(expression(italic(j)[ph]), expression(italic(j)[F]^o), expression(italic(j)[C]^o)),
#        bty ="n", col = c(cols[2], cols[4], cols[6]), lty = c(1,2,3), lwd = 3, cex = 1.5, y.intersp=1.4)
# dev.off()
####
##Photocurrent is generted
# Current    <-Photocurrent(k)  #provides the current flux in uA.cm^-2









