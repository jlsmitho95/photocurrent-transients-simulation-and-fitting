# Set working directory (this line is commented out; add setwd("path") if needed)

# load libraries
rm(list = ls())
library(deSolve)     # library for solving differential equations
library(ReacTran)    #Generates the grid for mass transport
library(FME)         #Performs nonlinear fitting
library(plotly)
library(RColorBrewer)#To add colors to the plot
library("dplyr")
# library("ggplot2")
library(shape)       #To draw arrows in the plots

# Clear the workspace
# rm(list = ls())

## =================================================
## Load experimental data
## =================================================
Experimental <- read.table("Data_corrected_12-02-2020-3-1.94-5-DP-1.csv", header = TRUE)

# Determine dimensions of the dataset
n <- length(Experimental[, 1])  # Number of rows
m <- length(Experimental[1, ])  # Number of columns

# Time vector (from 0 to 10 seconds, step of 0.05s)
time <- seq(0, 10, by = 0.05)

# Standard deviation vector for weighting in fitting
sd <- seq(1, length(time))
sd[1:40] <- 0.5               # Higher uncertainty at early times
sd[40:length(time)] <- 0.8  # Lower uncertainty at later times

# Create observed dataset for a specific potential (column 14, rows 200–400)
Obs <- cbind(time, Experimental[200:400, 16], sd)
colnames(Obs) <- c("time", "j_ph", "sd")

# Plot observed photocurrent
plot(Obs, xlab = "Time / s", ylab = expression("j / uA"%.%"cm"^{-2}, lwd = 3))

## =================================================
## Define model parameters and initial conditions
## =================================================
# Load kinetic parameters from CSV
# kt <- read.table("kf.csv", header = TRUE, sep = ",")
# kt <- as.numeric(kt$value)

# Assign named kinetic rate constants
k   <- c(
  Sig  = 1.0e5,      # Excitation rate
  kd   = 5.0e5,#50569.53,   # Decay of excited state
  ket1 = 3.614610e+05,    #1.452293e+05# Electron transfer rate to donor
  kps1 = 4.428191e+00,   # Product separation from donor complex
  kb1  = 1.847824e+01,         # Back reaction rate for donor complex
  ket2 = 1.361810e+05,    # Electron transfer to acceptor
  kb2  = 8.154705e-02,      #1.57 # Back reaction for acceptor complex
  kps2 = 1.725764e+01     # Product separation from acceptor complex
)

##Note: The last values were found in the data fitting at a potential of 0.75 V
## To add the constant parameters ---
####### CONSTANTS ######################################
##Further constants
D.coeff   <-0.0000014         #Diffusion coefficient, units = cm^2.s^-1
Faraday   <-96485             #Faraday's constant, units = C.mol^-1
R         <-8.314462618       #molar Boltzmann constant R = 8.314462618 J.mol-1.K-1, taken from NIST
Temp      <-293.15            #Thermodynamic temperature in Kelvin 20 °C
f         <-Faraday/(R*Temp)  # f = F/RT in 
Galv.Pot  <- 0.85             #Galvani potential
Galv.Pot  <- Galv.Pot*f       #Dimensionless Galvani potential
Cb.W      <- 10e-3            #Aqueous phase supporting electrolyte concentration mol.L^-1. Ivan's thesis
E.0       <- 8.8541878128e-12 #Vacuum electric permittivity in F.m^-1.

##Debye-Huckel values for organic component##
E.TFT     <- 9.22             #Dielectric constant of TFT at 20 °C. 9.22 at 298.2 K Source: CRC Handbook 2004, pp 6-159
Cb.org    <- 5e-3             #Organic phase supporting electrolyte concentration mol.L-1. Ivan's thesis
kappa.TFT <- ((2*1000*Faraday^2*Cb.org)/(E.TFT*E.0*R*Temp))^(1/2) #Reciprocal Debye-Length of organic phase in m^-1.
C.Org.DH  <- kappa.TFT*E.TFT*E.0/100  #Debye-Huckel value for organic phase.


##Debye-Huckel values for the aqueous component##
E.Water   <- 81.39             #Dielectric constant of water at 290 K, at 295K is 79.55. Source: CRC Handbook 2004, pp 6-13
Cb.W      <- 10e-3             #Aqueous phase supporting electrolyte concentration mol.cm^-3. Ivan's thesis
kappa.W   <- ((2*1000*Faraday^2*Cb.W)/(E.Water*E.0*R*Temp))^(1/2) #Reciprocal Debye-Length of aqueous phase. Note: this is an approximation, the lithium citrate permittivity must be calculated or found.
C.W.DH    <- kappa.W*E.Water*E.0/100  #Debye-Huckel value in cm^-1 for aqueous phase.

# Diffusion and spatial grid setup
N       <-100        #Lenght of space partition
xgrid   <-setup.grid.1D(x.up=0,x.down=0.01,N=N) #Generating the grid  #x.down = 0.1, discretizacion del espacio, metodo de las lineas generar una particion del espacio dx, x.down= 0.1 cm  
x       <-xgrid$x.mid

#######***Initial concentrations at t=0#####
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

## =================================================
## Define the reaction-diffusion model
## =================================================

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


## =================================================
## Photocurrent simulation function
## =================================================
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
  
  #The initial concentration
  # yini2           <-out1[length(out1[,1]),2:length(out1[1,])]  # Is the concentration of all species at time = t
  # k[1]=0.0        #No light
  # out2            <-lsode(y=yini2,time=times,func=Model1,parms=k)               #ask Ivan 
  out             <-out1 #rbind(out1[1:(length(out1[,1])-1),],out2)                   #combine the data provided by ode.1D and lsode 
  # jph1          <-k[3]*out[,N+2]*out[,2*N+2] + k[6]*out[,4*N+2]*out[,6*N+2] #Rate of electrotransfer I and II
  # jph2          <--k[5]*out[,3*N+2] - k[7]*out[,7*N+2]                      #Rate of recombination I and II, we need to correct this to Jose's 
  # jph1          <-Faraday*jph1; jph2<-Faraday*jph2                          #Total current flux for et and b
  # jph           <-jph1+jph2                                                 #Sum of the current fluxes 
  time            <-times #seq(0,10.0,by=0.05)                                #Sequence of time, the first 20 seconds
  sigma           <- -Faraday*(out[,4*N+2]+out[7*N+2])
  
  ##Function to calculate the potentials of equation h12#####  
  Potential     <- function(x,sigma){    
    (f/2)*sigma - C.Org.DH*sinh(x/2) + C.W.DH*sinh((Galv.Pot-x)/2)
  }
  Phi           <- vector(length = length(time))
  for (i in 1:length(time)){
    Phi[i]      <- uniroot(Potential, c(-5*Galv.Pot, 5*Galv.Pot), sigma[i],tol= 1e-18, check.conv = TRUE, trace=1)$root
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
  j_c_or        <- j_ph - j_f_or                                  #Capacitive current of W phase
  J_ph          <- cbind(time, j_ph,j_f_aq,j_f_or,j_c_aq,j_c_or,Phi, Psi.O)  #Vector with time and photocurrents
  # Jph             <-cbind(time,jph,jph1,jph2)
  #Jph<-cbind(time,jph)
  return(J_ph)
}



Photocurrents <- Photocurrent(k)

plot(Photocurrents[,1], Photocurrents[,2], type = "l", col = "blue", lwd = 2,
     xlab = "Time / s", ylab = expression("j / uA"%.%"cm"^{-2}),
     main = "Photocurrent Simulation", ylim = c(0, 2*max(Obs[,2])))
par(new=TRUE)
lines(Obs[,1], Obs[,2], col = "red", lwd = 2, ylim = c(0, 2*max(Obs)))



## =================================================
## Model cost function (goodness-of-fit)
## =================================================
Jphcost <- function(ki) {
  out <- Photocurrent(ki)
  cost <- modCost(model = out, obs = Obs, err = "sd")  # Weighted residuals
  return(cost)
}

# Print model prediction output
print(Jphcost(k)$model)
Residuals <- Jphcost(k)  
plot(Residuals$residuals$res)

### What is the meaning that residuals are not toally zero and not a random trend###


# plot(Residuals, ylim= c(-10e-6,3e-6))


## =================================================
## Local sensitivity analysis
## =================================================
Sfun <- sensFun(Jphcost, k)
print(summary(Sfun))

## Collinearity analysis 
ident <- collin(Sfun)
head(ident, n = 20)
plot(ident, log = "y")
collin(Sfun, parset = c("Sig", "kd", "ket1", "kps1", "kb1", "ket2", "kb2", "kps2"))
collinear <- collin(Sfun,  N = 2)

stop("lala")

 
###  Parameter sampling using Latin Hypercube ----

# Define log-scale bounds for sampling
# low <- c(-3, -14, -9, -9, -14, -14, -14, -14)
# up  <- c(5,   5,   7,  7,  10,  10,  10,  10)
# low <- c(5, 5, 4, -7, -5, 4, -5, -7)
# up  <- c(18,   18,   17,  5,  6,  17,  7,  7)

# low <-  log(c(1.7e5, 4.9e4,1.2e5, 0.5, 5, 1e5, 0.9, 7.2))
# up  <-  log(c(2.1e5, 5.4e4,1.5e5, 1, 8.0, 1.6e5, 1.6, 9)) 
low <-  log(c(1,   #k5
              1e-3    #k7
              ))

up  <-  log(c(20,
              20
              )) 

parRange <- data.frame(min = low, max = up)

rownames(parRange) <- c("k5", "k7")   #c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8")

# Generate 50 parameter sets
SampleP   <- Latinhyper(parRange, 10)
FinalResults <- matrix(ncol = (length(low) + 2), nrow = 10)     

# i <- 1
# 
# print(system.time(
#   while (i <= 50) {
#     Modelled  <- Photocurrent(c(k[1], 
#                                 k[2],
#                                 exp(SampleP[i,1]),
#                                 exp(SampleP[i,2]),
#                                 k[5],
#                                 exp(SampleP[i,3]), 
#                                 k[7],
#                                 exp(SampleP[i,4]))) # Simulate the model using initial parameters
#     
#     i <- i + 1
#  


# ### Plot experimental data + modelled to find some boundaries ---- 
# plot(Obs, xlab = "Time / s", ylab = expression("j / uA"%.%"cm"^{-2}), 
#      main = "Observed Photocurrent")
# par(new=TRUE)
# lines(Modelled[,1], Modelled[,2], col = "red",lwd=2)
#   }
# ))




## =================================================
## Model fitting loop
## =================================================
Jphcost2 <- function(pars) {
  new_pars  <- c(k[1],
                 k[2],
                 k[3],
                 k[4],
                 exp(pars[1]),
                 k[6],
                 exp(pars[2]),
                 k[8])
  Jphcost(new_pars)  # Exponentiate for log-scale fitting
}

# Fit model for each sampled parameter set
# Fit <- modFit(f = Jphcost2, p = SampleP[1, ], method = "Marq", upper = up, lower = low)

i <- 1

print(system.time(
  while (i <= 10) {
    Fit <- modFit(f = Jphcost2, p = SampleP[i, ], method = "Marq", upper = up, 
                  lower = low, control = list(maxiter = 500))
    FinalResults[i, ] <- c(i, Fit$ssr, Fit$par)
    print(i)
    i <- i + 1
  }
))

summary(Fit)  # Print summary of the fit
# Plot log of sum of squared residuals
# plot(log(FinalResults[, 2]))

# S <- sensFun(Jphcost2, Fit$par)
# summary(S)

# Print the coefficients obtained from the most recent model fitting (in log-scale)
# print(coef(Fit))

# Convert the log-scale fitted parameters to linear scale using exponentiation
kf <- c(k[1], 
        k[2],
        k[3],
        k[4],
        exp(coef(Fit)[1]),
        k[6],
        exp(coef(Fit)[2]),
        k[8]
        )
        
  #exp(coef(Fit)[1:8]))

# Print the original parameter values (before fitting)
print(k)

# Print the fitted parameter values (after exponentiation)
print(kf)
min_Fit <- FinalResults[which.min(FinalResults[,2]), ]
print(min_Fit)
#### Simulate the model using initial parameters (before fitting)-----
ini <- Photocurrent(k)
# Simulate the model using the fitted parameters
final <- Photocurrent(kf)

## Plot the observed photocurrents and the fitted data ---
plot(Obs, xlab = "time/s", ylab = "Jph", ylim= c(0,1.1*max(final[,2])))  # Optionally, limits can be set using ylim/xlim
# Add a dashed line representing the model prediction using initial parameters
par(new=TRUE)
lines(ini[,1], ini[,2], ylim= c(0,1.1*max(final[,2])),lty = 2, lwd = 2)
# Add a solid red line representing the model prediction using fitted parameters
par(new=TRUE)
lines(final[,1], final[,2], ylim= c(0,1.1*max(final[,2])), col = "red", lwd = 2)
# Add a legend to the plot for clarity
legend("bottomright", 
       c("data", "initial", "fitted"), 
       lty = c(NA, 2, 1), 
       pch = c(1, NA, NA),
       col = c("black", "black", "red"), 
       lwd = c(NA, 2, 2))




# Residuals <- Jphcost(kf)
# 
# # Perform local sensitivity analysis for the fitted parameters
# Sfun <- sensFun(Jphcost, kf)
# # Check parameter collinearity (to assess identifiability)
# ident <- collin(Sfun)


stop("lala")

#### Lets make a new NLR with the best parameter found in the last step -----
Fit2   <- modFit(f = Jphcost2, p = FinalResults[min_Fit[1],3:4] , 
                 method = "Marq", 
                 upper = up, 
                 lower = low,
                 control = list(maxiter = 500))


final2 <- Photocurrent(c(k[1], 
                         k[2],
                         k[3],
                         k[4],
                         exp(coef(Fit2)[1]),
                         k[6],
                         exp(coef(Fit2)[2]),
                         k[8]
                         ))

### Print the summary of the last data fitting ----
summary(Fit2)

### Print the sum of squares due to regression
print(Fit2$ssr)


# Convert the log-scale fitted parameters to linear scale using exponentiation
kf2 <- c(k[1],
         k[2],
         k[3],
         k[4],
         exp(coef(Fit2)[1]),
         k[6],
         exp(coef(Fit2)[2]),
         k[8]
         )


names(kf2) <- names(kf)
# Print the original parameter values (before fitting)
# print(k)

# Print the fitted parameter values (after exponentiation)
print(kf2)


# S <- sensFun(Jphcost, kf2)
# summary(S)

####New collinearity analysis with the fitted data --
# ident2 <- collin(S, parset = c("Sig", "k2", "ket1", "kps1", "k5", "ket2", "kb2", "k8"))
# plot(ident2, log = "y")
# collinear2 <- collin(Sfun,  N= 3)


## Plot the observed photocurrents and the fitted data ---
plot(Obs, xlab = "time/s", ylab = "Jph", ylim= c(0,1.1*max(Obs[,2])))  # Optionally, limits can be set using ylim/xlim
# Add a dashed line representing the model prediction using initial parameters
par(new=TRUE)
lines(final2[,1], final2[,2], ylim= c(0,1.1*max(Obs[,2])), col = "red", lwd = 2)

stop("End of second fitting")


#### bind all the data in one data frame
Obs_final          <- as.data.frame(Obs[,1:2])
Obs_final$group    <- rep("Experimental", length(Obs_final[,2]))
final2             <- as.data.frame(final2[,1:2])
final2$group       <- rep("Fitted", length(Obs[,2]))
final_data         <- rbind(Obs_final, final2)

#### Save the fitted data in a .txt archive 
# Add the new column 'Potential' with constant value 0.85
final_data$Potential <- 0.85

## Save the data ---
write.table(final_data, file = "Fitting2/0.85V/Fitted-0.85V-k57-03112025-2ndtry.txt", sep ="\t")

### Plot the data --- 
p1 <- ggplot(aes(time, 1e6*j_ph, linetype =group, color = group, group = group, linewidth = group), data = final_data)+
  geom_line()+
  scale_x_continuous(name = "Time (s)", 
                     breaks=seq(0,10,1), expand = c(0,0.1), 
                     minor_breaks = seq(0, 10, 0.5),
                     limits = c(0,10),sec.axis = dup_axis(),
                     guide = guide_axis(minor.ticks = TRUE))+
  scale_y_continuous(name = expression(bold(italic(j)[photo]~(mu*A%.%cm^{-2}))), 
                     breaks = seq(0,10,2),
                     minor_breaks = seq(0,10, 1), 
                     limits = c(0,6.5),
                     sec.axis = dup_axis(),
                     expand = c(0,0), 
                     guide = guide_axis(minor.ticks = TRUE))+
  scale_linetype_manual(" ", values = c(3,1))+
  scale_colour_manual(" ",values = c("black",brewer.pal(n=5, name="Dark2")[1:2]))+
  scale_linewidth_manual(" ", values = c(2,1)) +
  theme(text=element_text(family="Arial"),
        plot.margin=unit(c(0.1,0.5,0.1,0.1), 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size =20, face="bold"),
        axis.title.x = element_text(size =18, face="bold"),
        axis.text=element_text(size=20, face= "bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks.length= unit(-0.15, "cm"), #Axis tick
        axis.text.x.top = element_blank(),       # do not show top / right axis labels
        axis.text.y.right = element_blank(),     # for secondary axis
        axis.title.x.top = element_blank(),      # as above, don't show axis titles for
        axis.title.y.right = element_blank(),   # secondary axis either
        legend.position = c(0.8,0.70), #"right",#c(.16,.76),
        legend.direction="vertical",
        legend.key.size = unit(0.5, "cm"), #change legend key size
        legend.key.width = unit(2, "line"),
        legend.background = element_rect(fill = NA),
        legend.key= element_blank(),
        legend.text = element_text(size = 15, face= "bold"),
        legend.title = element_text(size = 15, face= "bold"),# element_blank(),
        plot.title = element_text(size=14.8, hjust = 1))
  
p1

ggsave(plot = p1,
       width = 6,
       height = 5,
       dpi = 300,
       filename = "Fitting3/0.85V/Photocurrent-fitted-0.85V-k57-16122025.png")

