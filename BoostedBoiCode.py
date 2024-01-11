# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 08:15:41 2020

@author: Jacob Munnicha
"""
"BOOSTED CBA"

def cp_equ (T):
    #Calculate the equivalent specific heat of the air (working fluid) in the 
    #cylinder at a given temperature. 
    #From Table A.9 in Moran (Slide 28 in Slide Vault)
    theta = T / 100.
    cp0_N2 = (39.060 - 512.79 * theta ** -1.5 + 1072.7 * theta ** -2. -
             820.40 * theta ** -3.)
    cp0_O2 = (37.434 + 0.0201 * theta ** 1.5 - 178.57 * theta ** -1.5 +
              236.88 * theta ** -2.0)
    
    cpq = (0.79 / 28) * cp0_N2 + (0.21 / 32) * cp0_O2
    
    return cpq

#Full displacement in meters^3
Vd =(float(input("Engine Displacment [L]:\t")) * 10**-3) 

#Compression Ratio
"rc = 8.5 "
#Pressure Downstream of Aftercooler
P2=float(input("Pressure Downstream of Aftercooler [kPa]:\t"))

#Intake Temperature Downstream of aftercooler
T3=float(input("Intake Temp. Downstream of aftercooler[K]:\t"))

#Compressor intake Pressure kpa
P1=float(input("Pressure intake of Compressor [kPa]:\t"))

#Compressor intake Temp
T1=float(input("Intake Temp. of Compressor[K]:\t"))

#RPM of Max Power
N=float(input("Max RPM Speed:\t"))

#isentropic efficiency
N_sc =float(input("Compressor's Isentropic efficiency [Dec.]:\t"))
N_st =float(input("Turbine's Isentropic efficiency [Dec.]:\t"))

#Volumetric Efficiency
N_v=float(input("Volumetric efficiency [Dec.]:\t"))

#Mechanicals Efficiency
N_m=float(input("Mechanical Efficiency[Dec.]:\t"))

 # Density
cp = cp_equ(T3)
cp=1.005


#Gas Constant"
R= 287.05

#"density of air"
Row=input("Is density given Y/N?").upper()
if Row == 'Y':
    Row=float(input("What is it [kg/m^3]?\t"))
else:
    Row= P2*(10**3) / (R * T3)
    
print("\n\t\t\tUseful metrics")
print ("cp:%.4f [kJ/kg-K]"%cp)
print ("Density of air: %.3f [kg/m^3]" %Row )


#"Mass Flow rate for 4 stroke cycle"
M_a = N_v * Row * Vd * (N/60/2)
print ("Mass Flowrate: %.4f [kg/s]" %M_a)

#"Isentropic Leaving Temperature"
k=1.4
T_2s = T1* (P2/P1) ** ((k-1)/k)
T_2sc=T_2s - 273
print ("Isentropic Leaving Temperature: %.4f [K]"%T_2s)
print ("Isentropic Leaving Temperature: %.4f [C]"%T_2sc)

#"given Compressor Isentropic Efficiency Find T2 Actual K"
T_2a=  ((T_2s - T1) / N_sc) + T1
T_2ac=T_2a - 273
print ("Actual air Temp.of Compressor(T2a): %.4f [K]"%T_2a)
print ("Actual air Temp.or Compressor(T2a): %.4f [C]"%T_2ac)

#Turbine isentropic Efficiency
T_2aT= T1 - ((T1 - T_2s) * N_st) 
T_2acT=T_2aT - 273
print ("Actual air Temp. of Turbine(T2a): %.4f [K]"%T_2aT)
print ("Actual air Temp. of Turbine(T2a): %.4f [C]"%T_2acT)

# Heat Rejection from aftercooler -- Q = M_a * cp * Tout - T_cool
Q=input("Is the Mass FLow given Y/N?").upper()
if Q=='Y':
    M_a=float(input("What is the Mass Flow[Kg/s]?\t"))
    Tr=float(input("What is the reduced Temp [K]?\t "))
    Q= M_a * cp * ((T_2a) - (Tr))
    W_c = M_a *cp *(T_2a - T1) / N_m #Compressor Work:W= M_a * cp * (T_2a - T_ into comp) Change M_a if needed
else:
    Tr=float(input("What is the reduced Temp [K]?\t "))
    Q= M_a * cp * ((T_2a) - (Tr))
    W_c = M_a *cp *(T_2a -T1) / N_m  #Compressor Work:W= M_a * cp * (T_2a - T_ into comp) Change M_a if needed
print("\n\t\t\tUseful metrics")
print( "Heat work (Q): %.4f [kw]"%Q)
print("Compressor work (W-dot): %.4f [kw]"%W_c)

"""""""""""""""turbo"""""""""""""""""""""""""""""""""

#"Work of the Turbo
W_t  = W_c /N_m

#Mechanicals Efficiency
N_m= W_c/W_t

#Overall Turbo Efficiency
N_t= N_sc * N_st * N_m *100







