# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:27:20 2020

@author: caden
"""

# -*- coding: utf-8 -*-"""Created on Fri Jan 25 07:47:37 2019
#******************************        Air standard OTTO        *******************************
#J.Rhett Mayor, Ph.D@author: jmayor3
import numpy as np
import matplotlib.pyplot as plt
# cp calculator
def cp_equ (T):
#Calculate the equivalent specific heat of the air (working fluid) in the 
#cylinder at a given temperature. 
#From Table A.9 in Moran (Slide 28 in Slide Vault)    
    theta = T / 100.    
    cp0_N2 = (39.060 - 512.79 * theta ** -1.5 + 1072.7 * theta ** -2. -820.40 * theta ** -3.)    
    cp0_O2 = (37.434 + 0.0201 * theta ** 1.5 - 178.57 * theta ** -1.5 +236.88 * theta ** -2.0)    
    cpq = (0.79 / 28) * cp0_N2 + (0.21 / 32) * cp0_O2
    return cpq
# =============================================================================
# MAIN PROGRAM
# =============================================================================
# =============================================================================
# Cylce Inputs 
# =============================================================================
# Cylinder geometries
C=8                #Number of cyclinders**
Vd= 7.98     # Total Displacement of engine in liters
Vs = Vd / C * 1e-3   #change the liters and the number of cyclinders it converts to m^3 this is swept Vs**
rc = 20           #compression ratio**
Vc = Vs / (rc - 1.)    #Clearance Volume

N=5500  # **Change the Rpm For top horsepower**
HP=395 # Given in Horsepower
Torque=HP*5252/ N    # Converts to lb-ft
Torque= (1/0.73756214927727)*Torque    #converts to NewtonMeters
# Combustion and intake efficiencies
eta_c = 1
xr = 0.04#Exhaust residual percentage
eta_s = 1 - xr

# Fuel in this analysis is Gasoline
LHV = 42.5# MJ/kg   light diesel is 42.5 and heavy is 41.4
psy = 0.5# Stoichioemetric ratio
AF_stoich = 14.9#Stoich AF ratio for Diesel
AF = AF_stoich / psy

# intake state
p1 = 80#intake pressure [kPa]  if naturally aspirated 
T1 = 20. + 273. #converted to kelvin from celsius
V1 = Vc + Vs
cp = cp_equ(T1)
cp=.819
R = 0.287# Gas constant for air kJ/kg
cv = cp - R
#k = cp / cv
k = 1.35

# Calculate the total mass in the cylinder (based on volume)
m1 = p1 * V1 / R / T1
# Correct for residual exhaust products and calculate air and fuel mass

m_m = m1
m_a = eta_s * m_m * (AF / (1 + AF))
m_f = eta_s * m_m * (1 / (1 + AF))

# Start with Power cycle
# Proces 1-2 : Isentropic Compression
p2 = p1 * (rc) ** k
T2 = T1 * (rc) ** (k - 1)
V2 = Vc

w12 = cv * (T1 - T2)
W12 = m1 * R * (T2 - T1) / (1 - k)

# Instant. combustion Constant Pressure process
Qc = eta_c * m_f * LHV * 1e3
L=cp_equ(T2)
T3 = T2 + Qc / m1 / 1.109
p3 = p2
B=T3/T2 # Cutoff ratio 
V3=B*V2
W23 = p2 * (V3-V2)
w23= cp_equ(T2) * (T3-T2)

# Power stroke: Isentropic Expansion
p4 = p3 * (V3/V1) ** k
T4 = T3 * (V3/V1) ** (k-1)
w34 = cv * (T3 - T4)
W34= m1 * R * (T4 - T3) / (1 - k)
# Blowdown
Qout = m1 * cv * (T1 - T4)

# Exhaust
p5 = p1
w56 = p5 * (Vd - Vc)

# Intake
w61 = p5*.8 * (Vc - Vd)
# Net work
W_net = m1 * (w12 + w34) + W23
W_net2 = (W12 + (W23+W34))

W_dot = C * W_net2 * N / 60 / 2  #  Go back and fix this!!!  ** change the cyclinders and if 2/4 
W_dotB=((2*np.pi*N) / 60) * (Torque/1000)   # Uses Peak torque-at peak horsepower and RPM at peak horsepower converted in KW
OPD= W_dotB/Vd
bmep=((4*np.pi*Torque) / (Vs*C)) / 1000 # Converts to KPa full displacement used and Peak torque at Peak horsepower
imep=(W_dot / W_dotB) * bmep
fmep=imep-bmep
# Display states:
print("p1:\t %.2f [kPa]\tT1:\t %.2f [K]" %(p1, T1))
print("p2:\t %.2f\tT2:\t %.2f" %(p2, T2))
print("p3:\t %.2f\tT3:\t %.2f" %(p3, T3))
print("p4:\t %.2f\t\tT4:\t %.2f" %(p4, T4))

print("Net work output:\t%.2f [kJ per cycle]" % W_net)
print("Net work output from book formulas:\t%.2f [kJ per cycle]" % W_net2)

print("Total engine power output (indicated):\t%.2f [kW]" % W_dot)
print("Total engine power output (Brake):\t%.2f [kW]" % W_dotB)
print("Output per Displacement:\t%.2f [kW/L]" % OPD)
print("Brake Mean Effective Pressure:\t%.2f [kPa]" % bmep)
print("Indicated Mean Effective Pressure:\t%.2f [kPa]" % imep)
print("Friction Mean Effective Pressure:\t%.2f [kPa]" % fmep)

#Set up volume vector and pressure vectors for p-V plot
V12 = np.linspace(V2, V1, 101)
V13 = np.linspace(V3, V1, 101)
p12 = p1 * np.divide(V12, V2) ** k
V23 = [V2, V3]
p23 = [p2, p3]
p34 = p3 * np.divide(V13, V3) ** -k
V45 = [V1, V1]
p45 = [p4, p1]

plt.close('all')
fig1 = plt.figure()
plt.plot(np.linspace(V1, V2, 101), p12, '-g')
plt.plot(V23, p23)
plt.plot(np.linspace(V3, V1, 101), p34, '-r')
plt.plot(V45, p45)


rc2=np.linspace(1,20,20)

ThermalD = np.empty(len(rc2))

for i in range(0,len(rc2),1):
    ThermalD[i]=(1-(((1/(i+1))**(k-1))*(((B**k)-1)/(k*(B-1)))))
    ThermalD[i]=ThermalD[i]*100



fig2=plt.figure()
plt.plot(rc2,ThermalD, '-g', label='Diesel')
plt.legend(loc='best')

for i in range(0,len(rc2),1):
    print("compression ratio:\t%.2f" % rc2[i], "Effciency: \t%.4f" %ThermalD[i])