#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:13:48 2020

@author: JeremiahDavid
"""

import numpy as np
import math
import matplotlib.pyplot as plt

def cp_eq(T):
    theta = T/100
    cp0_N2 = 39.060 - 512.79*(theta)**-1.5 + 1072.7*(theta)**-2 - 820.40*(theta)**-3
    cp0_O2 = 37.434 + 0.0201*(theta)**1.5 - 178.57*(theta)**-1.5 + 236.88*(theta)**-2.0
    
    cpq = (0.79/28)*cp0_N2 + (0.21/32)*cp0_O2
    
    return cpq

# =============================================================================
# Main Program (values from example 1 in chapter 3)
# =============================================================================
    
# =============================================================================
# Variables
# =============================================================================

nc = 6 # number of cylinders
rc = 11.5 # compression ratio
N = 1850 # RPM
n = 2 # rev/cycle (4 stroke)
bore = .086
stroke = .086
turbo = 0 # [psi] forced induction amount from turbo charger

A_p = bore**2*(np.pi/4) # piston head area

Vs = A_p*stroke # swept volume for a cylinder [m^3]
Vc = (Vs/(rc-1)) # clearance volume for a cylinder

Vd_L = 3 # [L] total engine displacement
Vd_m3 = Vd_L/1000 # [m^3] total engine displacement 

eta_c = .95 # combustion efficiency (usually this is between .95 and .98)
xr = 0.04 # exhaust residual percentage
eta_s = 1-xr # scavenging effeciency 

# fuel analysis values
LHV = 120000 #[kJ/kg]
psy = 0.5 # equivilance ratio (<1: lean)(>1: rich)
AF_stoich = 34
AF = AF_stoich/psy
k = 1.35  # Ratio of specific heat
R = 0.287 # Gas constant for air [kJ/kg-K]
rho_a = 1.181 # [kg/m^3] air standard density of air

# Initial states
p_atm = 101.3
boost = turbo*6.895 # [kPa] pressure boost from the turbo charger
T_atm = 300 # [K]

# intake state
p1 = p_atm + boost # intake pressure [kPa]
T1 = T_atm*(p1/p_atm) # [K]
V1 = Vc + Vs


# =============================================================================
# Get mass values
# =============================================================================
m1 = p1*V1/(R*T1) # Ideal gas law
m_m = m1 # mass of the mixture
m_a = eta_s*m_m*(AF/(1+AF)) # [kg/cylinder] mass of the air per cylinder
m_f = eta_s*m_m*(1/(1+AF)) # [kg/cylinder] mass of the air per cylinder


# =============================================================================
# Start analysis
# =============================================================================

# Process 1-2: Compression Stroke (Isentropic Compression)
p2 = p1*rc**k
T2 = T1*rc**(k-1)
V2 = Vc

W12 = m_m*R*(T2-T1)/(1-k) # work during compression

# specific heat constants
cp = cp_eq(T1) # pressure
cv = cp - R # volume

# Process 2-3: Combustion (Contant Volume)
Qc = m_f*LHV*eta_c # heat during combustion
T3 = T2 + Qc/(m_m*cv)
p3 = p2 * (T3/T2)
V3 = V2

# Process 3-4: Power Stroke (Isentropic Expansion)
p4 = p3*(rc)**-k
T4 = T3*(rc)**-(k-1)
V4 = V1

W34 = m_m*R*(T4-T3)/(1-k) # Work during expansion


# Blowdown
Qout = m_m*cv*(T1-T4)

# Exhaust
p5 = p1
V5 = Vs
W56 = p5*(Vs-Vc)

# Intake
p6 = p1
V6 = Vc
W61 = p5*0.8*(Vc-Vs)

# =============================================================================
# Final metrics
# =============================================================================

# Net work (indicated)
W_indicated = (W12 + W34) # [kJ] indicated work per cylinder
W_indicated_total = W_indicated*nc #[kJ] indicated work for total engine

# Engine power (indicated)
Power_indicated = (W_indicated*N)/(60*n) #[kW/cylinder]
Power_indicated_total = Power_indicated*nc #[kW] for the whole engine

# Thermal efficiency
eta_t = 1-(1/rc)**(k-1)

# Volumetric efficiency
eta_v = m_a/(rho_a*(Vs+Vc))

# Mean piston speed
u_p = 2*stroke*(N/60) # [m/s]

# Mechanincal efficiency (using figure 11 pg. 54)
eta_m = float(input('enter the mechanical efficiency associated with a piston speed of %.2f: ' %(u_p)))

# mean effective pressure
imep = W_indicated_total/(Vd_m3) # indicated
bmep = eta_m*imep # brake
fmep = imep - bmep # friction

# Net work (brake)
W_brake_total = W_indicated_total*eta_m

# Net work (brake) from geometries
W_brake_total2 = (bmep*Vd_m3)

# Engine power (brake)
Power_brake_total = Power_indicated_total*eta_m

# Engine power (brake) from geometries
Power_brake_total2 = bmep*(A_p*nc)*u_p/(2*n)

# Torque
torque1 = 60*(Power_brake_total)/(2*np.pi*N) # [N-m] (eq.43)
torque2 = (bmep*Vd_m3)/(4*np.pi) # [N-m] (eq.38)

# test brake work
### from eq.38 **this equation is for one revolution** analysis is for 2 ( so multiply by n)
W_brake_total3 = n*2*np.pi*torque1
W_brake_total4 = n*2*np.pi*torque2

# Brake specific fuel consumption
m_f_rate = m_f*(N/n)*(1/60) # [kg/s]
m_a_rate = m_a*(N/n)*(1/60)*nc # [kg/s]
isfc = 3600*1000*m_f_rate*nc/Power_indicated_total # [g/kW-hr] for total engine
bsfc = isfc/eta_m # [g/kW-hr] for total engine

# Output per displacement
opd = Power_brake_total/(Vd_L) # [kW/L]

# Power lost to friction
Power_friction_total = fmep*(A_p*nc)*u_p/(2*n) # [kW]

# Specific Power
specific_power = Power_brake_total/(A_p*nc)

print("\n")
print("\tEngine analysis for an Engine of the following parameters:")
print("RPM:\t\t\t %.2f [rev per min]"%N)
print("Cylinders:\t\t %.2f"%nc)
print("Displacement:\t\t %.2f [L]"%Vd_L)
print("Strokes:\t\t %.2f"%(n*2))
if turbo != 0:
    print("Forced induction:\t %.2f psi boost"%turbo)
else:
    print("Natural induction:\t %.2f psi boost"%turbo)

print("\n")
print("\tState Table")
print("p1: %.2f [kpa]\t\tT1: %.2f [K]" %(p1,T1))
print("p2: %.2f [kpa]\t\tT2: %.2f [K]" %(p2,T2))
print("p3: %.2f [kpa]\t\tT3: %.2f [K]" %(p3,T3))
print("p4: %.2f [kpa]\t\tT4: %.2f [K]" %(p4,T4))

print("\n")
print("\tWork and power Metrics")
print("Indicated work output:\t\t\t%.2f [kJ per cycle]" %(W_indicated_total))
print("Indicated engine power output:\t\t%.2f [kW]" %(Power_indicated_total))
print("Brake work output:\t\t\t%.2f [kJ per cycle]" %(W_brake_total))
print("Brake engine power output:\t\t%.2f [kW]" %(Power_brake_total))

print("\n")
print("\tOther Interesting Metrics")
print("Brake specific fuel consumption:\t %.2f [g/kW-hr]" %(bsfc))
print("Output per displacement:\t\t %.2f [kW/L]"%opd)
print("Specifc Power:\t\t\t\t %.2f [kW/m^2]"%specific_power)
print("Thermal efficiency:\t\t\t %.2f"%eta_t)
print("Volumetric efficiency:\t\t\t %.2f"%eta_v)
print("Power lost to friction:\t\t\t %.2f [kW]"%Power_friction_total)


p = lambda v,T: m_m*R*T/v # function to get pressure
V12 = np.linspace(V1,V2,101) # range: [1-rc]
T12 = np.linspace(T1,T2,101)
p12 = p(V12,T12)
# p12 = p1*np.divide(V12,V2)**k # a vector of (p2 = p1*rc^k) 

V23 = [V2, V2] # constant volume 
p23 = [p2, p3] 

# p34 = p3*np.divide(V12,V2)**-k # a vector of (p2 = p1*rc^-k) 
V34 = np.linspace(V3,V4,101) # range: [1-rc]
T34 = np.linspace(T3,T4,101)
p34 = p(V34,T34)

V45 = [V1, V1] # constant volume 
p45 = [p4, p5]

plt.close('all')
fig1 = plt.figure()
plt.plot(np.linspace(V1,V2,101),p12,'-b')
plt.plot(V23,p23,1000,'-g')
plt.plot(np.linspace(V2,V1,101),p34,'-b')
plt.plot(V45,p45,'-g')


# =============================================================================
# Exhaust Valve
# =============================================================================
T3=float(input("Max Temp [K]?\t"))
p3=float(input("Max pressure [kpa]?\t"))
Bd=float(input("What is your Blowdown angle?\t"))
Bdrev=Bd/360
TimeBd=Bdrev / (N / 60)
Nv=float(input("Number of Exhaust valves:\t"))
S=input("Square engine Y/N?\t").upper()
if S == 'Y':
    B=(Vd_L * 1e-3 *(4/np.pi)) **(1/3) #Bore equals the stroke
    S=float(input("What is the Stroke or Bore in [m]? \t"))
    B=S
    Rr=float(input("What is the crank ratio? \t"))
else:
    a=float(input("Crank offset in [mm]?\t"))
    r=float(input("Rod length in [mm]?\t"))
    Rr=r/a
    S=2*a / 1000 # in {meters}
    B=math.sqrt(Vd_L/nc * 1e-3 *(4/np.pi) / S)# in meters [m]

Cos=math.cos(math.radians(180-Bd))
Sin=math.sin(math.radians(180-Bd))**2
Vevo=Vc * (1 + ((0.5 * (rc - 1))*(Rr + 1 - Cos - math.sqrt(Rr**2 - Sin))))
Tevo=T3 * (Vc/Vevo)**(k - 1)         #[K]
Pevo=p3 * (Vc/Vevo)**k         #[kpa]
Mevo=(Pevo * Vevo) / (R * Tevo)       # {kg}
P7=p1
T7=T3*(P7/p3)**((k-1)/k)
M7=(P7*Vs) / (R * T7)
Mblow=Mevo-M7
Mper=(Mblow/Mevo)*100

C=input("Speed of sound given Y/N?\t").upper()
if C =='Y':
    C=float(input("Speed of Sound given [m/s]:\t"))
else:
    Tevo=925
    c_ex= math.sqrt(k*287*Tevo)         #Exhaust/blowdown velocity[m/s]

Up=2 * S *N/60  #mean piston speed [m/s]
Ai=1.3* B**2 * (Up/c_ex) * 1e4   #Total intake valve area [cm^2]
Dv=math.sqrt(Ai  / Nv * (4 / np.pi))   # diameter of each valve in [cm]
l= Dv / 4     # upper limit to valve lift [cm]
print("\n Exhaust Pressure,Temperature, and Volume")
print("Pevo:\t %.3f[kPa]\t\tTevo:\t %.3f[K] \nVevo:\t %.8f[m^3]\tMevo: %.8f[kg]" %(Pevo, Tevo, Vevo, Mevo))
print("P7:\t\t %.3f[kPa]\t\t\tT7:\t\t %.3f[K] \nMevo:\t\t %.8f[kg]\t\tMblow: %.8f[kg]\nMass gone in Blowdown:\t\t %.8f [percent] " %(P7, T7, Mevo, Mblow, Mper))
print("\n")
print("\t\t\t Exhaust Valve Sizing")
print("Blowdown Time:\t\t\t\t %.8f [s]" %TimeBd)
print("Speed of Sound:\t\t\t\t %.2f [m/s]" %c_ex)
print("Total intake valve size:\t %.8f [cm^2]" %Ai)
print("Individual valve size:\t\t %.3f [cm]" %Dv)
print("Upper Lift:\t\t\t\t\t %.2f [cm]" %l)








