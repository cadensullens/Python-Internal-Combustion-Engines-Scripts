# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:53:10 2020

@author: caden
"""

# -*- coding: utf-8 -*-"""Created on Fri Jan 25 07:47:37 2019
#******************************        Air standard OTTO        *******************************
#J.Rhett Mayor, Ph.D@author: jmayor3
import numpy as np
import math 
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
C=float(input("Number of Cylinders:\t"))    # Cyclinders
Vd=float(input("Engine Volume[L]:\t"))    # Total Displacement of engine in liters
Vs = Vd / C * 1e-3   #converts to m^3 this is swept Vs**
rc =float(input("What is your compression ratio? \t"))        #compression ratio**
Vc= Vs / (rc - 1.)    #Clearance Volume [m^3]
N=float(input("RPM Speed:\t"))     # Speed in RPM   # Peak Horsepower RPM

HP=200 # Given in Horsepower
Torque=HP*5252/ N    # Converts to lb-ft
Torque= (1/0.73756214927727)*Torque 
Torque=258.8  #converts to NewtonMeters

# Combustion and intake efficiencies
eta_c = 1
xr = 0.04#Exhaust residual percentage
eta_s = 1 - xr

# Fuel in this analysis is Gasoline
LHV = 44.3# MJ/kg    Gas is 44.3 light diesel is 42.5 and heavy is 41.4 Hydrogen is 120
psy = 1.2# Stoichioemetric ratio
AF_stoich = 14.7#Stoich AF ratio for Gas
AF = AF_stoich / psy

# intake state
p1 = 180#intake pressure [kPa]  if naturally aspirated 
T1 = 325. #converted to kelvin from celsius
V1 = Vc + Vs
cp = cp_equ(T1)
R = 0.287# Gas constant for air kJ/kg
cv = cp - R
#k = cp / cv
#]
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

# Instant. combustion
Qc = eta_c * m_f * LHV * 1e3
T3 = T2 + Qc / m1 / cv
p3 = p2 * (T3 / T2)

# Power stroke: Isentropic Expansion
p4 = p3 * (rc) ** -k
T4 = T3 * (rc) ** -(k-1)
w34 = cv * (T3 - T4)
W34= m1 * R * (T4 - T3) / (1 - k)
# Blowdown
Qout = m1 * cv * (T1 - T4)

# Exhaust
p5 = p1
w56 = p5 * (Vc - V1)

# Intake
w61 = p5 * (V1 - Vc)



# Net work
W_net = m1 * (w12 + w34)
W_net2 = (W12 + W34)
# =============================================================================
# W_dot2=W_net/(5.7/1000)            
# =============================================================================
W_dot = C * W_net * N / 60 / 2 #  Go back and fix this!!!  ** change the cyclinders and if stroke is 2=1 4=2 
W_dotB=((2*np.pi*N) / 60) * (Torque/1000)   # Uses Peak torque-at peak horsepower and RPM at peak horsepower converted in KW
OPD= W_dotB/Vd
bmep=((4*np.pi*Torque) / (Vs*C)) / 1000 # Converts to KPa full displacement used and Peak torque at Peak horsepower
imep=(W_dot / W_dotB) * bmep
fmep=imep-bmep
Mech=W_dotB/W_dot*100
massfuelrate=m_f*C*N/60*(1/(C/2))
bsfc=(massfuelrate/W_dotB) * 3595360.825
# Display states:
print("p1:\t %.2f [kPa]\tT1:\t %.2f [K]" %(p1, T1))
print("p2:\t %.2f\tT2:\t %.2f" %(p2, T2))
print("p3:\t %.2f\tT3:\t %.2f" %(p3, T3))
print("p4:\t %.2f\t\tT4:\t %.2f \n" %(p4, T4))

print("Net work output:\t%.2f [kJ per cycle]" % W_net)
print("Net work output from otto formulas:\t%.2f [kJ per cycle]" % W_net2)

print("Total engine power output (indicated):\t%.2f [kW]" % W_dot)
print("Total engine power output (Brake):\t%.2f [kW]" % W_dotB)
print("Output per Displacement:\t%.2f [kW/L]" % OPD)
print("Brake Mean Effective Pressure:\t%.2f [kPa]" % bmep)
print("Indicated Mean Effective Pressure:\t%.2f [kPa]" % imep)
print("Friction Mean Effective Pressure:\t%.2f [kPa]" % fmep)
print("Mechanical efficiency:\t%.2f " % Mech)
print("Brake Specific Fuel Consumption:\t%.2f [gm/Kwhr] \n" % bsfc)


#Set up volume vector and pressure vectors for p-V plot
V12 = np.linspace(V2, V1, 101)
p12 = p1 * np.divide(V12, V2) ** k
V23 = [V2, V2]
p23 = [p2, p3]
p34 = p3 * np.divide(V12, V2) ** -k
V45 = [V1, V1]
p45 = [p4, p1]

plt.close('all')
fig1 = plt.figure()
plt.plot(np.linspace(V1, V2, 101), p12, '-g')
plt.plot(V23, p23)
plt.plot(np.linspace(V2, V1, 101), p34, '-r')
plt.plot(V45, p45)

rc2=np.linspace(1,20,20)

Thermal = np.empty(len(rc2))
for i in range(0,len(rc2),1):
    Thermal[i]=(1-(1/(i+1))**(k-1))
    Thermal[i]=Thermal[i]*100


fig2=plt.figure()
plt.plot(rc2,Thermal, '-g', label='Otto')
plt.legend(loc='best')


for i in range(0,len(rc2),1):
    print("compression ratio:\t%.2f" % rc2[i], "Effciency: \t%.4f" %Thermal[i])





# =============================================================================
# Exhaust Valve
# =============================================================================
T3=float(input("Max Temp [K]?\t"))
p3=float(input("Max pressure [kpa]?\t"))
Bd=float(input("What is your Blowdown angle?\t"))
Bdrev=Bd/360
TimeBd=Bdrev / (N / 60)
Nv=float(input("Number of Intake valves:\t"))
S=input("Square engine Y/N?\t").upper()
if S == 'Y':
    B=(Vd * 1e-3 *(4/np.pi)) **(1/3) #Bore equals the stroke
    S=float(input("What is the Stroke or Bore in [m]? \t"))
    B=S
    Rr=float(input("What is the crank ratio? \t"))
else:
    a=float(input("Crank offset in [mm]?\t"))
    r=float(input("Rod length in [mm]?\t"))
    Rr=r/a
    S=2*a / 1000 # in {meters}
    B=math.sqrt(Vd/C * 1e-3 *(4/np.pi) / S)# in meters [m]

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
print("P7:\t\t %.3f[kPa]\t\t\tT7:\t\t %.3f[K] \nM7:\t %.8f[kg]\t\tMblow: %.8f[kg]\nMass gone in Blowdown:\t\t %.8f [percent] " %(P7, T7, M7, Mblow, Mper))
print("\n")
print("\t\t\t Exhaust Valve Sizing")
print("Blowdown Time:\t\t\t\t %.8f [s]" %TimeBd)
print("Speed of Sound:\t\t\t\t %.2f [m/s]" %c_ex)
print("Total intake valve size:\t %.8f [cm^2]" %Ai)
print("Individual valve size:\t\t %.3f [cm]" %Dv)
print("Upper Lift:\t\t\t\t\t %.2f [cm]" %l)



























