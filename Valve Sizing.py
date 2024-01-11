# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 19:58:26 2020

@author: caden
"""

import numpy as np
import math 


# =============================================================================
# Intake Valve sizing and timing
# =============================================================================
N=float(input("RPM Speed:\t"))     # Speed in RPM
c=float(input("Number of Cylinders:\t"))    # Cyclinders
Vd=float(input("Engine Volume[L]:\t"))   # Volume in liters
T=float(input("Inlet Temperature [K]:\t"))    # Intake temp
Nv=float(input("Number of Intake valves:\t"))
S=input("Square engine Y/N?\t").upper()

if S == 'Y':
# =============================================================================
#     B=(Vd * 1e-3 *(4/np.pi)) **(1/3) #Bore equals the stroke
#     S=B
# =============================================================================
    S=float(input("What is the Stroke or Bore? \t"))
    B=S
else:
    a=float(input("Crank offset in [mm]?\t"))
    r=float(input("Rod length in [mm]?\t"))
    S=2*a / 1000 # in {meters}
    B=math.sqrt(Vd/c * 1e-3 *(4/np.pi) / S)# in meters [m]
C=input("Speed of sound given Y/N?\t").upper()
if C =='Y':
    C=float(input("Speed of Sound given [m/s]:\t"))
else:
    C=math.sqrt(1.4*287*(T))  # [speed of sound m/s]


Up=2 * S *N/60  #mean piston speed [m/s]
Ai=1.3* B**2 * (Up/C) * 1e4   #Total intake valve area [cm^2]
Dv=math.sqrt(Ai  / Nv * (4 / np.pi))   # diameter of each valve in [cm]
l= Dv / 4     # upper limit to valve lift [cm]

print("\n")
print("\t\t\t Valve Sizing")
print("Speed of Sound:\t\t\t\t %.2f [m/s]" %C)
print("Total intake valve size:\t %.8f [cm^2]" %Ai)
print("Individual valve size:\t\t %.3f [cm]" %Dv)
print("Upper Lift:\t\t\t\t\t %.2f [cm]" %l)


# =============================================================================
# Intake Valve timing
# =============================================================================
N=float(input("RPM Speed:\t"))     # Speed in RPM
Op=float(input("Open Angle:\t"))    # Open angle bTDC
Cl=float(input("Close Angle:\t"))   #Close angle aBDC
l=float(input("Max lift [mm]:\t"))  # Max valve lift
Delta=Op + Cl + 180     # Total theta with the dwell(180)
Deltat=Delta/((N/60)*360)*1000  # Valve timing [total] [ms]
Full=float(input("Fully open time [%]:\t")) # percentage 
E=input("Are the open and close time equal Y/N?\t").upper()
if E=='Y':
    Percentage= (100-Full) / 2
    DeltaTime=Percentage *Deltat /100     # Time opening and closing
    Disp=0.5*l / 1000   # Displacement [m]
    A=(2*Disp) / ((DeltaTime/1000/2) ** 2 )     # Acceleration [m/s^2]
else:
     PercentageO=float(input("Percentage opening[%]:\t"))
     PercentageC=float(input("Percentage closing[%]:\t"))
     DeltaTimeO=PercentageO *Deltat  /100    # Time opening {ms}
     DeltaTimeC=PercentageC *Deltat /100     # Time closing {ms}
     Disp=0.5*l / 1000   # Displacement [m]
     Ao=(2*Disp) / ((DeltaTimeO/1000/2) ** 2 )      # Acceleration for opening [m/s^2]
     Ac=(2*Disp) / ((DeltaTimeC/1000/2) ** 2 )      # Acceleration for closing [m/s^2]
MaxA=float(input("Max Accel. [g]:\t"))
Mass=float(input("Mass[kg]:\t"))
Force= MaxA*9.81*Mass
if E=='Y':
     print("\n")
     print("\t\t\t Valve Timing")
     print("Time valves are open (partially/fully): %.4f [ms]" %Deltat)
     print("Acceleration to Open/Close: %.2f [m/s^2]" %A)
     print("Force to close Valve:\t %.2f [N]" %Force)
else:
    print("\n")
    print("\t\t\t Valve Timing")
    print("Time valves are open(partially/fully): %.4f [ms]" %Deltat)
    print("Time valves are open (fully): %.2f [ms]" %DeltaTimeO)
    print("Time valves are closed (fully): %.2f [ms]" %DeltaTimeC)
    print("Acceleration to Open: %.2f [m/s^2]" %Ao)
    print("Acceleration to Close: %.2f [m/s^2]" %Ac)
    print("Force to close Valve:\t %.2f [N]" %Force)

    
    