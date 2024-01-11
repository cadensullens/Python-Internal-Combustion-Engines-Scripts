# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 09:16:26 2020

@author: caden
"""

# -*- coding: utf-8 -*-"""Created on Tue Jan 22 20:05:15 2019@author: J.Rhett Mayor, Ph.D
#***************************************************************   IC Engine Code Example 2 *****************************************************
import numpy as np
import matplotlib.pyplot as plt
#    Establish the basic dimensions of the crack, rod 3-bar
a = 90.88/2.    #crank offset [mm]
R = 3.86#Rod length to crank ratio
r = R * a    # From (Eq 2-6 PBK)
N = 5600# RPM
# Basic cylinder geometry and clearance volume
Vc = 0.059  #[L]
B = 99.49
# Calcuatle the stroke and mean piston speed
S = 2. * a    # Engine stroke (= twice crank throw)
U_piston_mean = 2. * S * N / 1000 / 60

# Set up vector for {0, 180} crank angle sweep
theta = np.linspace(0, np.pi, 361)
R = np.linspace(2, 10, 5)
norm_U_piston = np.zeros((len(theta), len(R)))

omega = 2 * np.pi * N / 60.
dtheta = (theta[1] - theta[0])
dt = (theta[1] - theta[0]) / omega

#  Calculate piston position, speed & accel. (speed from numerical derivative)
s_piston = a * np.cos(theta) + np.sqrt(r **2. - a ** 2. * np.sin(theta)**2.)
s_dot_piston = 1e-3 * np.gradient(s_piston, dt)  # Inst. piston speed in [m/s]
s_ddot = (1/9.81) * np.gradient(s_dot_piston, dt) # inst. piston [g]
# Cylinder volume as a fuction of theta in [cc]
Vd_theta = 1e-3 * (Vc + np.pi / 4. * B ** 2. * (a + r - s_piston)) 
V_dot_intake = np.gradient(Vd_theta, dt) #required intake flow rate in cc/s3./6
m_dot_intake = 1.225e-3 * V_dot_intake

# Calculate the effect | impact of crank:rod ratio
for i in range(len(theta)):
    for j in range(len(R)):
        # Calculate the normalized instantaneous piston speed (2-5 PBK)        
        norm_U_piston[i,j] = (np.pi/2 * np.sin(theta[i]) *                  
                     (1 + (np.cos(theta[i]) / np.sqrt(R[j] ** 2. -                  
                           (np.sin(theta[i])**2)))))
        print("Dimensionless Piston Speed:\t%.4f" % norm_U_piston[i,j])
        # Work per cycle: W [kj/kg]
W_dot = 150# Total engine power in kW
W = W_dot * 2 / N / 60.
# =============================================================================
# Plotting following# 
#=============================================================================
plt.close('all')
fig1, ax1 = plt.subplots()
ax1.set_xlabel('Crank Angle [degrees]')
ax1.plot(np.degrees(theta), s_piston,'-r', label='Position')
ax2 = ax1.twinx()
ax2.plot(np.degrees(theta), s_dot_piston, '-g', label='Mean Piston Speed')
ax2.grid()
ax2.set_ylabel('Piston speed [m/s]')
ax1.set_ylabel('Piston position [mm]')
fig2 = plt.figure()
plt.xlim((0,180))
plt.grid()
for i in range(len(R)):    
    plt.plot(np.degrees(theta), norm_U_piston[:,i], label=("R %.0f" %R[i]))
    plt.legend(loc='best')
    plt.xlabel('Crank Angle [degrees]')
    plt.ylabel('Normalized Piston Speed')
    
fig3 = plt.figure()
plt.xlim((0,180))
plt.xlabel('Crank Angle [degrees]')
plt.ylabel('Cylinder volume [cc]')
plt.plot(np.degrees(theta), Vd_theta, '-b')
fig4, ax3 = plt.subplots()
ax3.set_xlabel('Crank angle')
ax3.plot(np.degrees(theta),m_dot_intake, '-b', label='Intake flow rate [cc/s]')
ax3.set_ylabel('Intake flow rate [cc/s]')
ax4 = ax3.twinx()
ax4.plot(np.degrees(theta), Vd_theta, '-r', label='Cylinder Volume [cc]')
ax3.grid()
ax3.legend(loc='upper right')
ax4.legend(loc='lower right')