# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary

# Import plotting modules
import numpy as np
import matplotlib.pyplot as plt

# Define all parameters
d_tether = 15.2/25.4        # tether diameter 15.2mm to [in]
L_tether = (370*1000)/25.4  # tether length overall 370m to [in]
P = 224.809*51.264          # force P of 51.264 kN [lbs]
SF = 1.5                    # safety factor [ul]
S_y = 50763                 # Yield strength 350 [MPa]
E_y = 2.9*10**7             # Youngs modulus 200 [GPa]
v = 0.3                     # poissons ratio [ul]
rho = 0.284                 # density of steel [lb/in^3]
SG = 1.01                   # spool gap 1% [ul]

# Prelim calcs
Sig_allow = S_y / SF

# Vector limits for exploration [in] for all
D_min = 10
D_max = 50
D_step = 1
t_min = 0.5
t_max = 2
t_step = 1/16

# Length Calculation for Drum Length
D = np.arange(D_min, D_max+D_step, D_step)

# Create vector of from [t_min, t_max] (must add t_step arrange does NOT include final) t=t_min:t_step:t_max
# t = np.arange(t_min, t_max+t_step, t_step)

# Length calculation [in]
l = (SG * d_tether * L_tether) / (np.pi * (2 * D + 3 * d_tether))

# Pressure
q = (2*P)/(D*d_tether)

# Aspect ratio
AR = l/D

# Plot aspect ratio vs D
plt.figure(1)
plt.plot(D, AR, label='Aspect Ratio')
plt.title('Aspect Ratio vs Outer Diameter')
plt.xlabel('Outer Diameter D [in]')
plt.ylabel('l/D')
plt.savefig('Figures\AR_vs_OD.png')
plt.show()