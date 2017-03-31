# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary
# Thin shell analysis of 28" drum

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
D = 28                      # OD 28 in drum for thin shell
# Prelim calcs
Sig_allow = S_y / SF
l = np.ceil((SG * d_tether * L_tether) / (np.pi * (2 * D + 3 * d_tether)))

Bx = np.arange(0, 11, 0.1)

Phi_Bx = (np.exp(-Bx))*(np.cos(Bx)+np.sin(Bx))
Psi_Bx = (np.exp(-Bx))*(np.cos(Bx)-np.sin(Bx))
Theta_Bx = (np.exp(-Bx))*(np.cos(Bx))
Zeta_Bx = (np.exp(-Bx))*(np.sin(Bx))

plt.figure(1)
# to write math text must use r'$\greekletter text$'
plt.plot(Bx, Phi_Bx, label=r'$\phi ( \beta x)$')
plt.plot(Bx, Psi_Bx, label=r'$\psi ( \beta x)$')
plt.plot(Bx, Theta_Bx, label=r'$\theta ( \beta x)$')
plt.plot(Bx, Zeta_Bx, label=r'$\zeta ( \beta x)$')
plt.xlabel(r'$\beta x$')
plt.ylabel('Bending fxn')
plt.legend()
plt.show()
