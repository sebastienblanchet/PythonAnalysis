# Sebastien Blanchet 
# Altaeros Energies, Systems Intern
# March 29 2017

# Import plotting modules
import matplotlib.pyplot as plt
import numpy as np

# Define Parameters
d_tether = 15.3/25.4
L_tether = (410*1000)/25.4
P_metric = 51.264   # force P kN
SF = 1.5          # safety factor
S_y = 36259.4       # Yield strength 250 MPa
E_y = 2.9*10**7   # Youngs modulus 200GPa
v = 0.3           # poissons ratio
rho = 0.284     # density of steel lb/in^3
# D_min = 10
# D_max = 50
P = 224.81*P_metric

# Length Calculation for Drum Length
# len_vect = abs(D_max-D_min)+1
D = 22

D_p1 = D+d_tether
D_p2 = D_p1+d_tether
L_p1 = np.pi*D_p1
L_p2 = np.pi*D_p2

n_wraps = L_tether/(L_p1+L_p2)
SG = 1.01
l = SG*n_wraps*d_tether
l_ws = SG*d_tether
mu = 0.05
wraps = np.arange(0, 30, 0.1)

T_cap = P*np.exp(-mu*2*np.pi*wraps)

q = 2*T_cap/(D*d_tether)
x_0 = SG * d_tether * wraps

# Plot 1 : and l vs t
plt.figure(1)
plt.plot(x_0, q, label='Pressure q')
plt.legend()
plt.title('Capstan Equation')
plt.xlabel('Distance from applied force [in]')
plt.ylabel('Pressure q [psi]')
plt.savefig('Figures\Cap.png')
plt.show()

Q_Exp = np.asarray([x_0, q])
np.savetxt('Data\DistLoadCap.csv', np.transpose(Q_Exp), delimiter=",")







