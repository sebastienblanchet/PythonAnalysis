'''
Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
March 27 2017

Note: all calculations to be completed in US customary
'''

# Import plotting modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker

# Define all parameters
d_tether = 15.2/25.4        # tether diameter 15.2mm to in
L_tether = (370*1000)/25.4  # tether length overall 370m to in
P = 224.809*51.264          # force P kN
SF = 1.5                    # safety factor
S_y = 50763                 # Yield strength 350 MPa
E_y = 2.9*10**7             # Youngs modulus 200GPa
v = 0.3                     # poissons ratio
rho = 0.284                 # density of steel lb/in^3
SG = 1.01                   # spool gap 1%

# Vector limits for exploration
D_min = 10
D_max = 50
t_min = 0.5
t_max = 2
t_step = 1/16

# Length Calculation for Drum Length
len_vect = abs(D_max-D_min)+1
D = np.linspace(D_min, D_max, num=len_vect)

# Create vector of from [t_min, t_max] (must add t_step arrange does NOT include final)
t = np.arange(t_min, t_max+t_step, t_step)

# Create 3D mesh of diameter and thickness
D_x, t_y = np.meshgrid(D, t)

# Drum and tether relations
D_p1 = D_x + d_tether
D_p2 = D_p1+d_tether
L_p1 = np.pi*D_p1
L_p2 = np.pi*D_p2
n_wraps = L_tether/(L_p1+L_p2)
l_ws = SG*d_tether
l = n_wraps*l_ws

# Drum mass properties
R = D_x / 2
D_i = D_x - (2 * t_y)
A = (np.pi/4)*((D_x ** 2) - (D_i ** 2))
mass = rho*A*l
I_m = (0.5*mass*((R**2)+((D_i/2)**2)))
q = (2*P)/(D_x *d_tether)
lmb = ((3*(1-v**2))/((R**2)*(t_y**2)))**0.25
Dm = (E_y*(t_y**3))/(12*(1-v**2))
longshell = lmb*l
TWPV = R/t_y

#Vector Prop
numrows = len(D_x)
numcols = len(D_x[0])
shapeofD = np.shape(D_x)

Valid_DT = np.zeros(shapeofD)
Valid_longshell = np.zeros(shapeofD)

# Loop to check for assumptions
for i in range(numrows):
    for j in range(numcols):
        if TWPV[i][j] >= 10:
            Valid_DT[i][j] = 1
        if longshell[i][j] >= 6:
            Valid_longshell[i][j] = 1

# Plot 1 mass vs D vs t
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(D_x, t_y, Valid_DT, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# Customize the z axis.
ax.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.set_xlabel('Diameter OD [in.]')
ax.set_ylabel('Thickness t [in.]')
ax.set_zlabel('Lmbd*l')
plt.title('3D plot')
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()