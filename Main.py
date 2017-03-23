# Import plotting modules
import numpy as np
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

'''
Note: all calculations to be completed in US customary
Sebastien Blanchet
03/09/2017
'''

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
D = np.linspace(D_min, D_max, num=(len_vect))

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

numrows = len(D_x)
numcols = len(D_x[0])

shapeofD = np.shape(D_x)

Valid_DT = np.zeros(shapeofD)
Valid_longshell = np.zeros(shapeofD)

for i in range(numrows):

    for j in range(numcols):

        # Loop to check for assumptions
        if TWPV[i][j] >= 10:
            Valid_DT[i][j] = 1

        if longshell[i][j] >= 6:
            Valid_longshell[i][j] = 1


# # Displacement Case 28-1c
# DelR_c = (q*(R**2)/E_y)*(1-v/2)
# # Displacement Case 28-8
# y_a8 = -1/(2*Dm*lmb**3)
# phi_a8 = 1/(2*Dm*lmb**2)
# # Displacement Case 28-10
# y_a10 = 1/(2*Dm*lmb**2)
# phi_a10 = -1/(Dm*lmb)
# # Reaction momment and shear
# M_o = (((-phi_a8*y_a10/y_a8)+phi_a10)**-1)*((phi_a8*DelR_c)/y_a8)
# V_o = (DelR_c/-y_a8)+(y_a10*M_o/-y_a8)
# # Stress Calculations Case 28-1c
# Sig_1c = q*R/2
# Sig_2c = q*R
# # Stress Calculations Case 28-10
# Sig_28 = -2*V_o*lmb*R
# # Stress Calculations Case 28-10
# Sig_110pr = 6*M_o
# Sig_210 = 2*M_o*(lmb**2)*R
# Sig_210pr = v*Sig_110pr
# # Final stress calc
# Sm_Sig1 = Sig_1c+Sig_110pr


# # Intialize loop values
# max_itnum = 5
# itnum = []
# tol = 10  #10% convergence error
# err = []
# Iterative loop
# for i in range(max_itnum):
#     itnum[i] = i
#
#     # err[i] = max((abs(new-old)/old)*100)
#     if err[i] <= tol:
#         break

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

# Convergence plot
# plt.figure(2)
# plt.plot(itnum, err, label='')
