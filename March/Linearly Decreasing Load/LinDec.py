'''
Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
March 27 2017

Note: all calculations to be completed in US customary
'''

# Import plotting modules
import numpy as np
import matplotlib.pyplot as plt
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
D = 28
t_min = 0.5
t_max = 2
t_step = 1/16

#calculate length for given D
l = np.ceil((SG*d_tether*L_tether)/(np.pi*(2*D+3*d_tether)))

# create non dimensional vector zeta and thickness vector
zeta = np.arange(0, 1.1, 0.1)
t = np.arange(t_min, t_max+t_step, t_step)

# Solution matrix
ZETA, T = np.meshgrid(zeta,t)
# Row and col nums
numcols = len(zeta)
numrows = len(t)

Sig_h = np.zeros()

for i in range(numrows):
    for j in range(numcols):

# # 3D Plot
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(D_x, t_y, Valid_DT, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# # Customize the z axis.
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax.set_xlabel('Diameter OD [in.]')
# ax.set_ylabel('Thickness t [in.]')
# ax.set_zlabel('Lmbd*l')
# plt.title('3D plot')
# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.show()