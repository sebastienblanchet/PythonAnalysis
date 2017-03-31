'''
Sebastien Blanchet,
Altaeros Energies,
Systems Engineering Intern,
March 27 2017

Note: all calculations to be completed in US customary
'''

# Import plotting modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm

# Define all parameters from thesis paper
t = 70 #mm
l = 1775 #mm
OD = 1000 #mm
R = OD/2 #mm
E_y = 2*10**5 #MPa
S_y = 520 # MPa
v = 0.35
zeta = 0.5
d_rope = 65 #mm
T = 1.96*10**3 #kN
P = (2*T)/d_rope  #MPa

# k term
w_p = ((P*R**2)/(E_y*t))
k = ((3*(1-v**2))**0.25)/((R*t)**0.5)

# Begin solving for t
w_tot = w_p

Sig_hoop = (E_y/R)*w_tot
Sig_hoop_approx = P*R/t

delta = abs((Sig_hoop-Sig_hoop_approx)/Sig_hoop_approx)*100

print("Accurate Sig= %.1f MPa , Approx. Sig= %.1f MPa " % (Sig_hoop, Sig_hoop_approx))

if delta >= 10:
    # %% is to print percentage
    print("Results not valid: %% %.1f error" %delta)

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