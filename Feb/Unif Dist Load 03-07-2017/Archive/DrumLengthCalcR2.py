# Import plotting modules
import matplotlib.pyplot as plt
import numpy as np

# Define Parameters
d_tether = 15.3/25.4
L_tether = (410*1000)/25.4
P_metric = 51.264   # force P kN
SF = 1.5          # safety factor
S_y = 36260       # Yield strength 250 MPa
E_y = 2.9*10**7   # Youngs modulus 200GPa
v = 0.3           # poissons ratio
rho = 0.284     # density of steel lb/in^3

# Length Calculation for Drum Length
D = np.linspace(20, 40, num=21)
v_hi = np.linspace(0, 1, num=11)

D_p1 = D+d_tether
D_p2 = D_p1+d_tether
L_p1 = np.pi*D_p1
L_p2 = np.pi*D_p2

n_wraps = L_tether/(L_p1+L_p2)
SG = 1.01
l = SG*n_wraps*d_tether
l_ws = SG*d_tether

x_lw = 1000*(l_ws/L_p1)

# Calculated Parameters
R = D/2
P = 224.81*P_metric
Sig_all = S_y/SF
q = 100*P/(D**2)
lmb = ((3*(1-v**2))/(R**2))**0.25
Dm = E_y/(12*(1-v**2))

# Displacement Case 28-1c
DelR_c = (q*(R**2)/E_y)*(1-v/2)

# Displacement Case 28-8
y_a8 = -1/(2*Dm*lmb**3)
phi_a8 = 1/(2*Dm*lmb**2)

# Displacement Case 28-10
y_a10 = 1/(2*Dm*lmb**2)
phi_a10 = -1/(Dm*lmb)


M_o = (((-phi_a8*y_a10/y_a8)+phi_a10)**-1)*((phi_a8*DelR_c)/y_a8)
V_o = (DelR_c/-y_a8)+(y_a10*M_o/-y_a8)

# Stress Calculations Case 28-1c
Sig_1c = q*R/2
Sig_2c = q*R

# Stress Calculations Case 28-10
Sig_28 = -2*V_o*lmb*R

# Stress Calculations Case 28-10
Sig_110pr = 6*M_o
Sig_210 = 2*M_o*(lmb**2)*R
Sig_210pr = v*Sig_110pr

Sm_Sig1 = Sig_1c+Sig_110pr
t = Sm_Sig1/Sig_all

D_i = D-(2*t)
A = (np.pi/4)*((D**2)-(D_i**2))
mass = rho*A*l

T_max = P_metric*(((D_p2/2)*25.4)/1000)

# Stress calcs for thickness Calculated Parameters
# Plot thickness results worst case
# FigSave = "Drum D,L vs t " + time.strftime("%x") + " " + time.strftime("%X")

plt.figure(1)
plt.plot(t, D, label='Outer Diameter D')
plt.plot(t, l, label='Length l')
plt.legend()
plt.title('Drum Sizing Parameters')
plt.ylabel('Unit Length [in]')
plt.xlabel('Drum Thickness [in]')
plt.savefig('Figures/Drum_D_v_t', dpi='400', format='pdf')
plt.close()
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom='off',      # ticks along the bottom edge are off
#     top='off',         # ticks along the top edge are off
#     labelbottom='off') # labels along the bottom edge are off

plt.figure(2)
plt.plot(D, mass)
plt.title('Drum Mass vs OD')
plt.ylabel('Mass [lbs]')
plt.xlabel('Drum OD [in]')
plt.savefig('Figures/Drum_mass_v_D', dpi='400', format='pdf')

plt.figure(3)
plt.plot(D, T_max)
plt.title('Gearbox Max Torque')
plt.ylabel('Torque output [kN]')
plt.xlabel('Drum OD [in]')
plt.savefig('Figures/Drum_Tmax_v_t', dpi='400', format='pdf')
plt.close()
