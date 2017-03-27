# Import plotting modules
import matplotlib.pyplot as plt
import numpy as np

# Define Parameters
d_tether = 15.3/25.4
L_tether = (410*1000)/25.4
P_max_metric = 51.264   # force P kN
P_op_metric = 26
SF = 1.5          # safety factor
S_y = 36260       # Yield strength 250 MPa
E_y = 2.9*10**7   # Youngs modulus 200GPa
v = 0.3           # poissons ratio
rho = 0.284     # density of steel lb/in^3

# Length Calculation for Drum Length
D_min = 20
D_max = 50
D = np.linspace(D_min, D_max, num=(abs(D_max-D_min)+1))
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
P = 224.81 * P_max_metric
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

T_max = P_max_metric * (((D_p2 / 2) * 25.4) / 1000)
n_max = (60/(2*np.pi))*((2*(1))/(25.4*D_p1/1000))

T_cont = (26)*(((D_p2/2)*25.4)/1000)
n_cont = (60/(2*np.pi))*((2*(0.5))/(25.4*D_p1/1000))


# Plotting operational points
plt.figure(1)
plt.title('Gearbox Output vs Drum Diameter D')
plt.plot(D, T_max, 'b', label='max. T')
plt.plot(D, T_cont, 'b--', label='Cont. T')
plt.ylabel('Torque T [kN]')
plt.xlabel('Drum OD [in]')
plt.legend()
plt.savefig('Figures/GB_torque.png')

# Showcasing Worst Case
plt.figure(2)
plt.title('Gearbox Output vs Drum Diameter D')
plt.plot(D, n_max, 'b-o', label='max. n')
plt.plot(D, n_cont, 'b--', label='Cont. n')
plt.ylabel('Speed n [kN]')
plt.xlabel('Drum OD [in]')

# for xy in zip(D[::3], n_max[::3]):
#     plt.annotate('%s, %.0f RPM' % xy, xy=xy, textcoords='data')

plt.legend()
plt.savefig('Figures/GB_speed.png')

# export to excel
exportTorque = np.asarray([D, T_max, T_cont])
exportspeed = np.asarray([D, n_max, n_cont])
np.savetxt("Data\Torquedata.csv", np.transpose(exportTorque), delimiter=",")
np.savetxt("Data\speeddata.csv", np.transpose(exportspeed), delimiter=",")


