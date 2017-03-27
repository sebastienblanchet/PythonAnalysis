# Import plotting modules
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# Define Parameters
d_tether = 15.3/25.4
L_tether = (410*1000)/25.4
P_metric = 30  # force P kN
SF = 1.5          # safety factor
S_y = 50763      # Yield strength 350 MPa
E_y = 2.9*10**7   # Youngs modulus 200GPa
v = 0.3           # poissons ratio
rho = 0.284     # density of steel lb/in^3
D_min = 10
D_max = 50
SG = 1.01   # spool gap 1%
P = 224.809*P_metric

# Length Calculation for Drum Length
len_vect = abs(D_max-D_min)+1
D = np.linspace(D_min, D_max, num=(len_vect))

# Drum and tether relations
D_p1 = D+d_tether
D_p2 = D_p1+d_tether
L_p1 = np.pi*D_p1
L_p2 = np.pi*D_p2
n_wraps = L_tether/(L_p1+L_p2)
l_ws = SG*d_tether
l = n_wraps*l_ws

# Calculated Parameters
R = D/2
Sig_all = S_y/SF
q = (2*P)/(D*d_tether)
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
# Reaction momment and shear
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
# Final stress calc
Sm_Sig1 = Sig_1c+Sig_110pr
t = Sm_Sig1/Sig_all

# Drum mass properties
D_i = D-(2*t)
A = (np.pi/4)*((D**2)-(D_i**2))
mass = rho*A*l
I_m = (0.5*mass*((R**2)+((D_i/2)**2)))*0.000292641

# Gearbox output torques and speeds
T_max = (51.264)*(((D_p2/2)*25.4)/1000)
n_max = (60/(2*np.pi))*((2*(1))/(25.4*D_p1/1000))
T_cont = (26)*(((D_p2/2)*25.4)/1000)
n_cont = (60/(2*np.pi))*((2*(0.5))/(25.4*D_p1/1000))


# Plot 1 : and l vs t
fig0, ax1 = plt.subplots()
plt.title('Drum Mass Properties')
ax1.plot(D, l, 'b', label='Length l')
ax1.set_ylabel('Length l [in]')
ax1.set_xlabel('Drum OD [in]')
ax2 = ax1.twinx()
ax2.plot(D, t, 'g', label='Thickness t')
ax2.set_ylabel('Thickness t [in]')
# Get legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)
# Format secondary thickness axis
y_labels = ax2.get_yticks()
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
plt.savefig('Figures\DrumDlt.png')

# Plot 2: Inertia and mass of drum
fig, ax1 = plt.subplots()
plt.title('Drum Mass Properties')
ax1.plot(D, mass, 'b', label='mass m')
ax1.set_ylabel('Mass [lbs]')
ax1.set_xlabel('Drum OD [in]')
ax2 = ax1.twinx()
ax2.plot(D, I_m, 'g', label='Inertia I')
ax2.set_ylabel('Inertia [kg*m^2]')
# Get legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)
plt.savefig('Figures\DrumMassProps.png')

# Plot 3
fig2, ax1 = plt.subplots()
plt.title('Gearbox Output vs Drum Diameter D')
ax1.plot(D, T_max,'b', label='Max. T')
ax1.plot(D, T_cont, 'b--', label='Cont. T')
ax1.set_ylabel('Torque T [kN*m]')
ax1.set_xlabel('Drum OD [in]')
ax2 = ax1.twinx()
ax2.plot(D, n_max, 'g', label='Max. n')
ax2.plot(D, n_cont, 'g--', label='Cont. n')
ax2.set_ylabel('Speed n [RPM]')
# Get legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)
plt.savefig('Figures\GearboxOutput.png')

for i in range(0, len_vect, round(D_max/10)):
    print('For a %.0f" OD by %.2f" long drum, a thickness of %.3f" is required' % (D[i], l[i], t[i]))

CalcExp = np.asarray([D, l, t, q])
np.savetxt('Data\DrumParams.csv', np.transpose(CalcExp), delimiter=",")
MassExp = np.asarray([D, l, t, A, mass, I_m])
np.savetxt('Data\Mass.csv', np.transpose(MassExp), delimiter=",")






