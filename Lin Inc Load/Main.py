# Sebastien Blanchet
# Altaeros Energies, Systems Inern
# Solution as per thesis paper for linearly incraesing load on drum

# Import plotting modules
import matplotlib.pyplot as plt
import numpy as np

#  Initial Terms
t = 70                  # Thick [mm]
l = 1775              # Len [mm]
R = 500                 # Rad [mm]
E_y = 2*10**5           # Yield [MPa]
v = 0.35                # Poissons ratio [ul]
zeta = 0.8              # Non dim'l param [ul]
T_rope = 1.96*10**3     # Tension [kN]
d_rope = 65              # Rope diam [mm]

# Hoop pressure, [MPa]
P = 2*T_rope/d_rope

# k term
# w_p = ((P*(R**2))/(E_y*t))*(1-zeta)
# k = ((3*(1-v**2))**0.25)/((R*t)**0.5)

w_p=0.125
k=0.0071

# F terms
F_1 = (np.sinh(k*l)**2)-(np.sin(k*l)**2)
F_2 = (np.sinh(k * l) ** 2) + (np.sin(k * l) ** 2)
F_3 = (np.sinh(k * l)) *(np.cosh(k*l)) + (np.sin(k * l)) * (np.cos(k*l))
F_4 = (np.sinh(k * l)) * (np.cosh(k * l)) - (np.sin(k * l)) * (np.cos(k * l))
F_5 = (np.sin(k*l))**2
F_6 = (np.sinh(k * l)) ** 2
F_7 = np.cosh(k*l*zeta)*np.cos(k*l*zeta)
F_8 = np.sinh(k * l * zeta) * np.sin(k * l * zeta)
F_9 = (np.cosh(k*l*zeta)*np.sin(k*l*zeta))-(np.sinh(k*l*zeta)*np.cos(k*l*zeta))
F_15 = (np.cosh(k * l * zeta) * np.sin(k * l * zeta))
F_16 = (np.sinh(k * l * zeta) * np.cos(k * l * zeta))

# S terms
S_1 = -(P*R**2)/(E_y*t)
S_2 = (-(P*R**2)/(E_y*t))*((F_3/F_1)-((1/k*l)*((F_6/F_1)+(F_8/F_1))))
S_3 = ((P * R ** 2) / (E_y*t)) * ((F_3 / F_1) - ((1 / k * l) * ((F_5 / F_1) + (F_8 / F_1))))
S_4 = ((P * R ** 2) / (E_y*t)) * ((F_2 / F_1) - ((1 / k * l) * ((F_4 / F_1) + (F_9 / F_1))))

# Begin solving for t
w_tot = w_p+(S_1*F_7)+(S_2*F_15)+(S_3*F_16)+(S_4*F_8)
Sig = (E_y/R)*w_tot

