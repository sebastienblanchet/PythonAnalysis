# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# March 27 2017
# Note: all calculations to be completed in US customary

# Import plotting modules
import numpy as np
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Define all parameters
d_tether = 15.2 / 25.4  # tether diameter 15.2mm to in
L_tether = (370 * 1000) / 25.4  # tether length overall 370m to in
P = 224.809 * 51.264  # force P kN
SF = 1.5  # safety factor
S_y = 50763  # Yield strength 350 MPa
E_y = 2.9 * 10 ** 7  # Youngs modulus 200GPa
v = 0.3  # poissons ratio
rho = 0.284  # density of steel lb/in^3
SG = 1.01  # spool gap 1%
D = 28
rt = 2 ** 0.5
c = 1

# Calculate length
l = np.ceil((SG * d_tether * L_tether) / (np.pi * (2 * D + 3 * d_tether)))

# Create vectors
theta = np.linspace(0, 2 * np.pi, 16)
x_l = np.arange(0, l)
THETA, X = np.meshgrid(theta, x_l)
row_num = len(x_l)
col_num = len(theta)
stack_num = 50
w = np.zeros((row_num, col_num, stack_num))
A_n = np.zeros(stack_num)
B_n = np.zeros(stack_num)
k_n = np.zeros(stack_num)
dif = np.zeros(stack_num)
itnum = np.zeros(stack_num)

# Convergence params
tol = 5
n_max = stack_num

# Values for 0th index of infinite sum
A_n[0] = -P / (rt * np.pi) * np.sin((np.pi * c) / l)
itnum[0] = 1
dif[0] = 100

# For and While loop info
# cycle through rows
# for i in range(row_num):
#     # cycle through cols
#     for j in range(col_num):
#         ROWS X COLS !!!!!!!
#         A[i][j]
# While loops only check condition after control returns back to it (i.e. for loop is done)
# while diffcheck <= tol:

# Sum from n=1 to max interation number : range(start,stop,step)
for n in range(1, n_max):
    # set iteration number (due to index start at 0)
    n_eval = n + 1

    # A_n coefficient
    A_n[n] = -P / (rt * np.pi * n_eval) * np.sin((n_eval * np.pi * c) / l)
    A_n[n] = A_n[n] + A_n[n - 1]

    # Error calculation % and summation term n
    dif[n] = abs(((A_n[n] - A_n[n - 1]) / A_n[n - 1])) * 100
    itnum[n] = n_eval

    # Check if difference is less than tol
    if dif[n] <= tol:
        # get out of for loop
        break

# Convergence plot, all values up to n_end (i.e. final iteration) + show
plt.plot(itnum[:n], dif[:n])
plt.show()
