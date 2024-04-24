import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


ec2 = 1.44 # charge electron e^2
a = 0.55 #const
b = float((np.pi/math.sqrt(3))*a) #const
Z1 = 20 # atomic number for first nuclei
Z2 = 28 # atomic number for second nuclei
A1 = 48 # mass number for first nuclei
A2 = 72 # mass number for second nuclei
R1 = 1.28*(A1**(1/3)) - 0.76 + 0.8*(A1**(-1/3)) # the effective sharp radius for first nuclei
R2 = 1.28*(A2**(1/3)) - 0.76 + 0.8*(A2**(-1/3)) # the effective sharp radius for second nuclei
C1 = R1*(1 - (b/R1)**2) # the Süssmann central radius of the target and projectile 
C2 = R2*(1 - (b/R2)**2) # the Süssmann central radius of the target and projectile 
RC = (C1*C2)/(C1+C2) # the mean curvature radius
#s = r - C1 - C2
gamma0 = 0.9517
k = 1.7826
N1 = A1 - Z1 # the number of neutron of first nuclei
N2 = A2 - Z2 # the number of neutron of second nuclei
A = ( (N1+N2) - (Z1+Z2) )/(A1+A2) # A is the asymmetry parameter for the compound nucleus, which means drastic reduction in the magnitude of the potential with asymmetry of the colliding pair
gamma_const = gamma0*(1 - k*A) # The surface energy coefficient gamma  is defined as a function of the neutron/proton excess
rmassa = A1 * A2 / (A1 + A2) # reduce massa


def gamma(R):
    return (-0.67754 + 223358.47648 * np.exp(-0.90292 * R))/10

def U_n(r):
    if r < (1.2511 * b + (C1 + C2)):
        Un = 4 * math.pi * gamma_const * b * RC * ((-1 / 2) * ((r - (C1 + C2)) / b - 2.54) ** 2 - 0.0852 * ((r - (C1 + C2)) / b - 2.54) ** 3)
        dUn = 4 * math.pi * gamma_const * b * RC * (
                    -((-2.54 + (-C1 - C2 + r) / b) / b) - (0.2556 * (-2.54 + (-C1 - C2 + r) / b) ** 2) / b)
    else:
        Un = 4 * math.pi * gamma_const * b * RC *(-3.437)*math.exp(-((r - (C1 + C2)) / b)/0.75)
        dUn = 4 * math.pi * gamma_const * b * RC * (4.58267 * math.exp(-((1.33333 * (-C1 - C2 + r)) / b))) / b
    return Un, dUn

def Coul(r, Z1, Z2):
    Uc = ec2 * Z1 * Z2 / r
    dUc = -ec2 * Z1 * Z2 / (r ** 2)
    return Uc, dUc


def calculate_R(R, V, dt, gamma, rmassa, dUdR):
    return R + V * dt**2 - (gamma(R) * V / rmassa) * dt**2 - (dUdR / rmassa) * dt**2

def calculate_V(V, dt, gamma, rmassa, dUdR):
    return V - (gamma(R) * V / rmassa) * dt - (dUdR / rmassa) * dt


# Начальные условия
dt = 0.05
R0 = 20
R = R0
Ecm = 250
Un, dUn = U_n(R)
U_c, dUc = Coul(R, Z1, Z2)
U = Un + Uc
V0_qvadrat = 2 * (Ecm - U) / rmassa
V0 = math.sqrt(V0_qvadrat)
V = V0
num_steps = 2500
R_values = [R]
V_values = [V]

for _ in range(num_steps+1):
    Un, dUn = U_n(R)
    Uc, dUc = Coul(R, Z1, Z2)
    U = Un + Uc
    dUdR = dUn + dUc
    V = calculate_V(V, dt, gamma, rmassa, dUdR)
    R = calculate_R(R, V, dt, gamma, rmassa, dUdR)
    R_values.append(R)
    V = V - dt


for _ in range(num_steps):
    Un, dUn = U_n(R)
    Uc, dUc = Coul(R, Z1, Z2)
    U = Un + Uc
    dUdR = dUn + dUc
    V = calculate_V(V, dt, gamma, rmassa, dUdR)
    V_values.append(V)
    

ntimes = np.linspace(0, num_steps*dt, num_steps+1)
R_values = np.array(R_values[1:])
V_values = np.array(V_values)

plt.plot(ntimes, R_values)

plt.xlabel("t")
plt.ylabel("R(t) ")
plt.title("The plot of the equation R(t) for $E_{cm}$=250 MeV")
plt.show()









