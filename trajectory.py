#privet

#import math
# import numpy as np
# import matplotlib.pyplot as plt
#
# ec2 = 1.44
# a = 0.55
# b = float((np.pi/math.sqrt(3))*a)
# Z1 = 40
# Z2 = 60
# A1 = 90
# A2 = 144
# R1 = 1.28*(A1**(1/3)) - 0.76 + 0.8*(A1**(-1/3))
# R2 = 1.28*(A2**(1/3)) - 0.76 + 0.8*(A2**(-1/3))
# C1 = R1*(1 - (b/R1)**2)
# C2 = R2*(1 - (b/R2)**2)
# RC = (C1*C2)/(C1+C2)
# #s = r - C1 - C2
# gamma0 = 0.9517
# k = 1.7826
# N1 = A1 - Z1
# N2 = A2 - Z2
# A = ( (N1+N2) - (Z1+Z2) )/(A1+A2)
# gamma_const = gamma0*(1 - k*A)
# #s = r - C1 - C2
# #e = s/b
# e = 1.21
# e2 = 1.25
# #r1 = b*e1 + C1 + C2
# #r2 = b*e2 + C1 + C2
# rmassa = A1 * A2 / (A1 + A2)
#
# def gamma(R0):
#     return -0.67754 + 223358.47648 * math.exp(-0.90292 * R0)
#
#
#
# def calculate_next_R(R0, V, dt, gamma, rmassa, dUdR):
#     return R0 - V * dt**2 - (gamma(R0) * V / rmassa) * dt**2 - (dUdR * rmassa) * dt**2
#
#
#
# def U_n(r):
#     if r < (1.2511 * b + (C1 + C2)):
#         un = 4.0 * math.pi * gamma_const * b * RC * (
#                     -((-2.54 + (-C1 - C2 + r) / b) / b) - (0.2556 * (-2.54 + (-C1 - C2 + r) / b) ** 2) / b)
#         dUn = 4 * math.pi * gamma_const * b * RC * (
#                     -((-2.54 + (-C1 - C2 + r) / b) / b) - (0.2556 * (-2.54 + (-C1 - C2 + r) / b) ** 2) / b)
#     else:
#         un = 4.0 * math.pi * gamma_const * b * RC * (4.58267 * math.exp(-((1.33333 * (-C1 - C2 + r)) / b))) / b
#         dUn = 4 * math.pi * gamma_const * b * RC * (4.58267 * math.exp(-((1.33333 * (-C1 - C2 + r)) / b))) / b
#     return un, dUn
#
# def Coul(r, Z1, Z2):
#     U_c = ec2 * Z1 * Z2 / r
#     dU_c = -ec2 * Z1 * Z2 / (r ** 2)
#     return U_c, dU_c
#
# def calculate_trajectory():
#     V0 = 0.0
#     Ecm = 250.0
#     dt = 0.05
#     rmax = 20
#     R0 = rmax
#     dmr = 0.0
#     Ekin = Ecm
#     times = []
#     positions = []
#
#     Rt = R0
#     V = V0
#
#     for itr in range(201):
#         t = itr * dt
#         times.append(t)
#         positions.append(Rt)
#         un , dUn = U_n(Rt)
#         U_c, dU_c = Coul(Rt, Z1, Z2)
#         U = un + U_c
#         dUdR = dUn + dU_c
#
#         Rt = calculate_next_R(R0, V, dt, gamma, rmassa, dUdR)
#
#     return times, positions
#
# times, positions = calculate_trajectory()
#
# plt.plot(times, positions)
# plt.xlabel("t")
# plt.ylabel("R(t)")
# plt.title("Траектория")
# plt.show()



import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

ec2 = 1.44
a = 0.55
b = float((np.pi/math.sqrt(3))*a)
Z1 = 40
Z2 = 60
A1 = 90
A2 = 144
R1 = 1.28*(A1**(1/3)) - 0.76 + 0.8*(A1**(-1/3))
R2 = 1.28*(A2**(1/3)) - 0.76 + 0.8*(A2**(-1/3))
C1 = R1*(1 - (b/R1)**2)
C2 = R2*(1 - (b/R2)**2)
RC = (C1*C2)/(C1+C2)
#s = r - C1 - C2
gamma0 = 0.9517
k = 1.7826
N1 = A1 - Z1
N2 = A2 - Z2
A = ( (N1+N2) - (Z1+Z2) )/(A1+A2)
gamma_const = gamma0*(1 - k*A)
#s = r - C1 - C2
#e = s/b
e = 1.21
e2 = 1.25
#r1 = b*e1 + C1 + C2
#r2 = b*e2 + C1 + C2
rmassa = A1 * A2 / (A1 + A2)
def gamma(R):
    return -0.67754 + 223358.47648 * np.exp(-0.90292 * R)

def U_n(r):
    if r < (1.2511 * b + (C1 + C2)):
        un = 4.0 * math.pi * gamma_const * b * RC * (
                    -((-2.54 + (-C1 - C2 + r) / b) / b) - (0.2556 * (-2.54 + (-C1 - C2 + r) / b) ** 2) / b)
        dUn = 4 * math.pi * gamma_const * b * RC * (
                    -((-2.54 + (-C1 - C2 + r) / b) / b) - (0.2556 * (-2.54 + (-C1 - C2 + r) / b) ** 2) / b)
    else:
        un = 4.0 * math.pi * gamma_const * b * RC * (4.58267 * math.exp(-((1.33333 * (-C1 - C2 + r)) / b))) / b
        dUn = 4 * math.pi * gamma_const * b * RC * (4.58267 * math.exp(-((1.33333 * (-C1 - C2 + r)) / b))) / b
    return un, dUn

def Coul(r, Z1, Z2):
    U_c = ec2 * Z1 * Z2 / r
    dU_c = -ec2 * Z1 * Z2 / (r ** 2)
    return U_c, dU_c

def calculate_next_R(R, V, dt, gamma, rmassa, dUdR):
    return R + V * dt**2 + (gamma(R) * V / rmassa) * dt**2 - (dUdR * rmassa) * dt**2

import math
import numpy as np
import matplotlib.pyplot as plt

import scipy as sp

ec2 = 1.44
a = 0.55
b = float((np.pi/math.sqrt(3))*a)
Z1 = 20
Z2 = 28
A1 = 48
A2 = 72
R1 = 1.28*(A1**(1/3)) - 0.76 + 0.8*(A1**(-1/3))
R2 = 1.28*(A2**(1/3)) - 0.76 + 0.8*(A2**(-1/3))
C1 = R1*(1 - (b/R1)**2)
C2 = R2*(1 - (b/R2)**2)
RC = (C1*C2)/(C1+C2)
#s = r - C1 - C2
gamma0 = 0.9517
k = 1.7826
N1 = A1 - Z1
N2 = A2 - Z2
A = ( (N1+N2) - (Z1+Z2) )/(A1+A2)
gamma_const = gamma0*(1 - k*A)
#s = r - C1 - C2
#e = s/b
e = 1.21
e2 = 1.25
#r1 = b*e1 + C1 + C2
#r2 = b*e2 + C1 + C2
rmassa = A1 * A2 / (A1 + A2)
def gamma(R):
    return (-0.67754 + 223358.47648 * np.exp(-0.90292 * R))/10

def U_n(r):
    if r < (1.2511 * b + (C1 + C2)):
        un = 4 * math.pi * gamma_const * b * RC * ((-1 / 2) * ((r - (C1 + C2)) / b - 2.54) ** 2 - 0.0852 * ((r - (C1 + C2)) / b - 2.54) ** 3)
        dUn = 4 * math.pi * gamma_const * b * RC * (
                    -((-2.54 + (-C1 - C2 + r) / b) / b) - (0.2556 * (-2.54 + (-C1 - C2 + r) / b) ** 2) / b)
    else:
        un = 4 * math.pi * gamma_const * b * RC *(-3.437)*math.exp(-((r - (C1 + C2)) / b)/0.75)
        dUn = 4 * math.pi * gamma_const * b * RC * (4.58267 * math.exp(-((1.33333 * (-C1 - C2 + r)) / b))) / b
    return un, dUn

def Coul(r, Z1, Z2):
    U_c = ec2 * Z1 * Z2 / r
    dU_c = -ec2 * Z1 * Z2 / (r ** 2)
    return U_c, dU_c


def calculate_next_R(R, V, dt, gamma, rmassa, dUdR):
    return R + V * dt**2 - (gamma(R) * V / rmassa) * dt**2 - (dUdR / rmassa) * dt**2

def calculate_next_V(V, dt, gamma, rmassa, dUdR):
    return V - (gamma(R) * V / rmassa) * dt - (dUdR / rmassa) * dt
# Начальные условия
R0 = 20
R = R0
Ecm = 250

un, dUn = U_n(R)
U_c, dU_c = Coul(R, Z1, Z2)
U = un + U_c
V0_qvadrat = 2 * (Ecm - U) / rmassa
V0 = math.sqrt(V0_qvadrat)
dt = 0.05


V = V0
num_steps = 2500
R_values = [R]


for _ in range(num_steps+1):
    Un, dUn = U_n(R)
    U_c, dU_c = Coul(R, Z1, Z2)
    U = Un + U_c
    dUdR = dUn + dU_c
    V = calculate_next_V(V, dt, gamma, rmassa, dUdR)
    R = calculate_next_R(R, V, dt, gamma, rmassa, dUdR)
    R_values.append(R)
    V = V - dt

V_values = [V]
for _ in range(num_steps):
    Un, dUn = U_n(R)
    U_c, dU_c = Coul(R, Z1, Z2)
    U = Un + U_c
    dUdR = dUn + dU_c
    V = calculate_next_V(V, dt, gamma, rmassa, dUdR)
    V_values.append(V)


#R = np.linspace(2, 20, 100)
#gamma_plot = [gamma(i) for i in R]


# Построение графика

ntimes = np.linspace(0, num_steps*dt, num_steps+1)
R_values = np.array(R_values[1:])
V_values = np.array(V_values)

plt.plot(ntimes, V_values)

plt.xlabel("t")
plt.ylabel("V(t) ")
#plt.title("The plot of the equation R(t) for $E_{cm}$=250 MeV")

#plt.axis([0, 60, 5, 21])


plt.show()

# # Начальные условия
# R0 = 20
# V0 = 0
# dt = 0.05
#
# R = R0
# V = V0
# num_steps = 200
# R_values = [R]
#
#
# for _ in range(num_steps+1):
#     un, dUn = U_n(R)
#     U_c, dU_c = Coul(R, Z1, Z2)
#     U = un + U_c
#     dUdR = dUn + dU_c
#     R = calculate_next_R(R, V, dt, gamma, rmassa, dUdR)
#     R_values.append(R)
#
#
#
#
# # Построение графика
# times = np.linspace(0, num_steps*dt, num_steps+1)
# R_values = np.array(R_values[1:])
# plt.plot(times, R_values)
# plt.xlabel("t")
# plt.ylabel("R(t)")
# plt.title("График уравнения R(t), скорость при Ecm=")
# plt.show()








