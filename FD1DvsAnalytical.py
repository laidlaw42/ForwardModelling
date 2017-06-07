import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# Physical parameters

L = 1000
T0 = 300
Ti = 1200

k = 1e-6
W = 30

dx = 1
nx = L+1
x = np.linspace(-L/2, L/2, num=nx)

t = 500000000


dt = 0.2*dx*dx / (2*k)

nsteps = int(t /dt)

A = 1 / (2 * np.sqrt(k * t))
B = np.multiply((W/2 - x),  A)
C = np.multiply((W/2 + x), A)
Ta = T0 + np.multiply((Ti-T0) / 2, (erf(B)+erf(C)))

T = T0*np.ones(nx)

for index, item in enumerate(x):
    if abs(item) < W/2:
        T[index] = Ti


Ts = T

time = 0
for n in range(1,nsteps):
    T_new = np.zeros(nx)
    for i in range(1,nx-1):
        T_new[i] = T[i] + (k*dt/(dx*dx))*(T[i-1]+T[i+1]-2*T[i])

    T_new[0] = T[0]
    T_new[nx-1] = T[nx-1]
    T = T_new
    time = time + dt


diff = T-Ta
plt.figure(1)
plt.plot(x,Ts)
plt.plot(x,Ta)
plt.plot(x,T)
plt.plot(x,diff)
plt.ylabel('Temperature')
plt.xlabel('Distance')

plt.show()