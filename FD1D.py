import numpy as np
import matplotlib.pyplot as plt

L = 2000                # length of the spatial domain
X1 = -50                # position dyke 1
X2 = 50               # position dyke 2

T1 = 1200               # temperature dyke 1
T2 = 1200               # temperature dyke 2
T0 = 200                # temperature host rock

W1 = 20                 # width dyke 1
W2 = 20                 # width dyke 2

t1 = 0                  # time intrusion dyke 1
t2 = 6e8                # time intrusion dyke 2

k = 1e-6

dx = 1
nx = L+1
x = np.linspace(-L/2, L/2, num=nx)                            # resolution of the space vector
print(x[int(nx/2)])
dt = 0.2*(dx*dx / (2*k))

T = T0*np.ones(nx)                         # define host temperature everywhere

dX1 = x - X1                               # distances from dyke 1 center
for index, item in enumerate(dX1):         # if within dyke 1 then initial temp is T1
     if abs(item) < W1/2:
         T[index] = T1

Ts = T
plt.figure(1)
plt.plot(x,Ts)
plt.ylabel('Temperature')
plt.xlabel('Distance')


nt1 = int((t2-t1) / dt)
print(nt1)

time = 0
for n in range(1,nt1):
    T_new = np.zeros(nx)
    for i in range(1,nx-1):
        T_new[i] = T[i] + k*dt*(T[i-1]-2*T[i]+T[i+1])/dx**2
    T_new[0] = T[0]
    T_new[nx-1] = T[nx-1]
    T = T_new
    time = time + dt


dX2 = x - X2                               # distances from dyke 1 center
for index, item in enumerate(dX2):         # if within dyke 1 then initial temp is T1
     if abs(item) < W2/2:
         T[index] = T2
Ts2 = T
plt.plot(x,Ts2)

nt2 = nt1
print(nt2)

for n in range(1,nt2):
    T_new = np.zeros(nx)
    for i in range(1,nx-1):
        T_new[i] = T[i] + k*dt*(T[i-1]-2*T[i]+T[i+1])/dx**2
    T_new[0] = T[0]
    T_new[nx-1] = T[nx-1]
    T = T_new
    time = time + dt

plt.plot(x,T)


plt.show()