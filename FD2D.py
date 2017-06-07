import numpy as np
import matplotlib.pyplot as plt
import math

# Geometrical parameters
L = 50001                            # length of domain (m)
k = 1e-6                             # thermal diffusivity (m/s2)

Xp = -10000                          # X coordinate of centre of intrusion (m)
Yp = 10000                           # Y coordinate of center of intrusion (m)
R = 5000                             # radius of intrusion (m)

# Numerical parameters
dx = 200                            # resolution x vector (m)
nx = int((L+1)/dx)                  # number of points
x = np.linspace(-L/2, L/2, num=nx)  # x vector

dy = 200                            # resolution y vector (m)
ny = int((L+1)/dx)                  # number of points
y = np.linspace(-L/2, L/2, num=ny)  # y vector

xv, yv = np.meshgrid(x, y)          # grid x,y

# temperature
T0 = 200                            # initial temperature everywhere (C)
Ti = 1200                           # initial temperature intrusion (C)
T = T0*np.ones((nx,ny))             # initial temperature

# time
tt = 6e12                           # time (s)
dt = 0.1*(dx*dx / (2*k))            # time step (s)
nt = int((tt) / dt)                 # number of time steps


print("number of timesteps is %d" % nt)

for i in range(0,nx):
    for j in range(0,ny):
        D = math.hypot(xv[i,j]-Xp, yv[i,j]-Yp) # distance from center
        if D < R:
            T[i,j] = Ti             # if within radius T = Ti

t = np.zeros(nt)
time = 0
t[0] = time
for n in range(1,nt):
    print(n)
    T_new = np.zeros((nx,ny))
    for i in range(1, nx - 1):
        for j in range(1,ny-1):
            T_new[i,j] = T[i,j] + (k*dt/(dx*dx))*(T[i+1,j] - 2*T[i,j] + T[i-1,j]) + (k*dt/(dy*dy))*(T[i,j+1] - 2*T[i,j] + T[i,j-1])

    T_new[0,:] = T0                # apply boundary conditions along four edges
    T_new[nx-1] = T0
    T_new[:,0] = T0
    T_new[:,ny-1] = T0
    T = T_new                      # update temperature field
    time = time + dt               # update time
    t[n] = time


plt.contourf(xv, yv, T)
plt.axis('equal')
plt.colorbar()
plt.show()