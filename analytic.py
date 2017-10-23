import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

L = 400             # model length (m)
T_0 = 100           # temperature host rock (C)
T_i = 1200          # temperature intrusion (C)
k = 1e-6            # thermal diffusivity rocks (m2/s)
W = 20              # width of intrusion (m)
D = 50              # distance from center of dyke where we want temperature evolution
t = 365*24*3600     # time for which we want a temperature profile (s)
                    # and evolution of temperature at pt D until time t
tend = 1*365*24*3600      # time end of simulation is seconds
dt = 1*24*3600             # frequence of measure of temperature
nt = int(tend / (dt+1))          # number of points in time evolution
t = np.linspace(0,tend,nt)  # time vector
Te = np.zeros(t.shape)
Te[0] = T_0

for k in range(1,nt):
    A = 1 / (2 * np.sqrt(k * t[k]))
    B = np.multiply((W/2 - D),  A)
    C = np.multiply((W/2 + D), A)
    Te[k] = T_0 + np.multiply((T_i-T_0) / 2, (erf(B)+erf(C)))

plt.figure(2)
plt.plot(t,Te)
plt.ylabel('Temperature')
plt.xlabel('Time (s)')
plt.show()