# TASK 1.2
# Run the simulation in order to get temperature evolution over a year at 10, 50, 100 and 150 m.
# Save the figures and place them in the report.
# Using the mouse get the time coordinate of the peak in each figure.
# Plot the position of peak versus time. Describe and comment the plot.


import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')
import numpy as np
from scipy.special import erf

def temperatureEvolution(D):
    #Physical parameters
    L = 400             # model length (m)
    T_0 = 100           # temperature host rock (C)
    T_i = 1200          # temperature intrusion (C)
    k = 1e-6            # thermal diffusivity rocks (m2/s)
    W = 20              # width of intrusion (m)


    month = (365 * 24 * 60 * 60) / 12       # defines the length of a month in seconds
    howManyMonths = 12
    t = month * howManyMonths     # time for which we want a temperature profile (s) and evolution of temperature at pt D until time t


    #Numerical parameters
    dx = 1      # point every 1
    nx = L+1    # number of points
    x = np.linspace(-L/2, L/2, num=nx)  # defines the position in the profile

    # Temperature profile at time t
    A = 1 / (2 * np.sqrt(k * t))
    B = np.multiply((W/2 - x),  A)
    C = np.multiply((W/2 + x), A)
    T = T_0 + np.multiply((T_i-T_0) / 2, (erf(B)+erf(C)))


    #Temperature evolution at point D
    tend = 1*365*24*3600      # time end of simulation in seconds
    dt = 24*3600             # time interval between measure
    nt = int(tend / (dt+1))          # number of points in time evolution
    t = np.linspace(0,tend,nt)  # time vector
    Te = np.zeros(t.shape)
    Te[0] = T_0

    for k in range(1,nt):
        A = 1 / (2 * np.sqrt(k * t[k]))
        B = np.multiply((W / 2 - D), A)
        C = np.multiply((W / 2 + D), A)
        Te[k] = T_0 + np.multiply((T_i-T_0) / 2, (erf(B)+erf(C)))

    # plots the maximum y value
    plt.ylabel('Temperature')
    plt.xlabel('Time (s)')
    fig, ax = plt.subplots()


    # def annot_max(T, ax=None):
    #     ymax = Te.max()
    #     xmax = t[np.argmax(T)]
    #     text= " y = {:.3f}".format(ymax)
    #     if not ax:
    #         ax = plt.gca()
    #     bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    #     arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    #     kw = dict(xycoords='data',textcoords="axes fraction",
    #               arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    #     ax.annotate(text, xy=(xmax,ymax))
    #
    # annot_max(Te)

    plt.figure(2)
    plt.title(label=('Temperature evolution at dyke width of ' + str(D) + ' metres'))
    plt.plot(t,Te)
    plt.ylabel('Temperature')
    plt.xlabel('Time (s)')
    plt.show()

    ymaxTe = Te.max()

    print('-----------------------------------------------------------------------------------------------------------')
    print('')
    print('Temperature evolution over a year with a dyke width of ' + str(D) + ' m. Peak temperature value is ' + str(ymaxTe) + ' C')
    print('')


temperatureEvolution(10)
temperatureEvolution(50)
temperatureEvolution(100)
temperatureEvolution(150)
