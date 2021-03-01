# # Heat equation - Cooling of intrusions
# ## Part 2: Finite difference - 1D
# The analytical solution for the cooling of a dyke only works if we make the asumptions that the cooling is one
# dimensional and there is only one heat source. Here we employ a numerical technique for a similar setup to demonstrate
# that numerical techniques allow investigating more complex situations. We will calculate the temperature profile
# across two parallel dykes. The situation is still one-dimensional but the temperature evolution is more complex.
# The heat equation and the employed finite difference scheme (Forward Euler) are detailed in the hand-out.


# TASK 2.1
# Run the simulation to calculate the temperature profile at years 3, 5, 10 and 20.
# Save the plots and place them in the report.
# Describe and comment on the temperature profiles.

import numpy as np
import matplotlib.pyplot as plt

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')

from ipywidgets import FloatProgress
from IPython.display import display


def part2task1(numberOfYears):
    # ### Physical parameters
    # Here we indicate the physical parameters of the simulation
    X1 = -50       # position dyke 1
    X2 = 50        # position dyke 2

    T1 = 1200      # initial temperature dyke 1
    T2 = 1200      # initial temperature dyke 2
    T0 = 200       # temperature host rock

    W1 = 20        # width dyke 1
    W2 = 20        # width dyke 2

    t1 = 0
    t2 = 2*365*24*3600      # time intrusion dyke 2

    tend = numberOfYears*365*24*3600  # time end of simulation

    k = 1e-6      # thermal diffusivity of rocks


    # ### Numerical parameters
    # Here we define how the physical parameters are going to be discretized.

    L = 2000                               # length of profile
    dx = 1                                 # resolution in profile
    nx = L + 1                             # number of points in the profile
    x = np.linspace(-L/2,L/2,nx)           # vector coordinates of points in the profile
    dt = 0.4*(dx*dx/(2*k))                 # time step
    Tc = np.zeros((int(tend/dt),1))        # vector temperature in central point
    tc = np.linspace(0,tend,int(tend/dt))  # vector time


    # ### Initial conditions
    # Defines the initial temperature profile
    T = T0 * np.ones(nx)                   # defines host temperature everywhere
    dX1 = x - X1                           # distance from center of dyke 1
    for index,item in enumerate(dX1):
        if abs(item) < W1/2:
            T[index] = T1                  # if within dyke 1 then temperature is T1


    # ### First time loop
    # We calculate the evolution of the temperature along the profile from t1 to t2, using the Forward Euler Scheme

    nt1 = int((t2-t1)/dt)                  # number of time steps to produce to get to t2
    ctime = 0

    f = FloatProgress(min=0,max=100)       # create a progress bar to see how the simulation goes
    display(f)

    for n in range(0,nt1):                 # main loop for time
        T_new = np.zeros(nx)               # new temperature profile
        f.value = int(n/nt1*100)

        for i in range(1,nx-1):            # main loop for space (along the profile)
            T_new[i] = T[i] + (k*dt/(dx*dx))*(T[i-1]-2*T[i]+T[i+1])

        T_new[0] = T[0]                    # apply the boundary condition (T remains constant at edges of profile)
        T_new[nx-1] = T[nx-1]
        T = T_new
        Tc[n] = T[int(nx/2)]               # get temperature at center point
        ctime = ctime+dt                   # update time


    # ### Change initial condition at time t2
    # We take the output of the first part of the simulation and modify it to ccount for the second intrusion

    dX2 = x-X2
    for index, item in enumerate(dX2):
        if abs(item) < W2/2:
            T[index] = T2

    nt2 = int((tend-t2)/dt)                # number of time steps to get from t2 to tend
    f = FloatProgress(min=0,max=100)       # create a progress bar to see how the simulation goes
    display(f)

    for n in range(0,nt2):                 # main loop for time
        T_new = np.zeros(nx)               # new temperature profile
        f.value = int(n/nt1*100)

        for i in range(1,nx-1):            # main loop for space (along the profile)
            T_new[i] = T[i] + (k*dt/(dx*dx))*(T[i-1]-2*T[i]+T[i+1])

        T_new[0] = T[0]                    # apply the boundary condition (T remains constant at edges of profile)
        T_new[nx-1] = T[nx-1]
        T = T_new
        Tc[nt1+n-1] = T[int(nx/2)]               # get temperature at center point
        ctime = ctime+dt

    def annot_max(T, ax=None):
        ymax = T.max()
        xmax = x[np.argmax(T)]
        text = " y = {:.3f}".format(ymax)
        if not ax:
            ax = plt.gca()
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        arrowprops = dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=60")
        kw = dict(xycoords='data', textcoords="axes fraction",
                  arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
        ax.annotate(text, xy=(xmax, ymax))

    annot_max(T)


    # Now we plot the temperature profile at tend

    plt.figure(1)
    plt.title(label='Temperature profile over '+ str(numberOfYears) + ' years')
    plt.plot(x,T)
    plt.xlabel("Position")
    plt.ylabel("Temperature")
    plt.show()


    Tc[0] = T0
    plt.figure(2)
    plt.plot(tc[:-2],Tc[:-2])
    plt.title(label='Temperature evolution over ' + str(numberOfYears) + ' years')
    plt.xlabel("Time")
    plt.ylabel("Temperature")
    plt.show()

    ymax = T.max()

    print('-----------------------------------------------------------------------------------------------------------')
    print('')
    print('Temperature profile over' + str(numberOfYears) + 'maximum temperature is '+ str(ymax))
    print('')
    print('-----------------------------------------------------------------------------------------------------------')

part2task1(3)
part2task1(5)
part2task1(10)
part2task1(15)