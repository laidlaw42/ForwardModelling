# TASK 3.1
# Calculate temperature field at t_end = 5000 and 50,000 years.
# Describe the changes.


import time, sys
import numpy as np
import matplotlib.pyplot as plt

# get_ipython().run_line_magic('matplotlib', 'notebook')

from IPython.display import set_matplotlib_formats

set_matplotlib_formats('png')

from ipywidgets import FloatProgress
from IPython.display import display

def task3part1(numberofyears):
    # ## Physical parameters
    Lx = 10000  # length of domain along x (m)
    Ly = 10000  # length of domain along y
    T_0 = 500  # host rock temperature
    T_1 = 900  # initial temperature of pluton 1
    T_2 = 900  # initial temperature of pluton 2
    k = 1e-6  # thermal diffusivity of rocks

    R_1 = 1200  # radius of pluton 1
    R_2 = 1200  # radius of pluton 2

    X_1 = 4000  # x coordinate of pluton 1
    Y_1 = 4000  # y coordinate of pluton 1

    X_2 = 6000  # x coordinate of pluton 2
    Y_2 = 6000  # x coordinate of pluton 2

    t_2 = 1000 * 365.25 * 24 * 3600  # time of intrusion 2 (after intrusion 1 at t=0)
    t_end = numberofyears * 365.25 * 24 * 3600  # time end of simulation

    # ## Numerical parameters
    # These are parameters that do not exist in nature, but are required for the discretisation of the problem and solving using numerical methods.
    nx = 100  # number of points along x
    ny = 100  # number of points along y

    x = np.linspace(0, Lx, nx + 1)
    y = np.linspace(0, Ly, ny + 1)

    dx = x[1] - x[0]  # resolution along x
    dy = y[1] - y[0]  # resolution along y

    mindx = min(dx, dy)  # minimum resolution
    dt = int(0.49 * (mindx * mindx / (2 * k)))  # time step

    nt_2 = int(t_2 / dt)  # number of steps to get to t_2
    nt_end = int((t_end - t_2) / dt)  # number of steps after intrusion 2 to get to the end of simulation

    # ## Initial conditions
    # Set the temperature everywhere to T_0, then set the temperature to T_1 if the distance to the center of the pluton is less than the radius of the pluton 1
    X, Y = np.meshgrid(x, y)  # make grids of x and y from vectors
    T = np.zeros(X.shape) + T_0  # make T=T_0 everywhere

    DD = np.hypot((X - X_1), (Y - Y_1))  # distance from center of pluton 1
    T[DD < R_1] = T_1  # set T to T_1 inside the pluton 1

    # ### Main loops for the first part of evolution
    # From intrusion of pluton 1 to intrusion of pluton 2.
    t = np.zeros(nt_2)  # initialise a time vector
    time = 0
    t[0] = time

    T0 = T  # save initial temperature distribution
    f = FloatProgress(min=0, max=100)  # progress bar
    display(f)

    for n in range(1, nt_2):
        T_new = np.zeros(X.shape)
        f.value = int(n / nt_2 * 100)

        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                T_new[i, j] = T[i, j] + (k * dt / (dx * dx)) * (T[i + 1, j] - 2 * T[i, j] + T[i - 1, j]) + (
                            k * dt / (dy * dy)) * (T[i, j + 1] - 2 * T[i, j] + T[i, j - 1])

        T_new[0, :] = T0[0, :]  # apply boundary conditions along four edges
        T_new[nx - 1] = T0[nx - 1]
        T_new[:, 0] = T0[:, 0]
        T_new[:, ny - 1] = T0[:, ny - 1]
        T = T_new  # update temperature field
        time = time + dt  # update time
        t[n] = time

    # ### Temperature field when second pluton intrudes
    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)
    plt.title(label='Temperature field when second pluton intrudes, t_end = ' + str(numberofyears))
    im = ax.pcolor(x, y, T, cmap='coolwarm')
    ax.set_aspect('equal')
    cax = fig.add_axes([0.85, 0.13, 0.02, 0.75])
    fig.colorbar(im, cax=cax)
    plt.show()

    # ### Next evolution
    # From intrusion of second pluton to the end of the simulation
    DD = np.hypot((X - X_2), (Y - Y_2))
    T[DD < R_2] = T_2

    t = np.zeros(nt_end)  # initialise a time vector
    t[0] = time

    T0 = T  # save initial temperature distribution

    f = FloatProgress(min=0, max=100)  # progress bar
    display(f)

    for n in range(1, nt_end):
        f.value = int(n / nt_end * 100)
        T_new = np.zeros(X.shape)
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                T_new[i, j] = T[i, j] + (k * dt / (dx * dx)) * (T[i + 1, j] - 2 * T[i, j] + T[i - 1, j]) + (
                        k * dt / (dy * dy)) * (T[i, j + 1] - 2 * T[i, j] + T[i, j - 1])

        T_new[0, :] = T0[0, :]  # apply boundary conditions along four edges
        T_new[nx - 1] = T0[nx - 1]
        T_new[:, 0] = T0[:, 0]
        T_new[:, ny - 1] = T0[:, ny - 1]
        T = T_new  # update temperature field
        time = time + dt  # update time
        t[n] = time

    # ### Temperature field at the end of simulation
    fig = plt.figure(2)
    ax = plt.subplot(1, 1, 1)
    im = ax.pcolor(x, y, T, cmap='coolwarm')
    ax.set_aspect('equal')
    cax = fig.add_axes([0.85, 0.13, 0.02, 0.75])
    fig.colorbar(im, cax=cax)
    plt.title(label='Temperature field at the end of simulation, t_end = ' + str(numberofyears))
    plt.show()

task3part1(50000)
task3part1(5000)