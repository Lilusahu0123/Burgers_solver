import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *
parser = argparse.ArgumentParser()
parser.add_argument('-nx', type=int, help='Number of gridpoints', default=201)
parser.add_argument('-cfl', type=float, help='for stability', default=0.5)
parser.add_argument('-tf', type=float,help='final time',default=1)
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',default=1)
parser.add_argument('-ic', choices=('hat'), help='Init cond', default='hat')
parser.add_argument('-scheme', choices=('Gudonov', 'Lax_friedrich', 'Murman_roe', 'Eo', 'Lax_wendroff'), help='which scheme', default='Lax_friedrich')
args = parser.parse_args()
scheme = args.scheme
nx = args.nx

# Parameters
  # Number of grid points
x = np.linspace(-2, 2, nx)  # Spatial grid
dx = x[1] - x[0]  # Grid spacing
cfl = args.cfl 
t=0 # CFL number
tf = args.tf  # Final time
args = parser.parse_args()
plot_freq = args.plot_freq
initial_condition = args.ic 
if initial_condition == 'hat':  
    G = hat(x)
    Ge = hat(x)




if plot_freq >0:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1,line2 = ax.plot(x, G, 'ro',x, Ge, 'b')
    #line1, = ax.plot(x, u, 'o')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.title('ng='+str(nx)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
    plt.legend(('Numerical','Exact'))
    plt.xlim(-2,2)
    plt.ylim(-0.5,2)
    plt.grid(True); plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

# Flux function


# Initialize solution
if initial_condition =='hat':
    u=hat(x) # Initial condition: u=1 for x <= 0, u=0 for x > 0


def Gudonov():
    nx, cfl, tf, t =args.nx, args.cfl, args.tf, 0, 
    while t <= tf:
        un = np.copy(u)  # Store the old solution
        m = np.abs(np.max(un))  # Maximum wave speed
        dt = cfl * dx / (m )  # Time step with safety to avoid division by zero

    # Update solution using the numerical flux
        for i in range(1, nx - 1):  # Exclude boundary points
            fluxR = max(f(max(un[i], 0)), f(min(un[i + 1], 0)))  # Flux at i+1/2
            fluxL = max(f(max(un[i - 1], 0)), f(min(un[i], 0)))  # Flux at i-1/2
            u[i] = un[i] - dt / dx * (fluxR - fluxL)  # Update rule
            

    # Apply boundary conditions
    #u[0] = 1 # Left boundary
    #u[nx - 1] = 0  # Right boundary

        t += dt  # Advance time
        if plot_freq > 0:
            if initial_condition == 'hat':
                ue = exact_solution(x,t)
            else:
                ue=0
            
            line1.set_ydata(u)
            line2.set_ydata(ue)
            plt.xlim(-2,2)
            plt.ylim(-0.5,2)
            plt.title('ng='+str(nx)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
            plt.draw(); plt.pause(0.1)
    plt.show()

def Murman_roe():
    nx, cfl, tf, t =args.nx, args.cfl, args.tf, 0, 
    while t <= tf:
        un = np.copy(u)  # Store the old solution
        m = np.abs(np.max(un))  # Maximum wave speed
        dt = cfl * dx / (m )  # Time step with safety to avoid division by zero

    # Update solution using the numerical flux
        for i in range(1, nx - 1):  # Exclude boundary points
            fluxR = (f(un[i])+f(un[i+1]))/2 - np.abs(Murman_Roe(un[i],un[i+1]))*(un[i+1]-un[i])  # Flux at i+1/2
            fluxL = (f(un[i-1])+f(un[i]))/2 - np.abs(Murman_Roe(un[i-1],un[i]))*(un[i]-un[i-1])  # Flux at i-1/2
            u[i] = un[i] - dt / dx * (fluxR - fluxL)  # Update rule
            

    # Apply boundary conditions
    #u[0] = 1 # Left boundary
    #u[nx - 1] = 0  # Right boundary

        t += dt  # Advance time
        if plot_freq > 0:
            if initial_condition == 'hat':
                ue = exact_solution(x,t)
            else:
                ue=0
            
            line1.set_ydata(u)
            line2.set_ydata(ue)
            plt.xlim(-2,2)
            plt.ylim(-0.5,2)
            plt.title('ng='+str(nx)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
            plt.draw(); plt.pause(0.1)
    plt.show()

def Eo():
    nx, cfl, tf, t =args.nx, args.cfl, args.tf, 0, 
    while t <= tf:
        un = np.copy(u)  # Store the old solution
        m = np.abs(np.max(un))  # Maximum wave speed
        dt = cfl * dx / (m )  # Time step with safety to avoid division by zero

    # Update solution using the numerical flux
        for i in range(1, nx - 1):  # Exclude boundary points
            fluxR =f(max(un[i],0))+f(min(un[i+1],0))-f(0)  # Flux at i+1/2
            fluxL = f(max(un[i-1],0))+f(min(un[i],0))-f(0)  # Flux at i-1/2
            u[i] = un[i] - dt / dx * (fluxR - fluxL)  # Update rule
            

    # Apply boundary conditions
    #u[0] = 1 # Left boundary
    #u[nx - 1] = 0  # Right boundary

        t += dt  # Advance time
        if plot_freq > 0:
            if initial_condition == 'hat':
                ue = exact_solution(x,t)
            else:
                ue=0
            
            line1.set_ydata(u)
            line2.set_ydata(ue)
            plt.xlim(-2,2)
            plt.ylim(-0.5,2)
            plt.title('ng='+str(nx)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
            plt.draw(); plt.pause(0.1)
    plt.show()


def Lax_wendroff():
    nx, cfl, tf, t =args.nx, args.cfl, args.tf, 0, 
    while t <= tf:
        un = np.copy(u)  # Store the old solution
        m = np.abs(np.max(un))  # Maximum wave speed
        dt = cfl * dx / (m )  # Time step with safety to avoid division by zero

    # Update solution using the numerical flux
        for i in range(1, nx - 1):  # Exclude boundary points
            fluxR =(f(un[i])+f(un[i+1]))/2 - (dt/dx)*(un[i]+un[i+1])/2*(f(un[i+1])-f(un[i]))  # Flux at i+1/2
            fluxL = (f(un[i-1])+f(un[i]))/2 - (dt/dx)*(un[i-1]+un[i])/2*(f(un[i])-f(un[i-1]))   # Flux at i-1/2
            u[i] = un[i] - dt / dx * (fluxR - fluxL)  # Update rule
            

    # Apply boundary conditions
    #u[0] = 1 # Left boundary
    #u[nx - 1] = 0  # Right boundary

        t += dt  # Advance time
        if plot_freq > 0:
            if initial_condition == 'hat':
                ue = exact_solution(x,t)
            else:
                ue=0
            
            line1.set_ydata(u)
            line2.set_ydata(ue)
            plt.xlim(-2,2)
            plt.ylim(-0.5,1.5)
            plt.title('ng='+str(nx)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
            plt.draw(); plt.pause(0.1)
    plt.plot(x,u,'-r')
    plt.show()

def Lax_friedrich():
    nx, cfl, tf, t =args.nx, args.cfl, args.tf, 0,
    while t <= tf:
        un = np.copy(u)  # Store the old solution
        m = np.abs(np.max(un))  # Maximum wave speed
        dt = cfl * dx / (m )  # Time step with safety to avoid division by zero

    # Update solution using the numerical flux
        for i in range(1, nx - 1):  # Exclude boundary points
            fluxR = (f(un[i])+f(un[i+1]))/2 - (dx/(2*dt))*(un[i+1]-un[i])  # Flux at i+1/2
            fluxL = (f(un[i-1])+f(un[i]))/2 - (dx/(2*dt))*(un[i]-un[i-1])  # Flux at i-1/2
            u[i] = un[i] - dt / dx * (fluxR - fluxL)  # Update rule

    # Apply boundary conditions
    #u[0] = 1 # Left boundary
    #u[nx - 1] = 0  # Right boundary

        t += dt  # Advance time
        if plot_freq > 0:
            if initial_condition == 'hat':
                ue = exact_solution(x,t)
            else:
                ue=0
            line1.set_ydata(u)
            line2.set_ydata(ue)
            plt.xlim(-2,2)
            plt.ylim(-0.5,2)
            plt.title('ng='+str(nx)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
            plt.draw(); plt.pause(0.1)
    plt.show()


#Run the solver
if scheme == 'Gudonov':
    Gudonov()
elif scheme == 'Lax_friedrich':
    Lax_friedrich()
elif scheme == 'Murman_roe':
    Murman_roe()
elif scheme == 'Eo':
    Eo()
elif scheme == 'Lax_wendroff':
    Lax_wendroff()


