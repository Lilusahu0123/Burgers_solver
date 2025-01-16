import numpy as np
from sympy import symbols, diff


def hat(x):
    u=np.zeros_like(x)
    u[abs(x)<1]=1
    u[abs(x)>=1]=0
    return u



def exact_solution(x, t):
    """Exact solution for u_t + (u^2 / 2)_x = 0 with 'hat' initial condition."""
    u_exact = np.zeros_like(x)
    for i, xi in enumerate(x):
        # Solve using the method of characteristics
        if xi  <= -1:  # Within the rarefaction wave region
            u_exact[i] = 0 
        elif -1 <= xi < t-1:  # Constant left state
            u_exact[i] = (xi+1)/t
        elif t-1<xi<t/2+1:  # Constant right state
            u_exact[i] = 1
        elif xi>=t/2+1:
            u_exact[i]=0
    return u_exact

def Murman_Roe(u,v):
    if u !=v :
        linear_function = (f(u)-f(v))/(u-v)
    elif u==v:
        linear_function = u
        
    return linear_function


def f(x):
    return x**2 / 2


    