#The initial condition
import numpy as np

left = 0
right = 2 * np.pi
s_time = 0
e_time = 20

def F(phi):
    temp = phi*phi -1
    return 25 * temp*temp

def f(phi):
    return 100 * phi * (phi*phi - 1)

def L1(v,dx):#the minus laplacian operator for periodic boundary condition
    B = 2 * v
    B[1:] -= v[:-1]
    B[:-1] -= v[1:]
    B[0] -= v[-1]
    B[-1] -= v[0]
    B /= dx * dx
    return B

def v_prod(u,v,dx):
    return dx*np.vdot(u,v)

def inverse(b,dt,N,dx): #inverse by rfft
    m = np.arange(N//2+1)

    l1 = 2 / dx * np.sin(m * np.pi / N)
    l1 = - l1 * l1

    c = 1/dt
    l = l1 - c #eigenvalue of (-L-1/dt)

    b_hat = np.fft.rfft(b)
    u_hat = b_hat/l
    return np.fft.irfft(u_hat)

def Energy(u,dx):
    return v_prod(u,L1(u,dx),dx)/2+dx*np.sum(F(u))

def Energy1(u,x,dx):
    return dx*np.sum(F(u))

def initial(X):
    return np.sin(X) 