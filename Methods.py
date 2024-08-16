from init import *
# Original solution
def Solution(phi_init,Nx,times,s_time, e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]

    #Solve by RFFT
    for t in T:
        b = f(phi0) - phi0/dt
        phi1 = inverse(b,dt,Nx,dx)
        
        if np.max(abs((phi1-phi0)/dt+L1(phi1,dx,)+f(phi0)))>1e-5:
            print(f"error! at t={t} in phi")
        
        E.append(Energy(phi1,dx))
        phi0 = phi1
    return (E,T)

def SPC(phi_init,Nx,times,s_time, e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    for t in T:
        b = 2*f(phi0)-2*phi0/dt+L2(phi0,x,dx)-L1(phi0,dx)
        phi_h = inverse(b,dt/2,Nx,dx)
        b = 2*f(phi_h)+L1(phi0,dx)-2*phi0/dt+2*(L2(phi_h,dx)-L1(phi_h,x,dx))
        phi1 = inverse(b,dt/2,Nx,dx)

        E.append(Energy(phi1,dx))
        phi0=phi1
    return (E,T)

def SAV(phi_init,Nx,times,s_time, e_time, left, right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E1 = np.empty(times+1)
    r = np.empty(times+1)
    now = 0
    E1[now] = dx*np.sum(F(phi0))
    E = [v_prod(phi0,L1(phi0,dx),dx)/2 + E1[now]]
    E_bar = E.copy()
    r[now] = np.sqrt(E1[now]+C0)

    #Solve
    for t in T:
        B = f(phi0)/np.sqrt(E1[now]+C0)
        C = B * r[now] - B/2 * v_prod(B,phi0,dx) - phi0/dt

        b = inverse(B,dt,Nx,dx)
        gamma = inverse(C,dt,Nx,dx)

        phi1 = gamma + v_prod(B,gamma,dx) / 2 / (1 - 0.5 * v_prod(B,b,dx)) * b
        
        now+=1
        r[now] = v_prod(B,phi1-phi0,dx)/2 + r[now-1]
        E0 = v_prod(phi1,L1(phi1,dx),dx)/2
        E1[now] = dx*np.sum(F(phi1))
        E.append(E0 + E1[now])
        E_bar.append(E0 + r[now] * r[now] - C0)
        phi0 = phi1
        
    return (E_bar,E,E1,r,T)

def SAV2(phi_init,Nx,times,s_time, e_time, left, right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E1 = np.empty(times+1)
    r = np.empty(times+1)
    now = 0
    E1[now] = Energy1(phi0,x,dx)
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    r[now] = np.sqrt(E1[now]+C0)

    #Solve
    for t in T:
        phi_bar = inverse(f(phi0)-2*phi0/dt,dt/2,Nx,dx)
        temp = Energy1(phi_bar,x,dx)
        B = f(phi_bar)/np.sqrt(temp+C0)
        v = -B/2
        C = 2*B*r[now] +v*v_prod(B,phi0,dx) - 2*phi0/dt +L1(phi0,dx)

        b = inverse(B,dt/2,Nx,dx)
        gamma = inverse(C,dt/2,Nx,dx)

        phi1 = gamma - v_prod(v,gamma,dx)/(1+v_prod(v,b,dx))*b
        
        now+=1
        r[now] = v_prod(B,phi1-phi0,dx)/2 + r[now-1]
        E0 = v_prod(phi1,L1(phi1,dx),dx)/2
        E1[now] = Energy1(phi1)
        E.append(E0 + E1[now])
        E_bar.append(E0 + r[now] * r[now] - C0)
        phi0 = phi1
    
    return (E_bar,E,E1,r,T)

def N_SAV(phi_init,Nx,times,s_time, e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi = phi_init
    E = [Energy(phi,dx)]
    E_bar = E.copy()
    R = E[0]
    theta = 1

    #Solve
    for t in T:
        b = f(phi)-phi/dt
        phi_temp = inverse(b,dt,Nx,dx)
        E_temp = Energy(phi_temp,dx)
        mu = f(phi_temp)+L1(phi_temp,dx)
        xi = R/(E_temp+dt*v_prod(mu,mu,dx))
        phi = (theta+(1-theta)*xi)*phi_temp

        E.append(Energy(phi,dx))
        R = xi * E_temp
        E_bar.append(R)
        theta = xi

    return (E_bar,E,T)

def N_SAV1(phi_init,Nx,times,s_time, e_time,left,right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    q0 = E[0] + C0
    xi = 1

    #Solve
    for t in T:
        phi1 = inverse(xi*f(phi0)-phi0/dt,dt,Nx,dx)

        """
        if np.max(abs((phi1-phi0)/dt+L1(phi1,dx,dy)+xi*f(phi0)))>1e-5:
            print(f"error! at t={t}, error in phi")
        """

        mu = L1(phi1,dx) +xi * f(phi0)
        q1 = q0 - xi * dt * v_prod(mu,mu,dx)
        E.append(Energy(phi1,dx))
        E_bar.append(q1-C0)
        xi = q1/(E[-1]+C0)
        q0 = q1
        phi0 = phi1

    return (E_bar,E,T)

def N_SAV1b(phi_init,Nx,times,s_time, e_time,left,right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    q0 = E[0] + C0
    xi = 1

    #Solve
    for t in T:
        phi1 = inverse(xi*f(phi0)-phi0/dt,dt,Nx,dx)

        mu = L1(phi1,dx) +xi * f(phi0)
        E.append(Energy(phi1,dx))
        q1 = q0/(1 + dt /E[-1] * v_prod(mu,mu,dx))
        xi = q1/(E[-1]+C0)
        E_bar.append(q1-C0)
        q0 = q1
        phi0 = phi1

    return (E_bar,E,T)

def N_SAV1c(phi_init,Nx,times,s_time, e_time,left,right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    q0 = E[0] + C0
    xi = 1

    #Solve
    for t in T:
        phi1 = inverse(xi*f(phi0)-phi0/dt,dt,Nx,dx)

        mu = L1(phi1,dx) + f(phi0)
        E.append(Energy(phi1,dx))
        q1 = q0/(1 + dt /E[-1] * v_prod(mu,mu,dx))
        xi = q1/(E[-1]+C0)
        E_bar.append(q1-C0)
        q0 = q1
        phi0 = phi1

    return (E_bar,E,T)

def N_SAV1d(phi_init,Nx,times,s_time, e_time,left,right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    q0 = E[0] + C0
    xi = 1

    #Solve
    for t in T:
        phi1 = inverse(xi*f(phi0)-phi0/dt,dt,Nx,dx)

        mu = L1(phi1,dx) + f(phi1)
        E.append(Energy(phi1,dx))
        q1 = q0/(1 + dt /E[-1] * v_prod(mu,mu,dx))
        xi = q1/(E[-1]+C0)
        E_bar.append(q1-C0)
        q0 = q1
        phi0 = phi1

    return (E_bar,E,T)

def N2_SAV1b(phi_init,Nx,times,s_time, e_time,left,right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    q0 = E[0] + C0

    #Solve
    for t in T:
        B = q0*f(phi0)/(E[-1]+C0)
        phi_bar = inverse(B-2*phi0/dt,dt/2,Nx,dx)
        mu = B+L1(phi_bar,dx)
        E_temp = Energy(phi_bar,dx)+C0
        q_bar = q0/(1+dt/2/(E_temp)*v_prod(mu,mu,dx))
        
        B = q_bar*f(phi_bar)/(E_temp)
        phi1 = inverse(L1(phi0,dx)+2*B-2*phi0/dt,dt/2,Nx,dx)
        mu = L1(phi1,dx)+B
        E.append(Energy(phi1,dx))
        E_temp = E[-1]+C0
        b = dt/2/E_temp*v_prod(mu,mu,dx)
        q1 = q0*(1-b)/(1 + b)
        E_bar.append(q1-C0)
        q0 = q1
        phi0 = phi1
    return (E_bar,E,T)

def N_SAV2(phi_init,Nx,times,s_time, e_time,left,right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)
    
    # initial value
    phi0 = phi_init
    E = [v_prod(phi0,L1(phi0,dx),dx)/2 + dx * np.sum(F(phi0))]
    E_bar = E.copy()
    q0 = E[0]+C0
    xi = 1

    #Solve
    for t in T:

        phi1 = inverse(-phi0/dt+xi*f(phi0),dt,Nx,dx)

        E.append(v_prod(phi1,L1(phi1,dx),dx)/2 + dx * np.sum(F(phi1)))
        q1 = E[-1]+C0
        xi = q1/q0
        E_bar.append(q1-C0)
        q0 = q1
        phi0 = phi1
        
    return (E_bar,E,T)

def gPAV1(phi_init,Nx,times,s_time, e_time,left,right, C0):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    q0 = E[0] + C0

    #Solve
    for t in T:

        u1 = inverse(-phi0/dt,dt,Nx,dx)
        u2 = inverse(f(phi0),dt,Nx,dx)
        phi_temp = u1 + u2
        L_temp = L1(phi_temp,dx)
        mu = f(phi_temp)+L_temp
        E_temp = v_prod(phi_temp,L_temp,dx)/2 + dx * np.sum(F(phi_temp))+C0
        q1 = q0/(1 + dt * v_prod(mu,mu,dx)/E_temp)
        phi1 = u1 + q1/E_temp * u2

        E.append(Energy(phi1,dx))
        E_bar.append(q1-C0)
        q0 = q1
        phi0 = phi1

    return (E_bar,E,T)

#gPAV second order
def gPAV2(phi_init,Nx,times,s_time, e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    q0 = E[0]

    u1 = inverse(-phi0/dt,dt,Nx,dx)
    u2 = inverse(-f(phi0),dt,Nx,dx)

    phi_a = u1 + u2
    E_a = Energy(phi_a,dx)
    mu_a = L1(phi_a,dx) + f(phi_a)
    V = dt*v_prod(mu_a,mu_a,dx)
    xi_a = q0/(E_a+V)
    q_a = xi_a * E_a
    Tilde_phi3_2 = (3*phi_a-phi0)/2
    q1_2 = (q_a+q0)/2
    E_bar = [q1_2]
    Tilde_E3_2 = Energy(Tilde_phi3_2,dx)
    xi = q1_2/(Tilde_E3_2+V)
    q3_2 = xi * Tilde_E3_2
    #q1 = (2 * q3_2+q0)/3
    phi1 = u1 + xi * u2
    E.append(Energy(phi1,dx))
    E_bar.append(q3_2)

    #Solve
    for t in T[1:]:
        gamma = 3/2
        phi_hat = 2 * phi1 - phi0 / 2
        phi_bar = 2 * phi1 - phi0
        u1 = inverse(-phi_hat/dt,dt/gamma,Nx,dx)
        u2 = inverse(f(phi_bar),dt/gamma,Nx,dx)
        phi2_temp = u1 + u2
        phi3_2_temp = (3*phi2_temp-phi1)/2
        E3_2 = Energy(phi3_2_temp,dx)
        mu2_temp = L1(phi2_temp,dx) + f(phi2_temp)
        xi = q3_2/(E3_2+dt*v_prod(mu2_temp,mu2_temp,dx))
        phi2 = u1 + xi * u2
        q3_2 = xi*E3_2
        #q1 = (2*q3_2+q1)/3

        E.append(Energy(phi2,dx))
        E_bar.append(q3_2)
        phi0 = phi1
        phi1 = phi2

    return (E_bar,E,T)

def gPAV2a(phi_init,Nx,times,s_time, e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = [E[0]]
    q0 = E[0]

    #Solve
    for t in T:
        phi_bar = inverse(2 * f(phi0)+L1(phi0,dx)-2*phi0 / dt,dt/2,Nx,dx)
        phi_half = (phi_bar+phi0)/2
        mu = f(phi_half)+L1(phi_half,dx)
        c = 0.5*v_prod(mu,mu,dx)/(Energy(phi_half,dx))
        q1 = q0 * (1/dt - c) / (1/dt + c)
        q_half = (q1 +q0)/2
        phi1 = inverse(2 * q_half * f(phi_half) / Energy(phi_half,dx) +L1(phi0,dx) - 2*phi0/dt,dt/2,Nx,dx)

        E.append(Energy(phi1,dx))
        E_bar.append(q1)
        phi0 = phi1
        q0 =q1

    return (E_bar,E,T)

def gPAV2b(phi_init,Nx,times,s_time, e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times)*dt
    x, dx = np.linspace(left,right,Nx,retstep = True, endpoint = False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = [E[0]]
    q0 = E[0]

    #Solve
    for t in T:
        phi_half = inverse(f(phi0) - 2 * phi0/dt,dt/2,Nx,dx)
        mu = f(phi_half)+L1(phi_half,dx)
        c = 0.5*v_prod(mu,mu,dx)/(Energy(phi_half,dx))
        q1 = q0 * (1/dt - c) / (1/dt + c)
        q_half = (q1 +q0)/2
        phi1 = inverse(2 * q_half * f(phi_half) / Energy(phi_half,dx) +L1(phi0,dx) - 2*phi0/dt,dt/2,Nx,dx)

        E.append(Energy(phi1,dx))
        E_bar.append(q1)
        phi0 = phi1
        q0 =q1

    return (E_bar,E,T)

def ESAV1(phi_init,Nx,times,s_time,e_time,left,right,C0):
    dt = (e_time-s_time)/times
    T = np.arange(times)*dt
    x,dx = np.linspace(0,2*np.pi,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E1 = np.empty(times+1)
    s = np.empty(times+1)

    E1[0] = dx * np.sum(F(phi0))
    i = 0
    E = [v_prod(phi0,L1(phi0,dx),dx)/2 + E1[0]]
    E_bar = E.copy()
    s[0] = np.exp(E1[i]/C0)
    b = f(phi0)

    #Solve
    for t in T:

        phi1 = inverse(-phi0/dt+b,dt,Nx,dx)
        
        i += 1
        s[i] = s[i-1]*np.exp(v_prod(b,phi1-phi0,dx)/C0)
        E1[i] = dx * np.sum(F(phi1))
        E0 = v_prod(phi1,L1(phi1,dx),dx)/2
        E.append(E0 + E1[i])
        E_bar.append(E0 + C0 * np.log(s[i]))
        
        b = s[i]/np.exp(E1[i]/C0) * f(phi1)
        phi0 = phi1
    
    return (E_bar,E,E1,s,T)  

def ESAVII(phi_init,Nx,times,s_time,e_time,left,right,C0):
    dt = (e_time-s_time)/times
    T = np.arange(times)*dt
    x,dx = np.linspace(0,2*np.pi,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E1 = np.empty(times+1)
    s = np.empty(times+1)

    E1[0] = Energy1(phi0)
    i = 0
    E = [Energy(phi0)]
    E_bar = E.copy()
    s[0] = np.exp(E1[i]/C0)
    
    phi_half = inverse(f(phi0)-2*phi0/dt,dt/2,Nx,dx)
    s_half = np.exp(Energy1(phi_half)/C0)
    b = f(phi_half)
    #Solve
    for t in T:
        phi1 = inverse(-2*phi0/dt+2*b+L1(phi0,dx),dt/2,Nx,dx)

        i += 1
        s[i] = s[i-1]*np.exp(v_prod(b,phi1-phi0,dx)/C0)
        E1[i] = Energy1(phi1)
        E0 = v_prod(phi1,L1(phi1,dx),dx)/2
        E.append(E0 + E1[i])
        E_bar.append(E0 + C0 * np.log(s[i]))

        phi_half = (3*phi1-phi0)/2
        phi0 = phi1
        s_half = (3*s[i]-s[i-1])/2
        b = s_half/np.exp(Energy1(phi_half)/C0) * f(phi_half)
        
    return (E_bar,E,E1,s,T)

def ESAV2(phi_init,Nx,times,s_time,e_time,left,right,C0):
    dt = (e_time-s_time)/times
    T = np.arange(times)*dt
    x,dx = np.linspace(0,2*np.pi,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    E_bar = E.copy()
    s0 = np.exp(E[0]/C0)
    xi = 1

    m = np.arange(Nx//2+1)

    #Solve
    for t in T:

        phi1 = inverse(-phi0/dt+xi*f(phi0),dt,Nx,dx)

        """
        if np.max(abs((phi1-phi0)/dt+L1(phi1,dx,dy)+xi*f(phi0)))>1e-5:
            print(f"error! at t={t}, error in phi")
        """

        mu = L1(phi1,dx) +xi * f(phi0)
        s1 = s0*np.exp(-dt * v_prod(mu,mu,dx)/C0)
        E.append(Energy(phi1,dx))
        E_bar.append(C0 * np.log(s1))
        xi = s1/np.exp(E[-1]/C0)
        s0 = s1
        phi0 = phi1
        
    return (E_bar,E,T)

def New_ESAV(phi_init,Nx,times,s_time,e_time,left,right,C0):
    dt = (e_time-s_time)/times
    T = np.arange(times)*dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi = phi_init
    E = [Energy(phi,dx)]
    E_bar = E.copy()
    S = np.exp(E[0]/C0)

    #Solve
    for t in T:

        phi_bar = inverse(f(phi)-phi/dt,dt,Nx,dx)
        mu = f(phi_bar)+L1(phi_bar,dx)
        S = S/(1+dt/C0*v_prod(mu,mu,dx))
        xi = S/np.exp(Energy(phi_bar,dx)/C0)
        U = xi*(2-xi)
        phi = U*phi_bar
        E.append(Energy(phi,dx))
        E_bar.append(np.log(S)*C0)
        
    return (E_bar,E,T)

def A1(phi_init,Nx,times,s_time, e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]

    #Solve by RFFT
    for t in T:
        b = f(phi0) - 2 * phi0/dt
        phi_bar = inverse(b,dt/2,Nx,dx)
        #phi_half = (phi_bar+phi0)/2
        #phi1 = inverse(2 * f(phi_half) +L1(phi0,dx) - 2*phi0/dt,dt/2,Nx,dx)
        phi1 = inverse(2 * f(phi_bar) +L1(phi0,dx) - 2*phi0/dt,dt/2,Nx,dx)
        
        E.append(Energy(phi1,dx))
        phi0 = phi1
    return (E,T)

def A2(phi_init,Nx,times,s_time, e_time,left,right, L):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]

    #Solve by RFFT
    for t in T:
        b = f(phi0) - 2 * phi0/dt - L1(phi0,dx)
        Q = 2/dt + M
        phi_bar = inverse(b,dt,N,dx)
        phi1 = inverse(2 * f(phi_bar) +L1(phi0,dx) - 2*phi0/dt,dt/2,Nx,dx)
        
        E.append(Energy(phi1,dx))
        phi0 = phi1
    return (E,T)

def SA(phi_init,Nx,times,s_time, e_time,left,right, L):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi = phi_init
    E = [Energy(phi,dx)]

    #Solve by RFFT
    for t in T:
        phi_bar = phi - dt * (L1(phi,dx)+f(phi))/2
        b = 2 * f(phi_bar) - 2 * phi/dt - 2 * L * phi_bar + L * phi + L1(phi,dx)
        phi = inverse(b,dt/(L*dt + 2),Nx,dx)
        
        E.append(Energy(phi,dx))
    return (E,T)

def SB(phi_init,Nx,times,s_time, e_time,left,right, L):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi = phi_init
    E = [Energy(phi,dx)]

    #Solve by RFFT
    for t in T:
        b = f(phi) - 2 * phi/dt
        phi_bar = inverse(b,dt/2,Nx,dx)
        b = 2 * f(phi_bar) - 2 * phi/dt - 2 * L * phi_bar + L * phi + L1(phi,dx)
        phi = inverse(b,dt/(L*dt + 2),Nx,dx)
        
        E.append(Energy(phi,dx))
    return (E,T)

def SC(phi_init,Nx,times,s_time, e_time,left,right,L):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]
    LL = L*L

    for t in T:
        b = 2*f(phi0)-2*phi0/dt-LL*phi0
        phi_h = inverse(b,dt/(2+LL*dt),Nx,dx)
        b = 2*f(phi_h)+L1(phi0,dx)-2*phi0/dt+LL*phi0-2*LL*phi_h
        phi1 = inverse(b,dt/(2+LL*dt),Nx,dx)

        E.append(Energy(phi1,dx))
        phi0=phi1
    return (E,T)

def SD(phi_init,Nx,times,s_time, e_time,left,right, L, M):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi = phi_init
    E = [Energy(phi,dx)]

    #Solve by RFFT
    for t in T:
        phi_bar = phi - dt * (L1(phi,dx)+f(phi)-L*phi)/2
        phi_bar /= (1+L*dt/2)
        b = (2*(f(phi_bar)-M*L1(phi_bar,dx)-L*phi_bar-phi/dt)+L*phi)/(1+M)+L1(phi,dx)
        phi = inverse(b,(1+M)*dt/(L*dt + 2),Nx,dx)
        
        E.append(Energy(phi,dx))
    return (E,T)

def SAttempt(phi_init,Nx,times,s_time,e_time,left,right):
    dt = (e_time-s_time)/times
    T = s_time + np.arange(times) * dt
    x,dx = np.linspace(left,right,Nx,retstep=True, endpoint=False)

    # initial value
    phi0 = phi_init
    E = [Energy(phi0,dx)]

    C = (np.exp(dt/(dt+0.01))-1)/(np.exp(1)-1)

    for t in T:
        #D = v_prod(f(phi0),f(phi0),dx)
        #D = D/(3*D+0.1)
        Amp = 1 - C
        b = Amp * f(phi0) - phi0/dt
        phi1 = inverse(b,dt,Nx,dx)

        E.append(Energy(phi1,dx))
        phi0 = phi1
    return (E,T)