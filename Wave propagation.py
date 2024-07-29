# -*- coding: utf-8 -*-
#%% libraries
import numpy as np
import matplotlib.pyplot as plt

save_path = 'C:/Users/lisea/OneDrive - TU Wien/Bachelorarbeit/Figures/'

#%% parameters
tc = 3600 # [s]time of concentration(peak flow)
Qmax = 1.8 # [m³/s] peak flow
Qn = 0.95 # [m³/s] normal flow
L = 2 # [m]width
K = 20 # Manning-Strickler coefficient
S0 = 0.001 # slope
#S0 = np.tan(theta) #if theta is given as a slope
nx = 21 # spatial resolution(uneven number)
tmax = 10 * tc #max time for calculation
t_steps = 21 #time resolution
nt = 6 #flood deformation timesteps
#%% initial condition floodwave
def floodwave_parameters(Qmax, Qn, L, K, S0):
    #maximum flowdepth and normal flowdepth
    hmax = (Qmax/(L*K*(S0)**(1/2)))**(3/5)
    hn = (Qn/(L*K*(S0)**(1/2)))**(3/5)
    # characteristic length of floodwave 
    lc = (tc/L)*((Qmax-Qn)/(hmax-hn))
    return [hmax, hn, lc]

hmax, hn, lc = floodwave_parameters(Qmax, Qn, L, K, S0)

# initial condition
def h_x_0(hmax, hn, lc, nx): 
    x0 = np.linspace(-2 *lc, 2*lc, nx)
    h = np.zeros_like(x0) 
    k = (hmax-hn)/lc #slope flood wave
    
    for i in range(len(x0)):
        if x0[i] <= -lc:
            h[i] = hn
        elif x0[i] <= 0:
            h[i] = x0[i] * k + hmax
        elif x0[i] <= lc:
            h[i] = -x0[i] * k + hmax
        else:
            h[i] = hn
    return x0, h

x0, h0 = h_x_0(hmax, hn, lc, nx)

fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(x0, h0, label='Water surface profile of the floodwave')
ax.axhline(hn, color='red', linestyle='--', label='Normal flow depth')
ax.axhline(hmax, color='green', linestyle='--', label='Peak flow depth')
ax.fill_between(x0, h0, hn, where=(h0 >= hn), color='blue', alpha=0.2)
ax.set_xlabel('x [m]')
ax.set_ylabel('h [m]')
ax.legend(loc ='upper left')
#ax.set_title('initial condition')
plt.savefig(save_path + 'wave_int_cond.png')

    
#%% Characteristics

def characteristics(K, S0, h0, x0, tmax, t_steps):
    t = np.linspace(0, tmax, t_steps)
    X = np.zeros((len(h0), t_steps))
    for i in range(len(h0)):
        c = (5/3)*K *np.sqrt(S0)*(h0[i]**(2/3))
        X[i,:] = x0[i] + c*t
    return X


X = characteristics(K, S0, h0, x0,tmax ,t_steps)

#plot characteristics
t = np.linspace(0, tmax, t_steps) #time values

fig, ax = plt.subplots(figsize=(12, 6))
for Xi in X:
    ax.plot(Xi, t, color='orange', alpha=0.5)
ax.set_xlabel('x [m]')
ax.set_ylabel('t [s]')
ax.set_title('characteristics')
plt.savefig(save_path + 'characteristics.png')


# plot characteristics and initial condition
fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(x0, h0, color='orange', alpha=0.5, label = 'Characteristics')
ax.plot(x0, h0, label='Initial flood wave')

ax1 = ax.twinx()
for Xi in X:
    ax1.plot(Xi, t, color='orange', alpha=0.5)
        
ax.set_xlabel('x [m]')
ax.set_ylabel('Water level of the initial condition [m]')
ax1.set_ylabel('Time [s]')
#ax.set_title('characteristics and initial condition')
ax.legend(loc ='lower right')
plt.savefig(save_path + 'characteristics_initial_condition.png')


#%% Deformation and propagation model
def deformation(x0, h0, nt, tmax, S0):
    X1 = x0.copy()
    X2 = np.zeros_like(x0)
    results = np.zeros((nt +1, len(x0)))
    results[0] = X1
    
    for time in range(nt):
        X2 = np.zeros_like(x0)
        for i in range(nx):
            X2[i] = X1[i] + (tmax/nt) * K * np.sqrt(S0) * (5/3) * (h0[i]**(2/3))
        
        results[time+1] = X2
        X1 = X2.copy()
    return results
        
results = deformation(x0, h0, nt, tmax, S0)

# 2D plot
fig, ax = plt.subplots(figsize=(12, 6))
for Xi in results:
    ax.plot(Xi, h0, color='darkblue', alpha=0.5)
ax.set_xlabel('x [m]')
ax.set_ylabel('Water level of the flood wave [m]')
#ax.set_title('Deformation of floodwave overtime')
plt.savefig(save_path + '2d_wave_prop.png')


# 2D plot with characteristics
fig, ax = plt.subplots(figsize=(12, 6))
for Xi in results[0:nt+1]:
    ax.plot(Xi, h0, color='cyan', alpha=0.5)
ax1 = ax.twinx()
for Xi in X:
    ax1.plot(Xi, t, color='orange', alpha=0.5)
        
ax.set_xlabel('x [m]')
ax.set_ylabel('h [m]')
ax1.set_ylabel('t [s]')
ax.set_title('floodwave deformation with characteristics')

# 3D plot

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

time_values = np.linspace(0, tmax, nt +1)
t = np.linspace(0, tmax, t_steps) #time values

for i in range(7):
    ax.plot(results[i], h0, zs=time_values[i], zdir='y', label=f'time={round(time_values[i]/3600, 1)} hours')

ax.set_xlabel('X [m]')
ax.set_ylabel('Time [s]')
ax.set_zlabel('H [m]')
ax.set_title('Deformation of floodwave over time')
ax.legend()
plt.show()

# 3D plot with characteristics
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')


for Xi in X:
    ax.plot(Xi, t, zs = hn,  color='orange', alpha=0.5)
    
for i in range(nt +1):
    ax.plot(results[i], h0, zs=time_values[i], zdir='y', label=f'Time = {round(time_values[i]/3600,1)} hours')
    
ax.set_xlabel('x [m]')
ax.set_ylabel('Time [s]')
ax.set_zlabel('H [m]')
#ax.set_title('Deformation of floodwave over time')
ax.legend()
plt.savefig(save_path + '3d_wave_prop.png')

plt.show()

#%% Shock initiation time and position

def shock_time(K, S0, hn, hmax, lc):
    k = (hmax - hn)/lc
    ts = 9/10 *((hn**(1/3))/(K*np.sqrt(S0)*k))
    return ts
def shock_position(K, S0, hn, ts):
    xs = 5/3 * K *np.sqrt(S0)* (hn**(2/3)* ts)
    return xs

ts = shock_time(K, S0, hn, hmax, lc)
xs = shock_position(K, S0, hn, ts)

# plots
# only characteristics
fig, ax = plt.subplots(figsize=(12, 6))
for Xi in X[0:10]:
    ax.plot(Xi, t, color='orange', alpha=0.5)
ax.scatter(xs, ts, label = 'shock')
ax.set_xlabel('x [m]')
ax.set_ylabel('t [s]')
ax.set_title('characteristics')

# 2D
fig, ax = plt.subplots(figsize=(12, 6))
for Xi in results[0:nt+1]:
    ax.plot(Xi, h0, color='cyan', alpha=0.5)
ax1 = ax.twinx()
for Xi in X:
    ax1.plot(Xi, t, color='orange', alpha=0.5)
ax1.scatter(xs, ts, label = 'shock')            
ax.set_xlabel('x [m]')
ax.set_ylabel('h [m]')
ax1.set_ylabel('t [s]')
ax.set_title('characteristics and initial condition')

# 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
time_values = np.linspace(0, tmax, nt+1)
for Xi in X:
    ax.plot(Xi, t, zs = hn,  color='orange', alpha=0.5)    
for i in range(nt+1):
    ax.plot(results[i], h0, zs=time_values[i], zdir='y', label=f'Time = {round(time_values[i]/3600,1)} hours')
ax.scatter(xs, ts, hn, label = 'Initial shock')    
ax.set_xlabel('X [m]')
ax.set_ylabel('Time [s]')
ax.set_zlabel('H [m]')
#ax.set_title('Deformation of floodwave over time')
ax.legend()
plt.savefig(save_path + 'shock.png')

plt.show()
