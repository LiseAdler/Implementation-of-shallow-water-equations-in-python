# -*- coding: utf-8 -*-
#%% libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

save_path = 'C:/Users/lisea/OneDrive - TU Wien/Bachelorarbeit/Figures/'

#%% 1. Required weir width 

def width(Q, C, h):
    l = Q/(C*(h**(3/2)))
    return l

Q = 0.25 #discharge[m³/s]
C = 1.4 #discharge coefficient
h = 0.1 #head above spillway
l = width(Q,C,h)

#%% 3. Differential equations for water level progression

def dhdt(h, t, a, b, d, Cr, Sr, C, l, S):
    ''' a,b = precipitation coefficients, d = duration of precipitation, Cr = runoff coefficient, Sr = runoff area,
    C = discharge coefficient, l = width of weir, h = head above spillway
    h[m], t[h]
    '''
    #Rain
    i = 0.001 * (a*(d**(-b)))
    #Runoff
    Qr = Cr * i * Sr
    if t <= d:
        R_t = Qr*(t/d)
        ip = i
    if d < t <= 2*d:
        R_t = Qr*(2-(t/d))
        ip = 0
    if t > 2*d:
        R_t = 0
        ip = 0
    # discharge
    Q = C * l* (h**(3/2))
    #differential equation
    dhdt = ip + (R_t/S)- (3600* Q/S)
    return dhdt

def solve_dhdt(h0, t, a, b, d, Cr, Sr, C, l, S):
    sol = odeint(dhdt, h0, t, args=(a, b, d, Cr, Sr, C, l, S))
    return sol

#%% 4. Plots water level progression example and water level progression with duration and width variation
T = 5000
a = 8.3 + 2.45*np.log(T)
b = 0.47
d = 5
Cr = 1
Sr = 20000
C = 1.4
l = 1
S = 5000

h0 = 0 
t = np.linspace(0, 24, 100)

#Variations
d_values = [0.1, 1, 10]
l_values = [0.1, 1, 5, 10]



sol = solve_dhdt(h0, t, a, b, d, Cr, Sr, C, l, S)
i = 0.001 * (a*(d**(-b)))
Qr = Cr * i * Sr
t_1 = np.linspace(0,d,100)
t_2 = np.linspace(d, 2*d, 100)
runoff_1 = 1000 * (Qr/S*(t_1/d)/S) 
runoff_2 = 1000 * (Qr/S*(2-(t_2/d))/S)
discharge = 3600* C * (l/S)* (sol[:,0]**(3/2))

#plots
fig, ax = plt.subplots( figsize=(12, 8))
ax.plot(t, sol, color = 'blue', label = 'Water height [m]')
ax.plot(t_1, runoff_1, color = 'pink', linestyle ='--')
ax.plot(t_2, runoff_2, color = 'pink', label = 'Runoff [m/h]', linestyle ='--')
ax.plot(t, discharge, color = 'cyan',label = 'Discharge[m/h]', linestyle ='--')
ax.hlines(i, 0, d, color = 'darkgreen', label = 'Precipitation[m]', linestyle ='--')
ax.legend()
#ax.set_title('duration = 5h & width = 1m')
ax.set_xlabel('t [h]')
plt.savefig(save_path + 'weir_water_progression.png')



fig, ax = plt.subplots(figsize= (12,8))
for d in d_values:
    for j, l in enumerate(l_values):
        sol = solve_dhdt(h0, t, a, b, d, Cr, Sr, C, l, S)
        ax.plot(t, sol, label=f'd = {d}h , l = {l}m')
ax.set_title('Variations in duration and width')
ax.legend()
ax.set_xlabel('t [h]')
ax.set_ylabel('water head [m]')
plt.savefig(save_path + 'weir_var_duration_width.png')

#%% 5. Minimum weir width for flood mitigation
d_values = np.arange(0.1, 3.1, 0.1)
l_values = np.arange(1,6)

maxh = np.zeros((len(d_values), len(l_values)))

for i, d in enumerate(d_values):
    for j, l in enumerate(l_values):
        p = (d, l)
        h0 = 0
        t = np.linspace(0, 24, 500)
        h = solve_dhdt(h0, t, a, b, d, Cr, Sr, C, l, S)
        maxh[i, j] = np.max(h)

# Plots
plt.figure(figsize=(10, 6))
for j in range(len(l_values)):
    plt.plot(d_values, maxh[:, j], label=f'Weir width = {l_values[j]}m')
plt.xlabel('Duration [h]')
plt.ylabel('Maximum water level [m]')
plt.legend()
plt.grid(True)
plt.savefig(save_path + 'weir_max_head.png')

plt.show()

def min_widht(h_max, d_values, l_values, a, b, Cr, Sr, C, S):
    m =[]
    for j, l in enumerate(l_values[::-1]):
        m.append(l)
        for i, d in enumerate(d_values):
            h0 = 0
            t = np.linspace(0, 24, 500)
            h = solve_dhdt(h0, t, a, b, d, Cr, Sr, C, l, S)
            if np.max(h) > h_max:
                return m[-2]

l_min = min_widht(0.1, d_values, l_values, a, b, Cr, Sr, C, S)
