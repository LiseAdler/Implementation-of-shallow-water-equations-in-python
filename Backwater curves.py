# -*- coding: utf-8 -*-
#%% libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

save_path = 'C:/Users/lisea/OneDrive - TU Wien/Bachelorarbeit/Figures/'

#%% Flow parameters

def normal_flow_depth(theta, K, q):
    #theta= slope angle, K= Manning-Strickler coefficient, q = Q/B specific discharge
    #Jf= Bottom slope, hn = Normal flow depth
    Jf = np.sin(theta)
    hn = (q/(K*Jf**(1/2)))**(3/5)
    return hn

def critical_flow_depth(theta, q):
    # theta= slope angle, q=specific discharge
    #hc= ciritical flow depth
    g = 9.81
    hc = (q**2/g*np.cos(theta))**(1/3)
    return hc

def froudenumber(h, q, theta):
     g = 9.81
     Fr = q/((g*np.cos(theta)*h**3)**(1/2))
     return Fr
 
def flow_regime(h, q, theta):
    g = 9.81
    Fr = q/((g*np.cos(theta)*h**3)**(1/2))
    if Fr>1:
        return 'supercritical/torrential'
    if Fr<1:
        return 'subcritical/fluvial'
    if Fr==1:
        return 'critical flow'
    
    
def velocity(q, h):
    u = q/h
    return u

#%% Bresse equation
def bresse(h,x,  theta, K, q):
    g = 9.81
    roh = 1000
    u=q/h
    tau_b= (roh * g/K**2)*u**2/h**(1/3)
    dhdx= (g*np.sin(theta)-tau_b/(roh*h))/(g*np.cos(theta)-q**2/h**3)
    return dhdx

#%% solving
def solve_bresse(h0, x, theta, K, q):
    sol = odeint(bresse, h0, x, args=(theta, K, q))
    return sol
#%% Example: Dam gate

#First case
x = np.linspace(0, 10, 500)
theta = 0.8 * np.pi / 180
K = 50
q = 0.1
h0 = 0.01

bed = x * -np.tan(theta)
sol = solve_bresse(h0, x, theta, K, q)
hc = critical_flow_depth(theta, q)
hn = normal_flow_depth(theta, K, q)
sol = bed + sol[:, 0]

plt.plot(x, bed, color='brown', linestyle='-', label='Bottom slope')
plt.plot(x, sol, label='Water surface profile')
plt.plot(x,hn + bed, color='red', linestyle='--', label='Normal flow depth')
plt.plot(x, hc + bed, color='green', linestyle='--', label='Critical flow depth')
plt.fill_between(x, bed, sol, color='blue', alpha=0.2)
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.legend(loc='lower left')
#plt.title('backwater curve small dam gate 1')
plt.savefig(save_path + 'backwater_curve_small_damgate_torr.png')
plt.show()

#Second case
# with hydraulic jump at around 3.4m, model no longer valid
x = np.linspace(0, 3.4, 500)
theta = 0.25 * np.pi / 180
K = 50
q = 0.1
h0 = 0.05
h0_1 = 0.1

bed = x * -np.tan(theta)
sol = solve_bresse(h0, x, theta, K, q)
hc = critical_flow_depth(theta, q)
hn = normal_flow_depth(theta, K, q)
sol = bed + sol[:, 0]


plt.plot(x, bed, color='brown', linestyle='-', label='Bottom slope')
plt.plot(x, sol, label='Water surface profile')
plt.plot(x,hn + bed, color='red', linestyle='--', label='Normal flow depth')
plt.plot(x, hc + bed, color='green', linestyle='--', label='Critical flow depth')
plt.fill_between(x, bed, sol, color='blue', alpha=0.2)
plt.xlabel('x [m]')
plt.ylabel('z [m]')
#plt.title('backwater curve small dam gate 2')
plt.legend(loc='lower left')
plt.savefig(save_path + 'backwater_curve_small_damgate_flu.png')

plt.show()
#%% Example: Dam
# Maximum height for flood mitigation
def h_max_solution(hn, initial_h0, x, theta, K, q):
    h0 = initial_h0
    max_iterations = 1000  
    iteration = 0
    
    while iteration < max_iterations:
        sol = solve_bresse(h0, x, theta, K, q)
        h_max = sol[-1]
        deviation = np.abs(h_max - hn)
        if deviation > 0.01:
            return h0
        h0 += 0.1
        iteration += 1
    print("no solution found")
    return None

#%% Example: First case

# Plot velocity and backwater curve
h0= 1
q = 0.1
K = 50
theta = 0.3* np.pi / 180

hn = normal_flow_depth(theta, K, q)
hc = critical_flow_depth(theta, q)
upstream_flow_1 = flow_regime(hn, q, theta)

#backwatercurve boundary conditions downstream
x = np.linspace( 0, -200, 500)

bed = x * -np.tan(theta)
sol = solve_bresse(h0, x, theta, K, q)
hc = critical_flow_depth(theta, q)
hn = normal_flow_depth(theta, K, q)
Sol = bed + sol[:, 0]
u = velocity(q, sol)

cmap = plt.get_cmap('jet')

plt.plot(x, bed, color='brown', linestyle='-', label='Bottom slope')
plt.plot(x, Sol, label='Water surface profile')
plt.plot(x,hn + bed, color='red', linestyle='--', label='Normal flow depth')
plt.plot(x, hc + bed, color='green', linestyle='--', label='Critical flow depth')
for i in range(len(x) - 1):
    plt.fill_between(x[i:i+2], bed[i:i+2], Sol[i:i+2], color=cmap((u[i] - np.min(u)) / (np.max(u) - np.min(u))), alpha=0.2)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=np.min(u), vmax=np.max(u)))
sm.set_array([])
cbar = plt.colorbar(sm, alpha=0.5)
cbar.set_label('Velocity [m/s]')
#plt.title('Dam backwater curve and velocity')
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.legend()
plt.savefig(save_path + 'backwater_curve_dam_flu.png')
plt.show()

# Maximum height of the dam for flood mitigation of a village 2,5km upstream
x = np.linspace(0, -2500, 500)

h_max = h_max_solution(hn, 3, x, theta, K, q)

h0= h_max
x = np.linspace(0, -2500, 500)
sol = solve_bresse(h0, x, theta, K, q)
Sol = bed + sol[:, 0]
u = velocity(q, sol)


plt.plot(x, bed, color='brown', linestyle='-', label='Slope')
plt.plot(x, Sol, label='bresse')
plt.plot(x,hn + bed, color='red', linestyle='--', label='hn')
plt.plot(x, hc + bed, color='green', linestyle='--', label='hc')
plt.fill_between(x, bed, Sol, color='blue', alpha=0.2)
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.legend()
plt.show()

#%% Example: Second case
h0= 1
q = 0.1
K = 50
theta = 0.8* np.pi / 180
hn = normal_flow_depth(theta, K, q)
hc = critical_flow_depth(theta, q)
upstream_flow_2 = flow_regime(hn, q, theta)

#backwatercurve boundary conditions downstream
x = np.linspace( 0, -62, 500)

bed = x * -np.tan(theta)
sol = solve_bresse(h0, x, theta, K, q)
hc = critical_flow_depth(theta, q)
hn = normal_flow_depth(theta, K, q)
Sol = bed + sol[:, 0]
u = velocity(q, sol)

cmap = plt.get_cmap('jet')
plt.plot(x, bed, color='brown', linestyle='-', label='Bottom slope')
plt.plot(x, Sol, label='Water surface profile')
plt.plot(x,hn + bed, color='red', linestyle='--', label='Normal flow depth')
plt.plot(x, hc + bed, color='green', linestyle='--', label='Critical flow depth')
for i in range(len(x) - 1):
    plt.fill_between(x[i:i+2], bed[i:i+2], Sol[i:i+2], color=cmap((u[i] - np.min(u)) / (np.max(u) - np.min(u))), alpha=0.2)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=np.min(u), vmax=np.max(u)))
sm.set_array([])
cbar = plt.colorbar(sm, alpha=0.5)
cbar.set_label('Velocity [m/s]')
#plt.title('Dam backwater curve and velocity')
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.legend()
plt.savefig(save_path + 'backwater_curve_dam_torr.png')
plt.show()

#%% HYDRAULIC JUMP
def conjugate_depth(h1, q, theta):
    g = 9.81
    Fr = q/((g*np.cos(theta)*h1**3)**(1/2))
    h2 = h1*(0.5*((1+8*Fr**2)**(1/2)-1))
    return h2

def length_hjump(h1, q, theta):
    g = 9.81
    Fr = q/((g*np.cos(theta)*h1**3)**(1/2))
    L = 160* np.tanh(Fr/20)-12
    return L

def length_hjump_1(h2, theta):
    Jf = np.tan(theta)
    L = 6.1 + 4* Jf
    return L
#%% Example: seperate backwater curve upstream and downstream
# backwatercurve upstream
x = np.linspace(0, 20, 500)
theta = 0.8 * np.pi / 180
K = 50
q = 0.1
h0 = 0.05
hc = critical_flow_depth(theta, q)
hn = normal_flow_depth(theta, K, q)

bed = x * -np.tan(theta)
sol = solve_bresse(h0, x, theta, K, q)
Sol = bed + sol[:, 0]
h2_con = conjugate_depth(sol, q, theta)

fig, axs = plt.subplots(2, 2, figsize=(16, 10))
x = np.linspace(-20, 0, 500)

axs[0,0].plot(x, bed, color='brown', linestyle='-', label='Bottom slope')
axs[0,0].plot(x, Sol, label='Water surface profile')
axs[0,0].plot(x,hn + bed, color='red', linestyle='--', label='Normal flow depth')
axs[0,0].plot(x, hc + bed, color='green', linestyle='--', label='Critical flow depth')
axs[0,0].fill_between(x, bed, Sol, color='blue', alpha=0.2)
axs[0,0].set_xlabel('x [m]')
axs[0,0].set_ylabel('z [m]')
axs[0,0].set_title('Dam gate')
axs[0,0].legend(loc = 'lower left')

axs[1,0].plot(x, sol, label='Water surface profile')
axs[1,0].axhline(hn, color='red', linestyle='--', label='Normal flow depth')
axs[1,0].axhline(hc, color='green', linestyle='--', label='Critical flow depth')
axs[1,0].set_xlabel('x [m]')
axs[1,0].set_ylabel('h [m]')
axs[1,0].legend(loc = 'lower left')

x_1 = np.linspace(0, -12 ,500)
h0_1 = 0.3
bed = x_1 * -np.tan(theta)
sol_1 = solve_bresse(h0_1, x_1, theta, K, q)
Sol_1 = bed + sol_1[:, 0]
axs[0,1].plot(x_1, bed, color='brown', linestyle='-', label='Bottom slope')
axs[0,1].plot(x_1, Sol_1, label='Water surface profile')
axs[0,1].plot(x_1,hn + bed, color='red', linestyle='--', label='Normal flow depth')
axs[0,1].plot(x_1, hc + bed, color='green', linestyle='--', label='Critical flow depth')
axs[0,1].fill_between(x_1, bed, Sol_1, color='blue', alpha=0.2)
axs[0,1].set_xlabel('x [m]')
axs[0,1].set_ylabel('z [m]')
axs[0,1].legend(loc = 'lower left')
axs[0,1].set_title('Dam')

axs[1,1].plot(x_1, sol_1, label='Water surface profile')
axs[1,1].plot(x, h2_con, label ='Conjugate depth', linestyle ='--')
axs[1,1].axhline(hn, color='red', linestyle='--', label='Normal flow depth')
axs[1,1].axhline(hc, color='green', linestyle='--', label='Critical flow depth')
axs[1,1].set_xlabel('x [m]')
axs[1,1].set_ylabel('h [m]')
axs[1,1].legend(loc = 'lower left')
plt.savefig(save_path + 'backwater_curve_dam_damgate.png')
plt.show()
#%% hydraulic jump
d = np.where(h2_con - sol_1 > 0)[0]
x2_con = x_1[d[0]]
H2 = h2_con[d[0]]
x1_con = 20+ x2_con
H1 = conjugate_depth(H2, q, theta)
H_jump = H2-H1


L = length_hjump(H1, q, theta)
L1 = length_hjump_1(H2, theta)
#%% complete backwater curve
x_sol = np.linspace(0,20,1000)
bed = x_sol* -np.tan(theta)
bed_h1 = bed[:500]
bed_h2 = bed[500:]
x_h1 = np.linspace(0, x1_con, 500)
sol_h1 = solve_bresse(h0, x_h1, theta, K, q)
Sol_h1 = bed_h1 + sol_h1[:, 0]

x_h2 = np.linspace(0, x2_con, 500)
sol_h2 = solve_bresse(h0_1, x_h2, theta, K, q)
Sol_h2 = bed_h1 + sol_h2[:, 0]
x_h2_1 = 20 + x_h2

ellipse = Ellipse((x1_con, H1 + H_jump / 2), 0.1, H_jump, color='red', alpha=0.5)

plt.figure(figsize=(12, 6))
plt.plot(x_h1, sol_h1, label ='Water surface profile', color = 'blue')
plt.plot(x_h2_1, sol_h2, color='blue')
plt.fill_between(x_h1, sol_h1[:,0], color='blue', alpha=0.2)
plt.fill_between(x_h2_1, sol_h2[:,0], color='blue', alpha=0.2)
plt.axhline(hn, color='red', linestyle='--', label='Normal flow depth')
plt.axhline(hc, color='green', linestyle='--', label='Critical flow depth')
plt.gca().add_patch(ellipse)
plt.text(x1_con , H1 + H_jump / 2 + 0.06 , 'Hydraulic Jump', color='red', ha='center', va='center')
plt.xlabel('x[m]')
plt.ylabel('h[m]')
#plt.title('Backwater curve between gate and dam')
plt.legend()
plt.savefig(save_path + 'backwater_curve_combined.png')


#%% specific energy and pressure loss
def specific_energy(h, q):
    g = 9.81
    Es = h + (q**2)/(2*g*h**2)
    return Es

def pressure_drop(dH, H1, H2):
    dP = dH**3/(4*H1*H2)
    return dP

h = np.linspace(0.07, 0.15, 500)
Es = specific_energy(h, q)
Es_h1 = specific_energy(H1, q)
Es_h2 = specific_energy(H2, q)
Es_d = Es_h1-Es_h2

fig, ax = plt.subplots( figsize=(12, 8))

ax.plot(Es, h, label = 'specific energy')
ax.scatter(Es_h1, H1, label = 'supercritical')
ax.scatter(Es_h2, H2, label = 'subcritical')
ax.plot([Es_h1, Es_h2], [H1, H2], 'k--') 
ax.axvline(x=Es_h1, color='gray', linestyle='--')
ax.hlines(H2, Es_h1, Es_h2, colors = 'red', label ='Energy loss')
ax.legend()
ax.set_title('specific energy curve')
ax.set_xlabel('Es')
ax.set_ylabel('h')
plt.savefig(save_path + 'specific_energy.png')

dP = pressure_drop(H_jump, H1, H2)    
#%% changing slope

def slope(theta, K,q):
    hn = normal_flow_depth(theta, K, q)
    hc = critical_flow_depth(theta, q)
    if hn> hc:
        return 'mild slope'
    if hc> hn:
        return 'steep slope'
    if hc == hn:
        return 'critical slope'
    
def slopechange(theta_up, theta_do, K,q):
    slope_up = slope(theta_up, K, q)
    slope_do = slope(theta_do, K, q)
    if slope_up == 'steep slope' and slope_do == 'mild slope':
        return 'hydraulic jump'
    if slope_up == 'mild slope' and slope_do == 'steep slope':
        return 'no hydraulic jump, but flowregimechange'
    if slope_up == 'steep slope' and slope_do == 'steep slope':
        return 'flowregime stays supercritical'
    if slope_up == 'mild slope' and slope_do == 'mild slope':
        return 'flowregime stays subcritical'
    else:
        return 'critical slope'      

def slopechange_hjump(theta_up, theta_do, K, q, x_stop, steps):
    hn_up = normal_flow_depth(theta_up, K, q)
    H2_up = conjugate_depth(hn_up, q, theta_up)
    hn_do = normal_flow_depth(theta_do, K, q)
    H1_do = conjugate_depth(hn_do, q, theta_do)
    if H2_up < hn_do:
        '''hydraulic jump before slope change'''
        x = np.linspace(0,x_stop,steps)
        sol = solve_bresse(H2_up, x, theta_up, K, q)
        d = np.where(sol - hn_do > 0.01)[0]
        distance = x[d[0]]
        X = np.linspace(0, distance, steps)
        Sol = solve_bresse(H2_up, X, theta_up, K, q)
        return [Sol, distance, hn_do]
    if H1_do > hn_up:
        '''hydraulic jump after slope change'''
        x = np.linspace(0, -x_stop, steps)
        sol = solve_bresse(H1_do, x, theta_do, K, q)
        d = np.where(sol - hn_up < 0.01)[0]
        distance = x[d[0]]
        X = np.linspace(0, distance, steps)
        Sol = solve_bresse(H1_do, X, theta_do, K, q)
        return [Sol, distance, hn_up]
    else:
        return 'error'

def sub_to_supercritical_slopechange(theta_up, theta_do, K, q, x_stop, steps):
    #upstream
    hc_up = critical_flow_depth(theta_up, q)
    hn_up = normal_flow_depth(theta_up, K, q)
    start_up = hc_up * 1.001
    x_up = np.linspace(0, -x_stop/2, int(steps/2))
    sol_up = solve_bresse(start_up, x_up, theta_up, K, q)
    d_up = np.where(hn_up - sol_up < 0.01)[0]
    distance_up = x_up[d_up[0]]
    X_up = np.linspace(0, distance_up, int(steps/2))
    Sol_up = solve_bresse(start_up, X_up, theta_up, K, q)
    
    #downstream
    hc_do = critical_flow_depth(theta_do, q)
    hn_do = normal_flow_depth(theta_do, K, q)
    start_do = hc_do * 0.999
    x_do = np.linspace(0, x_stop, int(steps/2))
    sol_do = solve_bresse(start_do, x_do, theta_do, K, q)
    d_do = np.where(sol_do - hn_do < 0.01)[0]
    distance_do = x_do[d_do[0]]
    X_do = np.linspace(0, distance_do, int(steps/2))
    Sol_do = solve_bresse(start_do, X_do, theta_do, K, q)
    
    return [Sol_up, Sol_do, distance_up, distance_do]

    
def subcritical_slopechange(theta_up, theta_do, K, q, x_stop, steps):
    hn_up = normal_flow_depth(theta_up, K, q)
    hn_do = normal_flow_depth(theta_do, K, q)
    x_up = np.linspace(0, -x_stop, steps)
    sol_up = solve_bresse(hn_do, x_up, theta_up, K, q)
    if hn_up < hn_do:#M1 curve
        d = np.where(sol_up - hn_up < 0.01)[0]
        distance = x_up[d[0]]
        X = np.linspace(0, distance, steps)
        Sol_up = solve_bresse(hn_do, X, theta_up, K, q)
        return [sol_up,distance, hn_do]
    if hn_up > hn_do:#M2 curve
        d = np.where(hn_up - sol_up < 0.01)[0]
        distance = x_up[d[0]]
        X = np.linspace(0, distance, steps)
        Sol_up = solve_bresse(hn_do, X, theta_up, K, q)
        return [Sol_up, distance, hn_do]
    
    
def supercritical_slopechange(theta_up, theta_do, K, q, x_stop, steps):
    hn_up = normal_flow_depth(theta_up, K, q)
    hn_do = normal_flow_depth(theta_do, K, q)
    x_do = np.linspace(0, x_stop, steps)
    sol_do = solve_bresse(hn_up, x_do, theta_do, K, q)
    if hn_up < hn_do:
        d = np.where(hn_do - sol_do < 0.01)[0]
        distance = x_do[d[0]]
        X = np.linspace(0, distance, steps)
        Sol_do = solve_bresse(hn_up, X, theta_do, K, q)
        return [Sol_do,distance, hn_up]
    if hn_up > hn_do:
        d = np.where(sol_do - hn_do < 0.01)[0]
        distance = x_do[d[0]]
        X = np.linspace(0, distance, steps)
        Sol_do = solve_bresse(hn_up, X, theta_do, K, q)
        return [Sol_do, distance, hn_up]
        
   
def slopechange_backwatercurve(theta_up, theta_do, K, q, x_stop, steps):
    s_c = slopechange(theta_up, theta_do, K, q)
    if s_c == 'hydraulic jump':
        sol = slopechange_hjump(theta_up, theta_do, K, q, x_stop, steps)
    if s_c == 'no hydraulic jump, but flowregimechange':
        sol = sub_to_supercritical_slopechange(theta_up, theta_do, K, q, x_stop, steps)
    if s_c == 'flowregime stays supercritical':
        sol = supercritical_slopechange(theta_up, theta_do, K, q, x_stop, steps)
    if s_c == 'flowregime stays subcritical':
        sol = subcritical_slopechange(theta_up, theta_do, K, q, x_stop, steps)
    return sol
    
#%% Example: Hydraulic jump before slope change

Jf_up = 0.08
Jf_do = 0.005
theta_up = np.arctan(Jf_up)
theta_do = np.arctan(Jf_do)
K = 17
B = 15
Q = 150
q = Q/B
slopechange1 = slopechange(theta_up, theta_do, K, q)
curve = slopechange_backwatercurve(theta_up, theta_do, K, q, 100, 500)
sol_up_a = curve[0]
distance = curve[1]


x_up_b = np.linspace(0, 20, 500)
x_up_a = np.linspace(20, 20 + distance, 500)
r = 50-(20+distance)
x_do = np.linspace(0,r, 500 )
x_do_1 = np.linspace(20+distance, 50, 500)

hc_up = critical_flow_depth(theta_up, q)
hn_up = normal_flow_depth(theta_up, K, q)
hc_do = critical_flow_depth(theta_do, q)
hn_do = normal_flow_depth(theta_do, K, q)

sol_up_b = solve_bresse(hn_up, x_up_b, theta_up, K, q)
sol_do = solve_bresse(hn_do, x_do, theta_do, K, q)

h2_con_up = conjugate_depth(sol_up_b, q, theta_up)
h1_con_do = conjugate_depth(sol_do, q, theta_do)


#plot without slope
x_hjump = 20
H_jump = sol_up_a[0][0] - hn_up
ellipse = Ellipse((x_hjump, hn_up + H_jump / 2), 1.5, H_jump, color='red', alpha=0.5)


fig, ax = plt.subplots(figsize=(12, 6))
x = np.linspace(0,50, 500)


ax.plot(x_up_a, sol_up_a, label='Backwater curve', color='cyan', linewidth = 2.5)
ax.plot(x_up_b, sol_up_b, label ='Water surface at normal flow depth upstream', color='lightblue', linewidth = 2.5)
ax.plot(x_do_1, sol_do, label ='Water surface at normal flow depth downstream', color='darkblue', linewidth = 2.5)
ax.fill_between(x_up_b, sol_up_b[:,0], color='blue', alpha=0.2)
ax.fill_between(x_up_a,sol_up_a[:,0], color='blue', alpha=0.2)
ax.fill_between(x_do_1 ,sol_do[:,0], color='blue', alpha=0.2)
ax.plot(x, h2_con_up, label ='Conjugate depth upstream', linestyle ='--', color = 'darkgreen')
ax.plot(x, h1_con_do, label ='Conjugate downstream', linestyle ='--', color = 'yellow')
ax.hlines(hc_up, 0, 50, label = 'Critical depth', linestyle ='--', color = 'magenta')
#ax.hlines(hc_do, 0, 100, label = 'critical height', linestyle ='--', color = 'purple')
ax.add_patch(ellipse)
ax.text(x_hjump , hn_up + H_jump / 2 + 1 , 'Hydraulic Jump', color='red', ha='center', va='center')
ax.axvline(20+distance, linestyle = '--', color = 'darkblue')
ax.text(20+distance, 1.5 , 'Change of slope', color='darkblue', ha='center', va='center')
ax.legend(loc = 'lower right')
plt.savefig(save_path + 'slopechange_hjump_before_wo.png')

#ax.set_title('backwater curve changing slope')


# plot with slope
fig, ax = plt.subplots(figsize=(12,6))
x_up_b = np.linspace(0,20,500) 
x_up_a = np.linspace(20, 20+distance, 500)
x_do = np.linspace(20+distance, 100, 500)
bed_up_b = x_up_b * -np.tan(theta_up)
bed_up_a = bed_up_b[-1] - (x_up_a - x_up_a[0]) * np.tan(theta_up)
bed_do = bed_up_a[-1] - (x_do - x_do[0]) * np.tan(theta_do)
Sol_up_b = bed_up_b + sol_up_b[:, 0]
Sol_up_a = bed_up_a + sol_up_a[:, 0]
Sol_do = bed_do + sol_do[:,0]
x_hjump = 20
H_jump = h2_con_up[0][0] - hn_up
ellipse = Ellipse((x_hjump, Sol_up_b[-1] + H_jump / 2), 1.5, H_jump, color='red', alpha=0.5)


ax.plot(x_up_b, Sol_up_b, label ='hn upstream', color='lightblue', linewidth = 2.5)
ax.plot(x_up_a, Sol_up_a, label='S1 curve', color='cyan', linewidth = 2.5)
ax.plot(x_do, Sol_do, label ='hn downstream', color='darkblue', linewidth = 2.5)
ax.plot(x_up_b, bed_up_b, color = 'brown', label = 'slope')
ax.plot(x_up_a, bed_up_a, color = 'brown')
ax.plot(x_do, bed_do, color = 'brown')
ax.fill_between(x_up_b, bed_up_b, Sol_up_b, color='blue', alpha=0.2)
ax.fill_between(x_up_a, bed_up_a, Sol_up_a, color='blue', alpha=0.2)
ax.fill_between(x_do, bed_do, Sol_do, color='blue', alpha=0.2)
ax.add_patch(ellipse)
ax.text(x_hjump , Sol_up_b[-1] + H_jump / 2 + 1 , 'Hydraulic Jump', color='red', ha='center', va='center')
ax.set_xlabel('x [m]')
ax.set_ylabel('z [m]')
ax.set_title('backwater curve changing slope: hydraulic jump upstream')
ax.legend()
plt.savefig(save_path + 'slopechange_hjump_before_w.png')


#%% Example: hydraulic jump after slopechange
Jf_up = 0.2
Jf_do = 0.009
theta_up = np.arctan(Jf_up)
theta_do = np.arctan(Jf_do)
K = 17
B = 15
Q = 150
q = Q/B
slopechange2 = slopechange(theta_up, theta_do, K, q)

hn_up = normal_flow_depth(theta_up, K, q)
hn_do = normal_flow_depth(theta_do, K, q)
hc_up = critical_flow_depth(theta_up, q)
hc_do = critical_flow_depth(theta_do, q)
h2_con_up_zahl = conjugate_depth(hn_up, q, theta_up)
h1_con_do_zahl = conjugate_depth(hn_do, q, theta_do)

curve_1 = slopechange_backwatercurve(theta_up, theta_do, K, q, 100, 500)
sol_do_b = curve_1[0]
distance = curve_1[1]

x_up = np.linspace(0, 20, 500)
x_do_b = np.linspace(20 - distance, 20, 500)
r = 100- (20-distance)
x_do_a_1 = np.linspace(0,r, 500 )
x_do_a = np.linspace(20-distance, 100, 500)


sol_up = solve_bresse(hn_up, x_up, theta_up, K, q)
sol_do_a = solve_bresse(hn_do, x_do_a_1, theta_do, K, q)

h2_con_up = conjugate_depth(sol_up, q, theta_up)
h1_con_do = conjugate_depth(sol_do_a, q, theta_do)

x_hjump = 20 - distance
H_jump = hn_do - h1_con_do_zahl 
ellipse = Ellipse((x_hjump, h1_con_do_zahl + H_jump / 2), 1.5, H_jump, color='red', alpha=0.5)

#plot without slope
fig, ax = plt.subplots(figsize=(12, 6))
x = np.linspace(0,100, 500)

ax.plot(x_do_b, sol_do_b, label='Backwater curve', color='cyan', linewidth = 2.5)
ax.plot(x_up, sol_up, label ='Water surface at normal flow depth upstream', color='lightblue', linewidth = 2.5)
ax.plot(x_do_a, sol_do_a, label ='Water surface at normal flow depth downstream', color='darkblue', linewidth = 2.5)
ax.fill_between(x_up, sol_up[:,0], color='blue', alpha=0.2)
ax.fill_between(x_do_b,sol_do_b[:,0], color='blue', alpha=0.2)
ax.fill_between(x_do_a ,sol_do_a[:,0], color='blue', alpha=0.2)
ax.plot(x, h2_con_up, label ='Conjugate depth upstream', linestyle ='--', color = 'darkgreen')
ax.plot(x, h1_con_do, label ='Conjugate depth downstream', linestyle ='--', color = 'yellow')
ax.hlines(hc_up, 0, 100, label = 'Critical depth', linestyle ='--', color = 'magenta')
#ax.hlines(hc_do, 0, 100, label = 'critical height', linestyle ='--', color = 'purple')
ax.add_patch(ellipse)
ax.text(x_hjump , hn_up + H_jump / 2 + 1.4 , 'Hydraulic Jump', color='red', ha='center', va='center')
ax.text(20, 0.5 , 'Change of slope', color='darkblue', ha='center', va='center')
ax.axvline(20, linestyle = '--', color = 'darkblue')
ax.legend()
#ax.set_title('backwater curve changing slope: hydraulic jump downstream')
plt.savefig(save_path + 'slopechange_hjump_after_wo.png')


# plot with slope
fig, ax = plt.subplots(figsize=(12, 6))

bed_up = x_up * -np.tan(theta_up)
bed_do_b = bed_up[-1] - (x_do_b - x_do_b[-1]) * np.tan(theta_do)
bed_do_a = bed_do_b[-1] - (x_do_b - x_do_b[-1]) * np.tan(theta_do)
Sol_up = bed_up + sol_up[:, 0]
Sol_do_b = bed_do_b + sol_do_b[:, 0]
Sol_do_a = bed_do_a + sol_do_a[:,0]


ellipse = Ellipse((x_hjump, Sol_do_b[0]+ abs(Sol_do_b[0] - Sol_do_a[-1]) / 2), 1.5, H_jump, color='red', alpha=0.5)


ax.plot(x_up, Sol_up, label ='hn upstream', color='lightblue', linewidth = 2.5)
ax.plot(x_do_b, Sol_do_b, label='S1 curve', color='cyan', linewidth = 2.5)
ax.plot(x_do_a, Sol_do_a, label ='hn downstream', color='darkblue', linewidth = 2.5)
ax.plot(x_up, bed_up, color = 'brown', label = 'slope')
ax.plot(x_do_b, bed_do_b, color = 'brown')
ax.plot(x_do_a, bed_do_a, color = 'brown')
ax.fill_between(x_up, bed_up, Sol_up, color='blue', alpha=0.2)
ax.fill_between(x_do_b, bed_do_b, Sol_do_b, color='blue', alpha=0.2)
ax.fill_between(x_do_a, bed_do_a, Sol_do_a, color='blue', alpha=0.2)
ax.add_patch(ellipse)
ax.text(x_hjump , Sol_do_b[0] + H_jump / 2 + 1 , 'Hydraulic Jump', color='red', ha='center', va='center')
ax.text(20, 0.5 , 'Change of slope', color='darkblue', ha='center', va='center')
ax.axvline(20, linestyle = '--', color = 'darkblue')
ax.set_xlabel('x [m]')
ax.set_ylabel('z [m]')
ax.set_title('changing slope backwater curve')
plt.savefig(save_path + 'slopechange_hjump_after_w.png')

ax.legend()

#%% Example: mild slopechange
K = 17
B = 15
Q = 150
q = Q/B
Jf_up = 0.005
Jf_do = 0.003
theta_up = np.arctan(Jf_up)
theta_do = np.arctan(Jf_do)
slopechange4 = slopechange(theta_up, theta_do, K, q)
hn_up = normal_flow_depth(theta_up, K, q)
hc_up = critical_flow_depth(theta_up, q)
hn_do = normal_flow_depth(theta_do, K, q)
hc_do = critical_flow_depth(theta_do, q)


curve_3 = slopechange_backwatercurve(theta_up, theta_do, K, q, 800, 500)
sol_up = curve_3[0]
distance = curve_3[1]
x_do_1 = np.linspace(0, 2000+distance, 500)
sol_do = solve_bresse(hn_do, x_do_1, theta_do, K, q)
x_up = np.linspace(-distance,0, 500)
x_do = np.linspace(-distance, 2000, 500)
x = np.linspace(0,2000,500)

#plot without slope
fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(x_up, sol_up, label ='back water curve', color='lightblue', linewidth = 2.5)
ax.plot(x_do, sol_do, label ='Water surface at normal flow depth downstream', color='darkblue', linewidth = 2.5)
ax.hlines(hn_up, 0, -distance,label = 'Normal flow depth upstream', linestyle ='--', color = 'purple')
ax.fill_between(x_up, sol_up[:,0], color='blue', alpha=0.2)
ax.fill_between(x_do,sol_do[:,0], color='blue', alpha=0.2)
ax.hlines(hc_up, 0, 2000, label = 'Critical depth', linestyle ='--', color = 'magenta')
ax.hlines(hc_do, 0, 2000, linestyle ='--', color = 'magenta')
ax.text(-distance, 0.5 , 'Change of slope', color='darkblue', ha='center', va='center')
ax.axvline(-distance, linestyle = '--', color = 'darkblue')
ax.legend()
#ax.set_title('backwater curve mild slopechange')
plt.savefig(save_path + 'slopechange_mild_wo.png')


#%% Example: steep slopechange
K = 17
B = 15
Q = 150
q = Q/B
Jf_up = 0.08
Jf_do = 0.15
theta_up = np.arctan(Jf_up)
theta_do = np.arctan(Jf_do)
slopechange5 = slopechange(theta_up, theta_do, K, q)
hn_up = normal_flow_depth(theta_up, K, q)
hc_up = critical_flow_depth(theta_up, q)
hn_do = normal_flow_depth(theta_do, K, q)
hc_do = critical_flow_depth(theta_do, q)

curve_4 = slopechange_backwatercurve(theta_up, theta_do, K, q, 100, 500)
sol_do = curve_4[0]
distance = curve_4[1]

x_up = np.linspace(0, 20, 500)
x_do = np.linspace(20, 20 + distance, 500)
sol_up = solve_bresse(hn_up, x_up, theta_up, K, q)


fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(x_do, sol_do, label ='', color='darkblue', linewidth = 2.5)
ax.plot(x_up, sol_up, label ='Water surface at normal flow depth upstream', color='lightblue', linewidth = 2.5)
ax.hlines(hn_do, 20, 20+distance,label = 'Normal flow depth downstream', linestyle ='--', color = 'purple')
ax.fill_between(x_up, sol_up[:,0], color='blue', alpha=0.2)
ax.fill_between(x_do,sol_do[:,0], color='blue', alpha=0.2)
ax.hlines(hc_up, 0, 20+distance, label = 'Critical depth', linestyle ='--', color = 'magenta')
ax.hlines(hc_do, 0, 20+distance, linestyle ='--', color = 'magenta')
ax.text(20, 0.5 , 'Change of slope', color='darkblue', ha='center', va='center')
ax.axvline(20, linestyle = '--', color = 'darkblue')
ax.legend(loc ='lower right')
#ax.set_title('backwater curve steep slopechange')
plt.savefig(save_path + 'slopechange_steep_wo.png')

#%% Example: slopechange with flow regime change from fluvial to torrential 
Jf_up = 0.005
Jf_do = 0.07
theta_up = np.arctan(Jf_up)
theta_do = np.arctan(Jf_do)
K = 17
B = 15
Q = 150
q = Q/B
slopechange6 = slopechange(theta_up, theta_do, K, q)
x_up = np.linspace(0,-640,500)
x_do = np.linspace(0, 100, 500)

hn_up = normal_flow_depth(theta_up, K, q)
hc_up = critical_flow_depth(theta_up, q)
hc_do = critical_flow_depth(theta_do, q)
hn_do = normal_flow_depth(theta_do, K, q)

sol_up = solve_bresse(2.17, x_up, theta_up, K, q)
sol_do = solve_bresse(2.17, x_do, theta_do, K, q)

#plot without slope
fig, ax = plt.subplots(figsize=(12, 6))
ax.hlines(hn_up, 0,-1000, label ='hn up', linestyle = '--')
ax.hlines(hn_do, 0,100, label ='hn up', linestyle = '--')
ax.hlines(hc_up, 0, -100, label = 'critical height', linestyle ='--', color = 'magenta')
ax.hlines(hc_do, 0, 100, linestyle ='--', color = 'magenta')
ax.legend()

Sol = slopechange_backwatercurve(theta_up, theta_do, K, q, 2000, 500)
sol_up = Sol[0]
sol_do = Sol[1]
d_up = Sol[2]
d_do = Sol[3]
fig, ax = plt.subplots(figsize=(12, 6))
x_up = np.linspace(0,d_up, 250)
x_do = np.linspace(0, d_do, 250)
ax.plot(x_up, Sol[0], label = 'Backwater curve upstream')
ax.plot(x_do, Sol[1], label = 'Backwater curve downstream')
ax.hlines(hn_up, d_up,0, label ='Normal flow depth upstream', linestyle = '--', color = 'green')
ax.hlines(hn_do, 0, d_do, label = 'Normal flow depth downstream', linestyle ='--')
ax.hlines(hc_up, d_up, d_do, label = 'Critical depth', linestyle ='--', color = 'magenta')
ax.hlines(hc_do, d_up, d_do, linestyle ='--', color = 'magenta')
ax.axvline(0, linestyle = '--', color = 'darkblue')
ax.text(0, 3 , 'Change of slope', color='darkblue', ha='center', va='center')
ax.fill_between(x_up, Sol[0][:,0], color='blue', alpha=0.2)
ax.fill_between(x_do, Sol[1][:,0], color='blue', alpha=0.2)
ax.legend(loc = 'lower right')
plt.savefig(save_path + 'slopechange_fchange_wo.png')


#%% plot with slope
fig, ax = plt.subplots(figsize=(12, 6))

bed_up = x_up * -np.tan(theta_up)
bed_do = - x_do * np.tan(theta_do)
Sol_up_1 = bed_up + sol_up[:, 0]
Sol_do_1 = bed_do + sol_do[:, 0]


ax.plot(x_up, Sol_up_1, label ='subcritical backwater curve', color='lightblue', linewidth = 2.5)
ax.plot(x_do, Sol_do_1, label = 'supercritical backwater curve', color='cyan', linewidth = 2.5)
ax.plot(x_up, bed_up, color = 'brown', label = 'slope')
ax.plot(x_do, bed_do, color = 'brown')
ax.fill_between(x_up, bed_up, Sol_up_1, color='blue', alpha=0.2)
ax.fill_between(x_do, bed_do, Sol_do_1, color='blue', alpha=0.2)
ax.text(0, 0.5 , 'change of slope', color='darkblue', ha='center', va='center')
ax.axvline(0, linestyle = '--', color = 'darkblue')
ax.set_xlabel('x [m]')
ax.set_ylabel('z [m]')
ax.set_title('changing slope backwater curve')
ax.legend()
plt.savefig(save_path + 'slopechange_fchange_w.png')


