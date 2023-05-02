import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

theta_e = 230
theta_d = 343
r = theta_e/theta_d

nu = np.linspace(0,2,100)

def dulong_petit(x):
    return np.ones(len(x))

def einstein(x):
    return ((1/x)**2) * ((np.exp(1/x))/(np.exp(1/x) - 1)**2)

einstein = np.vectorize(einstein)

def debye(x):
    u1 = -3*x/(np.exp(x)-1)
    u2 = 12/x**3
    ine = quad(lambda x : (x**3/(np.exp(x)-1)),0,x)
    return u1 + u2*ine[0]

debye = np.vectorize(debye)

def dos_einstein(nu):
    if abs(nu-1) <= 0.02:
       return 1
    else:
       return 0

def dos_debye(nu):
    if nu <= 1/r:
       return nu**2
    else:
       return 0

dos_einstein = np.vectorize(dos_einstein)
dos_debye = np.vectorize(dos_debye)
    
"Einstein distribution Law"
    
Temp_einstein = np.linspace(1e-8,2*theta_e,100)/theta_e

Specific_einstein = einstein(Temp_einstein)
    
"Debye distribution Law"

Temp_debye = np.linspace(1e-8,2*theta_d,100)/theta_d

Specific_debye = debye(1/Temp_debye)

"Dulong-Petit distribution Law"
    
Temp_dulong = np.linspace(Temp_debye[0],Temp_debye[-1],100) 

Specific_dulong = dulong_petit(Temp_dulong)
    
# Figure 1: Density of States 

fig,ax = plt.subplots(figsize = (16,8))
    
ax.plot(nu,dos_einstein(nu),label = 'Einstein', c = 'b')
ax.plot(nu,dos_debye(nu)/dos_debye(1/r),label = 'Debye', c = 'r')
ax.scatter(nu,dos_debye(nu)/dos_debye(1/r),marker = ".",c = "r")
ax.set_ylim([0,1.1])
ax.set(xlabel = '$\\nu$/$\\nu_x$',title = "Density of States",ylabel = "No. of States")
ax.legend()
ax.grid(ls = "--")
plt.show()
    
# Figure 2: Specific Heat Plot

fig1,ax1 = plt.subplots(figsize = (16,8))

ax1.plot(Temp_dulong,Specific_dulong,label = "Dulong-Petit distribution Law",ls = "--", c = 'r')
ax1.plot(Temp_einstein,Specific_einstein,label = "Einstein distribution Law",c = "b")
ax1.plot(Temp_debye,Specific_debye,label = "Debye distribution Law",c = "g")
ax1.scatter(Temp_einstein,Specific_einstein,marker = ".",c = "b")
ax1.scatter(Temp_debye,Specific_debye,marker = ".",c = "g")
ax1.set(xlabel = r"T/$\theta$",ylabel = "$\dfrac{C_v}{3R}$",title = "Specific Heat of Solids")
ax1.legend()
ax1.grid(ls = "--")
plt.show()
    
    
