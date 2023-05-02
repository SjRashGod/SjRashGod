import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

def F_rj(x):
    return (x**2)

def F_p(x):
    return x**3/(np.exp(x)-1)

x=np.linspace(0,12,151)
plt.scatter(x * np.pi,F_rj(x),label='Rayleigh Jeans')
plt.grid(ls='--')
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel("x")
plt.ylabel("G(x)")
plt.title("Density of States(Rayleigh Jeans)")
plt.legend()
plt.show()

plt.scatter(x * np.pi,F_p(x),label='Planck')
plt.grid(ls='--')
plt.xlabel("x")
plt.ylabel("G(x)")
#plt.xscale("log")
#plt.yscale("log")
plt.legend()
plt.title("Density of States (Planck)") 
plt.show()

h=6.62e-34
c=3e8
temp=[1200,1500,1800]
k=1.38*10**(-23)
l=1e-10

def  Ray(x,t):
    const = 8 * ((k)**4)/((h*c)**3)
    return (np.pi)*((t)**4)*const*F_rj(x)

for i in temp:
        plt.scatter(x*(k*i)/(2*h),Ray(x,i)*(2*h/k*i),marker='.',label='For Temp. ' + str(i) + 'K')
plt.grid(ls='--')
plt.xlabel("v")
plt.ylabel("U(v)")
plt.legend()
plt.title("Rayleigh Jeans ")
plt.show()


def Planck(x,t):
    const = 8*((k)**4)/((h*c)**3)    
    return (np.pi)*((t)**4)*const*F_p(x)

for i in temp:
        plt.scatter(x*(k*i)/(2*h),Planck(x,i)*(2*h/k*i),marker='.',label='For Temp. ' + str(i) + 'K')
plt.grid(ls='--')
plt.xlabel("v")
plt.ylabel("U(v)")
plt.legend()
plt.title("Planck ")
plt.show()

#Plotting with Solar Temp
sol=[6000]
for i in sol:
        plt.scatter(x*(k*i)/(2*h),Ray(x,i)*(2*h/k*i),marker='.',label='For Temp. ' + str(i) + 'K')
plt.grid(ls='--')
plt.xlabel("v")
plt.ylabel("U(v)")
plt.legend()
plt.title("Rayleigh Jeans (Solar Temp)")
plt.show()

for i in sol:
        plt.scatter(x*(k*i)/(2*h),Planck(x,i)*(2*h/k*i),marker='.',label='For Temp. ' + str(i) + 'K')
plt.grid(ls='--')
plt.xlabel("v")
plt.ylabel("U(v)")
plt.legend()
plt.title("Planck (Solar Temp)")
plt.show()

def density(v):
    return (8*np.pi)/(c**3)* (v**2)*(l**3)

v_complete=np.logspace(10,30,3000)

plt.scatter(v_complete,density(v_complete),marker='.')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Frequency(v)")
plt.ylabel("Density of states")
plt.title("Density of states with Frequency")
plt.grid(ls='--')
plt.show()

