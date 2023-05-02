import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int1
from scipy.stats import linregress

def F_p(x):
    return (x**3)/(np.exp(x)-1)

x=np.linspace(1e-8,12,151)
X=np.array(x)
P=np.array(F_p(x))
m =  np.max(F_p(x))
t = F_p(x)
i_val = np.where(t == m)

plt.scatter(x,F_p(x))
plt.grid(ls='--')
plt.xlabel("x")
plt.ylabel("U(x)")
plt.title("Spectral Energy Density")
plt.show()

Area=int1.quad(F_p,0,12)[0]
print("1.(a) Value of x_median is: \n" , Area/2)

h=6.62e-34
c=3e8
k=1.38*10**(-23)
T = 5800
b = (h*c)/(k*Area/2)

print("1.(b) Value of Wien's Displacement Constant,b is: \n",h*c/(k*Area/2), "metres Kelvin")
print("Hence, Wien's Displacement Law has been verified \n")
print("Corresponding Wavelength is: ", b/T, "metres")
      
I_1=int1.simps(P,X)
print("2. The value of integral I_p is: \n",I_1)

test = ((np.pi)**4)/15
print("Value of RHS: \n",test)
print("Hence, I_p = pi**4/15 is shown \n")

def U(T):
    return (np.pi**4)/(15) * 8*(np.pi)*((k*T)**4/(h*c)**3)

T=np.arange(100,10000,500)
values=[]
for i in T:
    value=U(i)*(15/(np.pi)**4)
    values.append(value)

def radiant_flux(T,C):
    return (c/4)*((np.pi)**4/(15)) *C

for i in T:
    value=U(i)*(15/(np.pi)**4)
    plt.scatter(i,radiant_flux(i,value) , label='T =' + str(i))
plt.xlabel("Temperature")
plt.ylabel("F(T)")
plt.title("2.(a) Plot of Radiant Flux vs Temperature")
plt.grid(ls="--")
plt.show()

temp=[]
for i in T:
    value=U(i)*(15/(np.pi)**4)
    plt.scatter(np.log(i),np.log(radiant_flux(i,value)))
    temp.append(radiant_flux(i,value))
slope=linregress(np.log(T) , np.log(temp))
print("slope",slope[0])
print("intercept",(slope[1]))
plt.xlabel("Temperature")
plt.ylabel("F(T)")
plt.title("2.(b) Log plot of Radiant Flux vs Temperature")
plt.grid(ls='--')
plt.show()
print ("Value of Stefan Boltzman law constant is: \n" , np.exp(slope[1]), "Joules K^-4 m^-2 s^-1")
