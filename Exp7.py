import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import linregress

h = 6.63*10**(-34)
k = 1.38*10**(-23)
N = 6.023*10**(23)
m = 1.67*10**(-27)
V = np.linspace(20*10**(-3),50*10**(-3),10)
T = np.linspace(150,450,10)

#Partition Function 

def Z_VT(n,V,T):
    const = h**2/(8*m*(V**(2/3))*k*T)
    return(np.pi/2)*(n**2)*np.exp((-const*(n**2)))
                         
def inte(n_i,n_f,V,T):
    I = quad(lambda n: Z_VT(n,V,T),n_i,n_f)
    return I

part_fn = np.zeros((len(V),len(T)))
for i in range(len(V)):
    for j in range(len(T)):
        part_fn[i,j] = inte(0,10**(11),V[i],T[j])[0]

part_fn_log = np.zeros((len(V),len(T)))
for i in range(len(V)):
    for j in range(len(T)):
        part_fn_log[i,j] = np.log(inte(0,10**(11),V[i],T[j]))[0]
        
plt.plot(np.log(V),part_fn_log[0,:])
plt.xlabel("Log V")
plt.ylabel("Log Z")
plt.title("Partition Function (at constant Temperature)")
plt.grid(ls="--")
plt.show()
    
plt.plot(np.log(T),part_fn_log[:,0])
plt.xlabel("Log T")
plt.ylabel("Log Z")
plt.title("Partition Function (at constant Volume)")
plt.grid(ls="--")
plt.show()

#Pressure

def pressure(V,T):
    for i in range(len(T)):
        Diff = ((part_fn_log[i,:][:-1]) - (part_fn_log[i,:][1:]))/(V[:-1]- V[1:])
        p = N*k*T[i]*Diff
        plt.plot(V[:-1],p,label='For T= '+str(T[i]))
        #return p

pressure(V,T)
plt.legend()
plt.xlabel("Volume")
plt.ylabel("Pressure")
plt.title("Pressure (at constant Temperature)")
plt.grid(ls="--")
plt.show()

#Internal Energy

def av_energy(T):
    E=[]
    for i in range(len(V)):
        Diff = ((part_fn_log[i,:][:-1]) - (part_fn_log[i,:][1:]))/(T[:-1]- T[1:])
        e = k*T[:-1]**2 *Diff
        E.append(e)
        plt.legend()
        plt.xlabel("Temperature")
        plt.ylabel("Energy")
        plt.title("Energy")
        plt.grid(ls="--")
        plt.plot(T[:-1],e,label='For V= '+str(V[i]))
    plt.show()
    return E
    
U = N*av_energy(T)[1]
plt.plot(T[:-1],U)
plt.xlabel("Temperature")
plt.ylabel("Energy")
plt.title("Internal Energy")
plt.grid(ls="--")
plt.show()

#Specific Heat

slope = linregress(T[:-1],U/N)[0]
print('Specific Heat is:' ,slope)

#Entropy

def entropy(U,T):
    for i in range(len(T)):
        Diff = ((part_fn_log[i,:][:-1]) - np.log(N)) #(part_fn_log[i,:][:-1]) - (part_fn_log[i,:][1:]))/(V[:-1]- V[1:])
        S = (U/T[:-1]) + (N*k*(Diff + 1))
        plt.plot(V[:-1],S,label='For T= '+str(T[i]))
        #return p
        

entropy(U,T)
plt.legend()
plt.xlabel("Volume")
plt.ylabel("Entropy")
plt.title("Entropy(at constant Temperature)")
plt.grid(ls="--")
plt.show()

