import numpy as np
import matplotlib.pyplot as plt

def partition_function(T,epsilon,g):
    Z_list = []
    for i in T:
        Z = 0
        for j,m in zip(epsilon,g):
            Z = Z + m * np.exp(-j/(k*i))
        Z_list.append(Z)
    return np.array(Z_list)

def fraction_population(g,T,epsilon,Z):
    frac_pop_list = []
    for j in range(len(epsilon)):                                                                                             
        frac_pop = (g[j] * np.exp(-epsilon[j]/(k*T)))/Z
        frac_pop_list.append(frac_pop)
    frac_pop_list = np.array(frac_pop_list)
    return frac_pop_list

def internal_energy(frac_pop,N,epsilon,T):
    N_j = N * frac_pop
    inte_energy = np.zeros([len(T)])
    for i in range(len(N_j)):
        inte_energy = inte_energy + N_j[i]*epsilon[i]
    return inte_energy

def entropy(Z,N,T,U):
    S = (N*k*np.log(Z/N)) + (U/T) + (N*k)
    return S

def free_energy(N,T,Z):
    F = -N*k*T*np.log(Z)
    return F

def graph(x1,x2,y1,y2,title,y_label,frac_pop_low,frac_pop_high,key):
    fig1,ax1 = plt.subplots(1,2)
    fig1.suptitle(title)
    if key == 0:
        ax1[0].scatter(x1,y1,label = "Low Temperature")
        ax1[0].set_xlabel("T")
        ax1[0].set_ylabel(y_label)
        ax1[0].grid(ls = "--")
        ax1[0].legend()
        ax1[1].scatter(x2,y2,label = "High Temperature")
        ax1[1].set_xlabel("T")
        ax1[1].set_ylabel(y_label)
        ax1[1].grid(ls = "--")
        ax1[1].legend()
        plt.show()
    elif key == 1:
        for i in range(len(frac_pop_low)):
            ax1[0].scatter(x1,frac_pop_low[i],label = "LowTemperature")
            ax1[1].scatter(x2,frac_pop_high[i],label = "High Temperature")
        ax1[0].set_xlabel("T")
        ax1[0].set_ylabel("$\\dfrac{N_i}{N}$")
        ax1[1].set_xlabel("T")
        ax1[1].set_ylabel("$\\dfrac{N_i}{N}$")
        ax1[0].grid(ls = "--")
        ax1[1].grid(ls = "--")
        ax1[0].legend()
        ax1[1].legend()
        plt.show()
                    
if __name__ == "__main__":
    k = 8.617 * 10**(-5)
    
    T_low = np.linspace(10**(-18),5000,50)
    T_high = np.linspace(5000,10**(5),50)

    g = [1,1] ; epsilon = [0,1]

    Z_1 = partition_function(T_low,epsilon,g)
    Z_2 = partition_function(T_high,epsilon,g)
    
    frac_pop_low = fraction_population(g,T_low,epsilon,Z_1)
    frac_pop_high = fraction_population(g,T_high,epsilon,Z_2)
    
    U_low = internal_energy(frac_pop_low,1,epsilon,T_low)
    U_high = internal_energy(frac_pop_high,1,epsilon,T_high)
    
    S_low = entropy(Z_1,1,T_low,U_low)
    S_high = entropy(Z_2,1,T_high,U_high)
    
    F_low = free_energy(1,T_low,Z_1)
    F_high = free_energy(1,T_high,Z_2)
    
    graph(T_low,T_high,Z_1,Z_2,"Partition Function","Z",None,None,0)
    graph(T_low,T_high,None,None,"Fraction Population",None,frac_pop_low,frac_pop_high,1)
    graph(T_low,T_high,U_low,U_high,"Internal Energy","U",None,None,0)
    graph(T_low,T_high,S_low,S_high,"Entropy","S",None,None,0)
    graph(T_low,T_high,F_low,F_high,"Helmholtz free energy","F",None,None,0)
