import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import diff
from scipy.stats import linregress

def partition_function(V,T):

    def integral(n):
        const = h**2/(8*m*(V**(2/3))*k_b*T)
        return n**(2)* np.exp(-n**2 * const)

    return np.pi/2 * quad(integral,0,10**(11))[0]

def result_matrix(V_array,T_array):
    result = np.zeros((len(V_array),len(T_array)))
    result_log = np.zeros((len(V_array),len(T_array)))

    for i in range(len(V_array)):
        for j in range(len(T_array)):
            result[i][j] = partition_function(V_array[i],T_array[j])
            result_log[i][j] = np.log(partition_function(V_array[i],T_array[j]))
    return result.T,result_log.T

def pressure(V,T,Z):
    Z_diff = diff(Z)
    V_diff = diff(V)
    diff_term = Z_diff/V_diff
    p = N_a*k_b*T*diff_term
    return p

def energy(Z,T):
    Z_diff = diff(Z)
    T_diff = diff(T)
    diff_term = Z_diff/T_diff
    E = k_b * T[:-1]**2 * diff_term
    return E,E*N_a

def entropy(U,T,Z):
    term = N_a * k_b * (Z - np.log(N_a) + 1)
    en = U/T + term
    return en

def energy_fluctation(T,C_v):
    e_fluc = k_b * T**2 * C_v
    return e_fluc
    
if __name__ == "__main__":
    k_b = 1.38 * 10**(-23)
    N_a = 6.022 * 10**(23)
    m = 3.32 * 10**(-27)
    h = 6.63 * 10**(-34)

    V_array = np.linspace(20 * 10**(-3),50 * 10**(-3),10)
    T_array = np.linspace(150,450,10)
    p_array = [] ; e_array = [] ; u_array = [] ; en_array = []

    result,result_log = result_matrix(V_array,T_array)

    for i in range(len(V_array)):
        plt.plot(V_array,result_log[i,:],label = "At T = "+str(T_array[i]))
    plt.xlabel("V")
    plt.ylabel("log(Z)")
    plt.title("Partition Function (at constant temperature)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()

    for i in range(len(V_array)):
        plt.plot(T_array,result_log[:,i],label = "At V = "+str(V_array[i]))
    plt.xlabel("T")
    plt.ylabel("log(Z)")
    plt.title("Partition Function (at constant volume)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()

    # PRESSURE

    for i in range(len(T_array)):
        p = pressure(V_array,T_array[i],result_log[i,:])
        p_array.append(p)
        plt.plot(V_array[:-1],p,label = "At T = "+str(T_array[i]))
    plt.xlabel("Volume (meter^3)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Pressure Vs Volume (At Constant Temperature)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()
        
    p_array = np.array(p_array)
    
    for i in range(len(V_array[:-1])):
        plt.plot(T_array[:-1],p_array[:,i][:-1],label = "At V = "+str(V_array[i]))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Pressure Vs Tempertaure (At Constant Volume)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()
        
    # ENERGY
  
    for i in range(len(T_array)):
        E,U = energy(result_log[:,i],T_array)
        e_array.append(E)
        u_array.append(U)
        #print(len(U))
        plt.plot(T_array[:-1],E)
    plt.xlabel("Temperature")
    plt.ylabel("Energy")
    plt.title("Energy Vs temperature (At constant Volume)")
    plt.grid(ls = "--")
    plt.show()

    slope = linregress(T_array[:-1],e_array[0])[0]
    print("\nCv :",slope)
        
    u_array = np.array(u_array)

    # ENTROPY

    for i in range(len(T_array[:-1])):
        en = entropy(u_array[i,:],T_array[:-1],result_log[:,i][:-1])
        en_array.append(en)
        plt.plot(T_array[:-1],en,label = "At V = "+str(V_array[i]))
    plt.xlabel("Temperature")
    plt.ylabel("Entropy")
    plt.title("Entropy Vs temperature (At constant Volume)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()

    # ENERGY FLUCTATION

    e_fluc = energy_fluctation(T_array,slope)
    print("\nEnergy Fluctation :",e_fluc)

    plt.plot(T_array,e_fluc)
    plt.xlabel("Temperature")
    plt.ylabel("Variance")
    plt.title("Variance Vs temperature")
    plt.grid(ls = "--")
    plt.show()

    
