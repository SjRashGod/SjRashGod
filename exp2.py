import numpy as np
import matplotlib.pyplot as plt

def max_bolt(x):
    return np.exp(-x)

def bose_einstein(x,alpha):
    return 1/(np.exp(x+alpha) - 1)

def fermi_dirac(x,alpha):
    return 1/(np.exp(x+alpha) + 1)


if __name__ == "__main__":
    x_range_max = np.linspace(-4,4,50)
    x_range_bose = np.linspace(0.1,4,50)
    x_range_fermi = np.linspace(-4,4,50)
    alpha = [0,1]
    k = 8.617333 * 10**(-5)
    T = np.array([10,100,1000,5000])
    fermi_total_fermi = [] ; fermi_total_bose = [] ; fermi_total_max = []
    
    f_max_bolt = max_bolt(x_range_max)
    f_bose_einstein = bose_einstein(x_range_bose,alpha[0])
    f_fermi_dirac = fermi_dirac(x_range_fermi,alpha[1])

    fig, ax1 = plt.subplots(1, 3, figsize=(10, 4))

    ax1[0].plot(x_range_max,f_max_bolt)
    ax1[0].scatter(x_range_max,f_max_bolt,marker = ".")
    ax1[1].plot(x_range_bose,f_bose_einstein,c = "g")
    ax1[1].scatter(x_range_bose,f_bose_einstein,c = "g",marker = ".")
    ax1[2].plot(x_range_fermi,f_fermi_dirac,c = "r")
    ax1[2].scatter(x_range_fermi,f_fermi_dirac,c = "r",marker = ".")
    for i in range(3):
        ax1[i].set(xlabel = "$\epsilon$/KT",ylabel = "f($\epsilon$)")
        ax1[i].grid(ls = "--")
    ax1[0].set(title = "Maxwell Boltzman Distribution")
    ax1[1].set(title = "Bose Einstein Distribution")
    ax1[2].set(title = "Fermi Dirac Distribution")
    plt.show()

    for i in range(len(T)):
        fermi_x = x_range_fermi * T[i]*k
        fermi_total_fermi.append(fermi_x)

    for i in range(len(T)):
        fermi_x = x_range_max * T[i]*k
        fermi_total_max.append(fermi_x)
        
    for i in range(len(T)):
        fermi_x = x_range_bose * T[i]*k
        fermi_total_bose.append(fermi_x)

    
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()

    for i in range(len(fermi_total_fermi)):
        ax2.plot(fermi_total_fermi[i],f_fermi_dirac,label = "At T = "+str(T[i])+" K")
        ax2.scatter(fermi_total_fermi[i],f_fermi_dirac,marker = ".")
        ax3.plot(fermi_total_max[i],f_max_bolt,label = "At T = "+str(T[i])+" K")
        ax3.scatter(fermi_total_max[i],f_max_bolt,marker = ".")
        ax4.plot(fermi_total_bose[i],f_bose_einstein,label = "At T = "+str(T[i])+" K")
        ax4.scatter(fermi_total_bose[i],f_bose_einstein,marker = ".")
    ax2.set(xlabel = "$\epsilon$",ylabel = "f($\epsilon$)",title = "Fermi Dirac distribution at constant temperature")
    ax3.set(xlabel = "$\epsilon$",ylabel = "f($\epsilon$)",title = "Maxwell Boltzman distribution at constant temperature")
    ax4.set(xlabel = "$\epsilon$",ylabel = "f($\epsilon$)",title = "Bose Einstein distribution at constant temperature")
    ax2.grid(ls = "--")
    ax3.grid(ls = "--")
    ax4.grid(ls = "--")
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()

    
