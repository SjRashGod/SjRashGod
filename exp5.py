import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.integrate import quad
from scipy.stats import linregress

def planck(x):
    return x**3/(np.exp(x) - 1)

def median_point(tol):
    area = simps(u,x)
    area_h = area/2
    for i in range(2,len(x)):
        area1 = simps((u[0:i]),(x[0:i]))
        
        if np.abs(area1 - area_h) <= tol:
            break
    
    return x[i]

def u_d(T):
    return (np.pi**4/15)*(8*np.pi*(k*T)**4/(h*c)**3)       

if __name__ == "__main__":
    k = 1.38 * 10**(-23)
    h = 6.624 * 10**(-34)
    c = 3 * 10**(8)
    
    x_in = 1e-2 ; x_fin = 12
    x = np.linspace(x_in,x_fin,5000)
    
    u = planck(x)
    
    median = median_point(1e-3)
    
    print('\nx_median:',median)

    # WEIN'S CONSTANT

    b = h*c/(k*median)
    print("\nWein's constant:",b)
    
    # STEFAN
    
    I_p = quad(planck,1e-15,20)
    
    print("\nValue of I_p is:",I_p[0],"and Standard value of I_p is:",np.pi**4/15)
    
    T = np.arange(100,10200,500)
    
    F = lambda T : (c/4)*u_d(T)
    
    res = linregress(np.log(T),np.log(F(T)))

    print('\nSlope:',res[0])
    print('\nIntercept:',res[1])

    print('\nStefan constant:',np.exp(res[1]),"\n")

    fig,ax = plt.subplots()
    ax.plot(x,u)
    ax.set_title('Planck law of radiation')
    ax.set_xlabel('x')
    ax.set_ylabel('Energy spectral density')
    plt.grid(ls = "--")
    plt.show()

    fig,ax = plt.subplots()
    ax.plot(T,F(T))
    ax.scatter(T,F(T))
    ax.set_title('Radiant flux (F) Vs Temperature (T)')
    ax.set_xlabel('Temperature (T)')
    ax.set_ylabel('Radiant flux (F)')
    plt.grid(ls = "--")
    plt.show()

    fig,ax = plt.subplots()
    ax.plot(np.log(T),np.log(F(T)))
    ax.scatter(np.log(T),np.log(F(T)))
    ax.set_title("Linear Regression")
    ax.set_xlabel('log(T)')
    ax.set_ylabel('log(F)')
    plt.grid(ls = "--")
    plt.show()
    
    