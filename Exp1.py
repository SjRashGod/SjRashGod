import numpy as np
import matplotlib.pyplot as plt
import random

def macrostates(N_t,N_c):
    result = []
    for i in range(N_t): # (N_t)
            outcomes = []
            for j in range(N_c):
                out = random.randint(0,1)
                outcomes.append(out)
            result.append(outcomes)
    return result

def macrostates_counter(result,N_c):
    result = np.array(result)
    total = [] ; hc_array  = []
    for i in range(len(result)):
        macro = []
        h = sum(result[i])
        hc_array.append(h)
        t = N_c - h
        macro.append(h)
        macro.append(t)
        total.append(macro)
    return total,hc_array

def head_count(hc_array,N_c):
    pass
    freq_h = []
    for i in range(N_c+1):
        f_h = hc_array.count(i)
        freq_h.append(f_h)
    return np.array(freq_h)
    
if __name__ == "__main__":
    # Random number generator
    # 1 = Heads
    # 0 = Tails
    
    N_c = 3
    N_t = [10,100,1000,10000]
    
    prob_nc = [] ; num_of_heads_nc = []
    prob_nt = [] ; num_of_heads_nt = []
    
    for k in N_t:
    
        num_of_heads = [i for i in range(0,N_c+1)]
        num_of_heads_nt.append(num_of_heads)
        
        result = macrostates(k,N_c)
        total,hc_array = macrostates_counter(result,N_c)
        freq_h = head_count(hc_array,N_c)
        probability = freq_h/k
        prob_nt.append(probability)
        
    fig1,ax1 = plt.subplots()
        
    for i in range(len(prob_nt)):
        ax1.plot(num_of_heads_nt[i],prob_nt[i],label = "No of trails = "+str(N_t[i]))
        ax1.scatter(num_of_heads_nt[i],prob_nt[i])
    ax1.set(xlabel = "No of Heads",ylabel = "Probabilty of macrostates",title = "Probility of macrostates V/s No. of Heads at constant no. of coins ("+str(N_c)+")")
    ax1.grid(ls = "--")
    ax1.legend()
    plt.show()
        
    N_c = [1,2,3,4,5,6,7,8,9,10]
    N_t = 1000
        
    for k in N_c:
        
        num_of_heads = [i for i in range(0,k+1)]
        num_of_heads_nc.append(num_of_heads)
        
        result = macrostates(N_t,k)
        total,hc_array = macrostates_counter(result,k)
        freq_h = head_count(hc_array,k)
        probability = freq_h/N_t
        prob_nc.append(probability)
    
    print("Probability (at constant no. of trails):\n",prob_nc)
    print("\nProbability (at constant no. of coins)\n",prob_nt)

    fig2,ax2 = plt.subplots()
    
    for i in range(len(prob_nc)):
        ax2.plot(num_of_heads_nc[i],prob_nc[i],label = "No of coins = "+str(N_c[i]))
        ax2.scatter(num_of_heads_nc[i],prob_nc[i])
    ax2.set(xlabel = "No of Heads",ylabel = "Probabilty of macrostates",title = "Probility of macrostates V/s No. of Heads at constant no. of trails ("+str(N_t)+")")
    ax2.grid(ls = "--")
    ax2.legend()
    plt.show()
    
    N_c = 3
    N_t = 50
    no_of_trails  = []
    cum_prob_head = [] ; cum_prob_tail = [] ; cum_freq_head = [] ; cum_freq_tail = []
    prob_head = 0 ;  prob_tail = 0 ; freq_head = 0 ;  freq_tail = 0
    
    for k in range(1,N_t+1):
        no_of_trails.append(k)
        total_outcomes = 3*k
        result = macrostates(k,N_c)
        freq_head = freq_head + (sum(result[k-1]))
        freq_tail = freq_tail + (N_c - sum(result[k-1]))
        prob_head = freq_head/total_outcomes
        prob_tail = freq_tail/total_outcomes
        cum_freq_head.append(freq_head)
        cum_freq_tail.append(freq_tail)
        cum_prob_head.append(prob_head)
        cum_prob_tail.append(prob_tail)
        freq_head = cum_freq_head[k-1]
        freq_tail = cum_freq_tail[k-1]
        
    fig3,ax3 = plt.subplots()
    
    ax3.plot(no_of_trails,cum_prob_head,label = "Heads")
    ax3.plot(no_of_trails,cum_prob_tail,label = "Tails")
    ax3.scatter(no_of_trails,cum_prob_head)
    ax3.scatter(no_of_trails,cum_prob_tail)
    ax3.set(xlabel = "No. of trails",ylabel = "Cumulative probability of head and tails",title = "Cumulative probability Vs no of trails")
    ax3.grid(ls = "--")
    ax3.legend()
    plt.show()
        
