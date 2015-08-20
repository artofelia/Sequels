
from numpy import *
from random import *
import matplotlib.pyplot as plt
plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt


def model():
    #parameters
    bet = .5 #contact rate
    sig = .5 #transfer rate
    gam = .5 #rec rate
    ne = 2 #number of latent compartments
    ni = 2 #number of inf compartments
    t = 0.0 #time
    interval = 10 #size of discrete time interval
    steps = 600 #numer of iterations per run
    
    ipop = 10 #inital population
    iinf = 2 #inital infection
    inf_time = ones( (ipop,1) )*-1 #records time of infection for each person
    inf_time[0:iinf] = 0
    
    sus_ppl = arange(iinf,ipop) # index of suseptiple ppl
    inf_ppl = arange(iinf) # index of infected ppl
    
    ppl = [] #holds states of all ppl
    ppl.append(sus_ppl)
    for i in range(0, ne):
       ppl.append(array([]))
     
    ppl.append(inf_ppl)
    for i in range(1, ni):
        ppl.append(array([]))
        
    print 'ppl', ppl
    
    def f(): #stochastic step
        s = len(ppl[0]) #numer of susptible
        lat = array( [len(ppl[i]) for i in range(1,1+ne) ] ) #number of latents
        inf = array( [len(ppl[i]) for i in range(1+ne,1+ne+ni) ] ) #number of infected
        ''' how does the actual s[t], i[t] update with the lamda changes
        since s[t], i[t] gets updaated per individual person '''
        
        #what do the lambdas stand for?
        lam_s = bet*s*inf.sum()/ipop
        lam_e = sig/ne * lat.sum()
        lam_i = gam/ni * inf.sum()
        lam = lam_s+lam_e+lam_i
        tao = expovariate(lam) #next time step
        #why is mean,the number of unrecovered individuals?
        event = uniform(0, lam) #event type
        
        if event < lam_s:
            sus_person = choice(ppl[0]) #pick suseptible person
            inft = t+tao #infection time
            inf_time[sus_person] = inft
            
            #sus_person is infected
            inf_person = choice( concatenate( [ppl[i] for i in range(1+ne,1+ne+ni) ] ) ) #pick infected person
            g = t+tao-inf_time[inf_person]
        else:
            change_person = choice( concatenate( [ppl[i] for i in range(1+ne,1+ne+ni) ] ) ) #pick person to move
            
            
        
        t = t+tao
        return np
        
    
    f()
    #arrx = kutta(x0, f, steps, .01)
    
    #sus = [e[0] for e in arrx]
    #lat1 = [e[1] for e in arrx]
    #inf1 = [e[ne+1] for e in arrx]
    #t = range(0,steps)
    
    #plt.plot(t, sus, 'r--', t, lat1, 'g--', t, inf1)
    #plt.ylabel('Suseptible')
    #plt.xlabel('Time')
    #plt.show()
      
model()
    
    
