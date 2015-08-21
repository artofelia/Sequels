#ARE PEOPLE MOLECULES???

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
    ent = 50; #end time
    inter_size = int(1.0) #size of discrete time interval
    inter = range(0, ent, inter_size) #discrete intervals
    
    
    ipop = 1000 #inital population
    iinf = 3 #inital infection
    inf_time = array([-1]*ipop)#records time of infection for each person
    inf_time[0:iinf] = 0
    
    sus_ppl = range(iinf,ipop) # index of suseptiple ppl
    inf_ppl = range(iinf) # index of infected ppl
    
    ppl = [] #holds states of all ppl
    ppl.append(sus_ppl)
    for i in range(0, ne):
       ppl.append([])
     
    ppl.append(inf_ppl)
    for i in range(1, ni):
        ppl.append([])
        
    #print 'ppl', ppl
    
    b_ppl = [0]*(iinf)+[None]*(ipop-iinf) #backwards interval for ppl
    f_ppl = [[]]*ipop #forwards interval for ppl
    
    b_dist = []
    
    while t < ent: #stochastic step
        s = len(ppl[0]) #numer of susptible
        lat = array( [len(ppl[i]) for i in range(1,1+ne) ] ) #number of latents
        inf = array( [len(ppl[i]) for i in range(1+ne,1+ne+ni) ] ) #number of infected
        if inf.sum()==0:
            print 'EPIDEMIC ENDED'
            break
        
        ''' how does the actual s[t], i[t] update with the lamda changes
        since s[t], i[t] gets updaated per individual person '''
        
        #lamdbas = respective chance that individaul in compartmanet a changes to b
        lam_s = bet*s*inf.sum()/ipop
        lam_e = sig/ne * lat.sum()
        lam_i = gam/ni * inf.sum()
        lam = lam_s+lam_e+lam_i
        #print 'LAM', lam
        tao = expovariate(lam) #next time step
        #why is mean,the number of unrecovered individuals?
        event = uniform(0, lam) #event type
        
        #print 'lams,lame,lami,lam,tao,event', lam_s,lam_e,lam_i,lam,tao,event
        
        if event < lam_s:
            #print 'event is SUS'
            sus_person = int(choice(ppl[0])) #pick suseptible person
            ppl[0].remove(sus_person)
            ppl[1].append(sus_person)
            inft = t+tao #infection time
            inf_time[sus_person] = inft
            
            #sus_person is infected
            inf_person = int(choice( concatenate( [ppl[i] for i in range(1+ne,1+ne+ni) ] ) )) #pick infected person
            g = t+tao-inf_time[inf_person]
            
            b_ppl[sus_person] = g
            f_ppl[inf_person].append(g)
            
        elif event <= lam_s+lam_e:
            
            #move person in one of latent compartments
            comp = 1+choice( [i for i, j in enumerate(map(bool, ppl[1:ne+1] )) if j] ) #pick nonempty latent compartment

            #print 'event is LATENT', ppl[comp]
            chosen = choice( ppl[comp] )
            ppl[comp].remove(chosen)
            ppl[comp+1].append(chosen)
        else:
            
            #move person in one of infected compartments
            comp = 1+ne+choice( [i for i, j in enumerate(map(bool, ppl[1+ne:1+ne+ni] )) if j] )
            #print 'event is INF',  ppl[1+ne:1+ne+ni]

            chosen = choice( ppl[comp] )
            ppl[comp].remove(chosen)
            if comp!=ne+ni:
                ppl[comp+1].append(chosen)
        t = t+tao
        #print 'end step', t
        #print 'ppl status', ppl
    
    batt = 5 #look at bdist at discrete time batt
    bdistt = [] #discrete backwards dist at time batt
    for i in range(ipop):
        if not b_ppl[i]:
            continue
        tdic = int( (inf_time[i]-(inf_time[i]%inter_size)) / inter_size ) #discrete time of infection
        if tdic==batt-1 : #pick out individuals who got sick at batt
            gdic = int( (b_ppl[i]-(b_ppl[i]%inter_size)) / inter_size ) #dicrete size of gen interval
            bdistt.append(gdic)
            
    fdistt = [] #discrete forward dist at time batt
    for i in range(ipop):
        if not f_ppl[i]:
            continue
        tdic = int( (inf_time[i]-(inf_time[i]%inter_size)) / inter_size ) #discrete time of infection
        if tdic==batt-1 :
            for e in f_ppl[i]:
                gdic = int( (e-(e%inter_size)) / inter_size )
                fdistt.append(gdic) 
        
    plt.hist(fdistt,bins=50)

    plt.show()
      
model()
    
    
