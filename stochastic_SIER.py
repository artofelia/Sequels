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
    ne = 3 #number of latent compartments
    ni = 10 #number of inf compartments
    t = 0.0 #time
    ent = 300.0; #end time
    inter_size = int(4.0) #size of discrete time interval
    bins = int(ent/inter_size)
    
    
    ipop = 5000 #inital population
    iinf = 30 #inital infection
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
            f_ppl[inf_person] = f_ppl[inf_person]+ [g]
            
        elif event <= lam_s+lam_e:
            #move person in one of latent compartments
            comp = 1+choice( [i for i, j in enumerate(map(bool, ppl[1:ne+1] )) if j] ) #pick nonempty latent compartment

            #print 'event is LATENT', ppl[comp]
            chosen = choice( ppl[comp] )
            ppl[comp].remove(chosen)
            ppl[comp+1].append(chosen)
        elif event <= lam_s+lam_e+lam_i:
            #move person in one of infected compartments
            comp = 1+ne+choice( [i for i, j in enumerate(map(bool, ppl[1+ne:1+ne+ni] )) if j] )
            #print 'event is INF',  ppl[1+ne:1+ne+ni]

            chosen = choice( ppl[comp] )
            ppl[comp].remove(chosen)
            if comp!=ne+ni:
                ppl[comp+1].append(chosen)
        else:
            print 'something bad'
        t = t+tao
        #print 'end step', t
        #print 'ppl status', ppl
   
   
    bdist = [[]]*(bins+1) #backward distributions
    fdist = [[]]*(bins+1) #forward distributions
    bmean = [0]*(bins+1) #backeard mean
    fmean = [0]*(bins+1) #forward mean
    
    for i in range(ipop): # i stands for person
        if inf_time[i]==-1: #didnt get infected during epidemic
            continue
        tdic = int( (inf_time[i]-(inf_time[i]%inter_size)) / inter_size ) #sort by discrete time of infection
        #if b_ppl[i]:
        #    gdic = int( (b_ppl[i]-(b_ppl[i]%inter_size)) / inter_size ) #dicrete size of gen interval
        #   bdist[tdic].append(gdic)
        if f_ppl[i]:
            for e in f_ppl[i]:
                gdic = int( (e-(e%inter_size)) / inter_size )
                fdist[tdic] = fdist[tdic] + [gdic]
    
    #extract means
    for i in range(0,len(bdist)):
        #bmean[i] = mean(bdist[i])
        if fdist[i]:
            fmean[i] = mean(fdist[i])
        
    #print bmean, fmean
    #print fdist
    #plt.plot(fmean)
    #plt.show()
    return array(fmean)
    
def alg():
    reps = 10
    mn = model()
    for i in range(1,reps):
        mn += model()
    mn /= reps
    print mn
    plt.plot(mn)
    plt.show()
      
alg()
    
    
