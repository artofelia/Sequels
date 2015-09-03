#ARE PEOPLE MOLECULES???

from numpy import *
from random import *
import pickle
import matplotlib.pyplot as plt
plt.rcdefaults()

def model():
    #parameters
    bet = .5 #contact rate
    sig = .5 #transfer rate
    gam = .5 #recovery rate
    ne = 3 #number of latent compartments
    ni = 3 #number of inf compartments
    t = 0.0 #time
    ent = 100; #end time
    inter_size = int(1.0) #size of discrete time interval
    bins = int(ent/inter_size) #number of discrete time intervals
    
    ipop = 10000 #inital population
    iinf = 100 #inital infection
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
    intervals_info = [] # all generation intervals, (time of infection of y, person y who got infected, person x who infected, time of x infection, generation interval)
    npp = [] #record number of people in each compartment per time (testing purposes)
   
    g_mean_time = [] #record evolving generating mean in list
    g_mean_current = 0;  #current evolving mean
    g = 0; #current generation time at this step
    intrinsic_gmean = ne*1/sig + ni*1/gam #average of intrensic interval
    
    steps = 0; #record number of steps
    
    while t < ent: #stochastic step
        steps += 1
        
        s = len(ppl[0]) #numer of susptible
        lat = array( [len(ppl[i]) for i in range(1,1+ne) ] ) #number of latents
        inf = array( [len(ppl[i]) for i in range(1+ne,1+ne+ni) ] ) #number of infected
        npp.append([s,lat[0],inf[0], t]) # keep data of few compartments for testing purposes
        
        if inf.sum()==0:
            print 'EPIDEMIC ENDED at', t
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
            
            intervals_info.append([t+tao, sus_person, inf_person, inf_time[inf_person], g])
            
        elif event <= lam_s+lam_e:
            #move person in one of latent compartments
            comp = 1+choice( [i for i, j in enumerate(map(bool, ppl[1:ne+1] )) if j] ) #pick nonempty latent compartment
            
            #print 'event is LATENT', ppl[comp]
            chosen = choice( ppl[comp] ) #choose person in compartment
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
        
        #update evolving generation mean
        g_mean_current +=g 
        g_mean_time.append([g_mean_current/steps, intrinsic_gmean])
   
    sus_data = [i[0] for i in npp]
    lat_data = [i[1] for i in npp]
    inf_data = [i[2] for i in npp]
    time_data = [i[3] for i in npp]
    plt.figure(1)
    susplot, = plt.plot(time_data, sus_data, 'r--')
    infplot, = plt.plot(time_data, inf_data, 'g--')
    latplot, = plt.plot(time_data, lat_data, 'b--')
    plt.legend([susplot, infplot, latplot], ["suseptible", "infected first cmpt", "latent first cmpt"])
    plt.title('Suseptibles, Latents, Infected evolutions')
    plt.show()

    gmean,   = plt.plot(time_data, [i[0] for i in g_mean_time], 'r')
    intrinsic_gmean,   = plt.plot(time_data, [i[1] for i in g_mean_time], 'b')
    plt.title('evolution of generation mean')
    plt.show()
    
    #calculate intervals
    bdist = [[]]*(bins+1) #backward distributions
    fdist = [[]]*(bins+1) #forward distributions
    bmean = [-1]*(bins+1) #backeard mean
    fmean = [-1]*(bins+1) #forward mean
    
    for i in range(ipop): # i stands for person
        if inf_time[i]==-1: #didnt get infected during epidemic
            continue
        tdic = int( (inf_time[i]-(inf_time[i]%inter_size)) / inter_size ) #sort by discrete time of infection
        if b_ppl[i]:
            gdic = int( (b_ppl[i]-(b_ppl[i]%inter_size)) / inter_size ) #dicrete size of gen interval
            bdist[tdic] = bdist[tdic] + [gdic] #add discrete gen interval to backwards distribution at discrete time
        if f_ppl[i]:
            for e in f_ppl[i]:
                gdic = int( (e-(e%inter_size)) / inter_size )  #dicrete size of gen interval
                fdist[tdic] = fdist[tdic] + [gdic]
                
    #save to file            
    #with open('backward_intervals_2.txt', 'wb') as f:
    #    pickle.dump(intervals_info, f)
   
    #extract means
    for i in range(0,len(bdist)):
        if bmean[i]:
            bmean[i] = mean(bdist[i])
        if fdist[i]:
            fmean[i] = mean(fdist[i])
        

   
    
    return [array(fmean), array(bmean)]
    
def alg():
    reps = 1 #number of times to repeat algorithm
    mn = model() #run simlutaion first time
    fmn = mn[0] #forward means
    bmn = mn[1] #backward means
    bins = range(len(fmn)) #number of time intervals
    for i in range(1,reps):
        mn = model() #run simulation again
        fmn += mn[0]
        bmn += mn[1]
    fmn /= reps
    bmn /= reps
    
    plt.figure(2)
    fmplot, = plt.plot(bins, fmn, 'g--')
    bmplot, = plt.plot(bins, bmn, 'b--')
    plt.title('forward and backward means')
    plt.legend([fmplot, bmplot], ["forward mean", "backward mean"])
    plt.show()
      
alg()
    
    
