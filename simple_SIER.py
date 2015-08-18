
from numpy import *
import matplotlib.pyplot as plt
plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt


#Runge-Kutta method
def kutta(x0, f, steps, dt): #initail vals, diff eq, number of steps, time step
    x = x0 #s,e1,e2,i1,in
    arrx = [];
    for t in range(0,steps):
        k1 = f(x)*dt
        k2 = f(x+.5*k1)*dt
        k3 = f(x+.5*k2)*dt
        k4 = f(x+k3)*dt
        x = x + (1/6.0)*(k1+2*k2+2*k3+k4)
        arrx.append(x)
    return arrx
    


def model():
    #paramets
    bet = .5
    sig = .5
    gam = .5
    ne = 3
    ni = 3
    
    def f(p): #differential equation
        s = p[0] #susptible
        lat = p[1:1+ne] #latents
        inf = p[ne+1:1+ne+ni] #infected
        np = []
        #print 'all', p
        #print 'lat', lat
        #print 'inf', inf
         
        np.append( -bet*s*inf.sum() ) #d suseptible
        
        np.append( bet*s*inf.sum()-sig*lat[0] ) #d latent1 
        for i in range(1, ne):
            np.append( sig*(lat[i-1]-lat[i]) ) #d latent i+1
         
        np.append( sig*lat[ne-1]-sig*inf[0] ) #d infected 1
        for i in range(1, ni):
            np.append( gam*(inf[i-1]-inf[i]) ) #d infected i+1
    
        return array(np)
        
    x0 = zeros((1+ne+ni, 1), dtype=float64) #initial vals
    x0[0] = 1000
    x0[ne+1] = 1
    #print 'init', x0
    
    steps = 600
    arrx = kutta(x0, f, steps, .01)
    
    sus = [e[0] for e in arrx]
    lat1 = [e[1] for e in arrx]
    inf1 = [e[ne+1] for e in arrx]
    t = range(0,steps)
    
    plt.plot(t, sus, 'r--', t, lat1, 'g--', t, inf1)
    plt.ylabel('Suseptible')
    plt.xlabel('Time')
    plt.show()
      
model()
    
    
