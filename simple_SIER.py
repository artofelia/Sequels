
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
    
    def f(p): #differential equation
        s = p[0]
        l1 = p[1]
        l2 = p[2]
        i1 = p[3]
        i2 = p[4]
         
        f1 = -bet*s*(i1+i2) #ds
        f2 = bet*s*(i1+i2)-sig*l1 #dl1
        f3 = sig*(l1-l2) #dl2
        f4 = sig*l2-sig*i1 #di1
        f5 = gam*(i1-i2) #di2
        
        return array([f1,f2,f3,f4,f5])
        
    x0 = array([100,0,0,1,0], dtype=float64) #initial vals
    steps = 50
    arrx = kutta(x0, f, steps, .1)
    s1 = [e[0] for e in arrx]
    s2 = [e[1] for e in arrx]
    print s1
    t = range(0,steps)
    plt.plot(t, s1, 'r--', t, s2, 'g--')
    plt.ylabel('Suseptible')
    plt.xlabel('Time')
    plt.show()
      
model()
    
    
