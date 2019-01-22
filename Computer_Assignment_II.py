from numpy import *
import numpy as np
import mpmath as mp
import sympy
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import *


t=0.4
b=0.96
d=0.1

def myFunction(z):
   kss = z[0]
   hss = z[1]
   css = z[2]

   F = empty((3))
   F[0] = (kss**t)*(hss**(1-t))-css-d*kss
   F[1] = t*(kss**(t-1))*(hss**(1-t))+(1-d)-(1/b)
   F[2] = css*(hss**(t))-(1-hss)*(1-t)*(kss**(t))
   return F

zGuess = array([2,0.5,1])
z = fsolve(myFunction,zGuess)
print(z)  #first value in the following order kss, hss, css

kss=z[0]
hss = z[1]
css = z[2]
print(kss)
print(hss)
print(css)

#%%
#Part c
runs = 0 # tracks how many times the loop runs
i = 0 #used later
cl = 0
ch = kss

k = []
h = []
c = []
k.append(kss/2)
h.append(0)
c.append(0)

k[0] = kss/2
done = False 

while(done == False):
    runs += 1
    c[0] = (cl + ch)/2  #setting the guess for initial consumption    
    for i in range(5000):
        if i == 0:
            def hls(z): #calculating h0 given c0 and k0
                F = (1-z)*(1-t)*(k[0]**(t))-c[0]*(z**(t))
                return F            
            zGuess = array([0.5])
            zsolution = fsolve(hls,zGuess)
            print(zsolution)
            h.append (zsolution)
            
        else:
            k.append(k[i-1]-d*k[i-1]+(k[i-1]**t)*(h[i-1]**(1-t))-c[i-1])  # finding this periods k
            h.append(((1-t)*k[i])/((1-t)*k[i] + b*t*c[i-1])) #finding this periods h
            c.append(c[i-1]*b*((t*(k[i]**(t-1))*(h[i]**(1-t)))+1-d)) #finding this periods c
            
    if np.absolute(c[i] - css) < 10**(-4): #checking if the consumption gets close to steady state
        ksolution = k
        hsolution = h
        csolution = c #saving solutions and ending the while loop
        done = True
        break

    if max(k) < kss:
        ch = c[0]  
        k = []
        h = []
        c = []
        k.append(kss/2)
        h.append(0)
        c.append(0)
        break
        
    if max(k) > kss: 
        cl = c[0] 
        k = []
        h = []
        c = []
        k.append(kss/2)
        h.append(0)
        c.append(0)
        break
#Part d
plt.plot(ksolution)
plt.xlabel("loops")
plt.ylabel("k value")
plt.title("Part c")
plt.show()

plt.plot(hsolution)
plt.xlabel("loops")
plt.ylabel("h value")
plt.title("Part c")
plt.show()

plt.plot(csolution)
plt.xlabel("loops")
plt.ylabel("c value")
plt.title("Part c")
plt.show()

#%%
## Part e
runs = 0 # tracks how many times the loop runs
i = 0 #used later
cl = 0
ch = kss

k = []
h = []
c = []
k.append(2*kss)
h.append(0)
c.append(0)

k[0] = 2*kss
done = False 

while(done == False):
    runs += 1
    c[0] = (cl + ch)/2  #setting the guess for initial consumption    
    for i in range(5000):
        if i == 0:
            def hls(z): #calculating h0 given c0 and k0
                F = (1-z)*(1-t)*(k[0]**(t))-c[0]*(z**(t))
                return F            
            zGuess = array([0.5])
            zsolution = fsolve(hls,zGuess)
            print(zsolution)
            h.append (zsolution)
            
        else:
            k.append(k[i-1]-d*k[i-1]+(k[i-1]**t)*(h[i-1]**(1-t))-c[i-1])  # finding this periods k
            h.append(((1-t)*k[i])/((1-t)*k[i] + b*t*c[i-1])) #finding this periods h
            c.append(c[i-1]*b*((t*(k[i]**(t-1))*(h[i]**(1-t)))+1-d)) #finding this periods c
            
    if np.absolute(c[i] - css) < 10**(-4): #checking if the consumption gets close to steady state
        ksolution = k
        hsolution = h
        csolution = c #saving solutions and ending the while loop
        done = True
        break

    if max(k) < kss:
        ch = c[0]  
        k = []
        h = []
        c = []
        k.append(kss*2)
        h.append(0)
        c.append(0)
        break
        
    if max(k) > kss: 
        cl = c[0] 
        k = []
        h = []
        c = []
        k.append(kss*2)
        h.append(0)
        c.append(0)
        break
#Part d
plt.plot(ksolution)
plt.xlabel("loops")
plt.ylabel("k value")
plt.title("Part e")
plt.show()

plt.plot(hsolution)
plt.xlabel("loops")
plt.ylabel("h value")
plt.title("Part e")
plt.show()

plt.plot(csolution)
plt.xlabel("loops")
plt.ylabel("c value")
plt.title("Part e")
plt.show()
