#!/usr/bin/env python
# coding: utf-8

# In[20]:


from sympy import solveset
from sympy.abc import x
from sympy import pprint
from sympy import re,im
import matplotlib.pyplot as plt
from matplotlib.axes import Axes



# In[31]:


# Problem 1-Part C) ROOT SOLVER for various Kc values
  
#Given values
Kc = [5,15,30,45,60,75]
Kp = 0.503 
tau = 12.4 #min

for cst in Kc:
    roots = solveset((tau*x+1)**3+(Kp**3)*cst,x)
    print("The roots for Kc = %d"%cst)
    for root in roots:
        print(root.round(4))
    print("")    


# In[22]:


#Problem 1-Part D) Plot root locus plot
fig, ax = plt.subplots()
ax.set_xlabel("real")
ax.set_ylabel("imaginary",rotation=0)
ax.spines["left"].set_position(("data", 0))
ax.spines["bottom"].set_position(("data", 0))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.set_label_coords(1,1.1)

for cst in Kc:
    roots = solveset((tau*x+1)**3+(Kp**3)*cst,x)
    for root in roots:
        #This is always true since all real numbers are complex numbers with 0j
        if root.is_complex:
            #cast root as complex since they are sympy objects
            val = complex(root)
            real = re(root)
            imaginary = im(root)
            plt.plot(val.real,val.imag,"bo")
            #plt.annotate(s,(val.real,val.imag),xytext=(0,cst-3),textcoords='offset points',arrowprops=dict(arrowstyle="->",
                            #connectionstyle="arc3"),)
        #Just in case the root is a different type   
        else:
            plt.plot(root,0,"bo")
            plt.plot(root,0,"bo")
            #plt.annotate(s,(val.real,val.imag), xytext=(-20, 20),
            #textcoords='offset points',ha="left", va='bottom',
            #arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
 


# In[23]:


#BODE STABILITY: Problem 2-PART C
from pylab import *

def G(w):
    Kp= 0.503  #process or disturbance gain here
    tau= 12.4  # time constant here
    Kc = 1
    Gol = Kp**3*Kc*(1-tau*1j*w)**3/(1+tau**2*w**2)**3
    return  Gol

f=logspace(-3,0)*2*pi # frequencies from 0.01 to 10, in angular frequencies
Gol=G(f)

AR=abs(Gol)
phi=angle(Gol)*180/pi


# Correct phi for postive angles
for i in range(len(phi)):
    if(phi[i]>0.0):
        phi[i]-=360
        
fig1, ax1 = plt.subplots(figsize=(10,5))
ax1.set_xlabel("Angular Frequency (rad/min)")
ax1.set_ylabel("Amplitude Ratio")
ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.plot(f,abs(Gol))

fig2, ax2 = plt.subplots(figsize=(10,5))
ax2.set_xlabel("Angular Frequency (rad/min)")
ax2.set_ylabel("Phase angle")
ax2.set_xscale('log')

ax2.set_yticks(np.arange(0,-300,-20))
ax2.plot(f,phi)


# In[24]:


#BODE STABILITY: Problem 2-Part F) Changing Kc until instability
from pylab import *


def G(w):
    Kp= 0.503  #process or disturbance gain here
    tau= 12.4  # time constant here
    Kc = 63
    Gol = Kp**3*Kc*(1-tau*1j*w)**3/(1+tau**2*w**2)**3
    return  Gol

f=logspace(-3,0)*2*pi # frequencies from 0.01 to 10, in angular frequencies
Gol=G(f)

AR=abs(Gol)
phi=angle(Gol)*180/pi


# Correct phi for postive angles
for i in range(len(phi)):
    if(phi[i]>0.0):
        phi[i]-=360
        
fig1, ax1 = plt.subplots(figsize=(10,5))
ax1.set_xlabel("Angular Frequency (rad/min)")
ax1.set_ylabel("Amplitude Ratio")
ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.plot(f,abs(Gol))

fig2, ax2 = plt.subplots(figsize=(10,5))
ax2.set_xlabel("Angular Frequency (rad/min)")
ax2.set_ylabel("Phase angle")
ax2.set_xscale('log')

ax2.set_yticks(np.arange(0,-300,-20))
ax2.plot(f,phi)

#Critical frequency is when phi (phase angle) is -180 degrees. Past -180 degrees, the system is instable for the process parameters
for angle,w,ar in zip(phi,f,AR):
    print("[phi,w,AR]=",angle,w,ar)


# In[25]:


#PROBLEM 4-8: Find AR for a single w
from pylab import *
def G(w):
    Kp= 1  #process or disturbance gain here
    tau= 7.9 # time constant here
    Kc = 2
    Gol = Kp*Kc/((2+tau*1j*w)*(1+tau*1j*w))
    return  Gol

freq = 5 #min
w = 2*pi/freq #rad/min
Gol = G(w)
AR=abs(Gol)
phi=angle(Gol)*180/pi
print("For w=",w,"rad/min")
print(AR,phi)


# In[26]:


#PROBLEM 4-8: C) Plot AR and phi vs. w for singular and double mixing tanks to determine which equipment should be used
from pylab import *

#Two tanks
def G(w):
    Kp= 1  #process or disturbance gain here
    tau= 7.9 # time constant here
    Kc = 1
    Gol = Kp*Kc/((tau*1j*w+1)**2)
    return  Gol

def singleTank(w):
    Kp= 1  #process or disturbance gain here
    tau =7.9 # time constant here
    Kc = 1
    Gol = 1/(tau*w*1+1)
    return  Gol

f=logspace(-3,0)*2*pi # frequencies from 0.01 to 10, in angular frequencies

Gol=G(f)
AR=abs(Gol)
phi=angle(Gol)*180/pi


# Correct phi for postive angles
for i in range(len(phi)):
    if(phi[i]>0.0):
        phi[i]-=360

Gol2=singleTank(f)
AR2=abs(Gol2)
phi2=angle(Gol2)*180/pi
for i in range(len(phi2)):
    if(phi2[i]>0.0):
        phi2[i]-=360        
        
fig1, ax1 = plt.subplots(figsize=(10,5))
ax1.set_xlabel("Angular Frequency (rad/min)")
ax1.set_ylabel("Amplitude Ratio")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title("Two Mixing Tanks in Series")
ax1.plot(f,abs(Gol),label="Two Tanks")
ax1.plot(f,abs(Gol2),label="Single Tank")
ax1.legend()
plt.show()
""""
fig2, ax2 = plt.subplots(figsize=(10,5))
ax2.set_xlabel("Angular Frequency (rad/min)")
ax2.set_ylabel("Phase angle")
ax2.set_xscale('log')

ax2.set_yticks(np.arange(0,-300,-20))
ax2.plot(f,phi)
"""

"""
fig3, ax3 = plt.subplots(figsize=(10,5))
ax3.set_xlabel("Angular Frequency (rad/min)")
ax3.set_ylabel("Amplitude Ratio")
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_title("Singular Mixing Tank")
ax3.plot(f,abs(Gol2))
"""

"""
fig4, ax4 = plt.subplots(figsize=(10,5))
ax4.set_xlabel("Angular Frequency (rad/min)")
ax4.set_ylabel("Phase angle")
ax4.set_xscale('log')

ax4.set_yticks(np.arange(0,-300,-20))
ax4.plot(f,phi2)
"""

