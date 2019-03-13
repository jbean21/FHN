# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 09:29:16 2018

@author: Martin Shahi
"""

#RK4 of HR eqns

#imports
import matplotlib.pyplot as plt


#define variables
a = 1
b = 3
c = 1
d = 5
I = 1
r = 0.001
s = 4
x_R = -8.0/5

X = [0,10,1] #will hold the current value of each variable
t = 0 #holds current time
t_end = 5

#define step size

h = 0.01 #for example

#initialise arrays to store results
x_values = []
y_values = []
z_values = []
time = []

#define equations and methods

def eqnx(x):
    return[x[1]+b*((x[0])**2)-a*((x[0])**3)-x[2]+I, x[1], x[2]]

def eqny(x):
    return[x[0], c-d*((x[0])**2)-x[1], x[2]]
    
def eqnz(x):
    return[x[0], x[1], (r*(s*(x[0]-x_R)-x[2]))]
    
def k1(x,f):
    value = f(x)
    return [h*i for i in value]

def k2(x,f):
    epsilon = [i/2 for i in k1(x,f)]
    value = f(x+epsilon)
    return [h*i for i in value]
    
def k3(x,f):
    epsilon = [i/2 for i in k2(x,f)]
    value = f(x+epsilon)
    return [h*i for i in value]
    
def k4(x,f):
    value = f(x+k3(x,f))
    return [h*i for i in value]
    
def RK4(x,f):
    return x + [i*(1/6) for i in (k1(x,f)+2*k2(x,f)+2*k3(x,f)+k4(x,f))]
    
#run

while (t<t_end):
    x_values.append(X[0])
    y_values.append(X[1])
    z_values.append(X[2])
    time.append(t)
    X = RK4(X,eqnx)
    X = RK4(X,eqny)
    X = RK4(X,eqnz)
    t += h
    print(len(X))
    

#plot results

plt.plot(time,x_values)
plt.show()










































