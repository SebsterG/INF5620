from matplotlib.pylab import *
import numpy as np
def differentiate(u, dt):
  n=len(u)
  d = np.zeros(n)
  for i in range (1,n-1):
    d[i] = (u[i+1]-u[i-1])/(2*dt) 
  d[0] = (u[1]-u[0])/dt 
  d[-1] = (u[-1]-u[-2])/dt
  return d
tid = linspace(0,4,4/0.1+1)
def func(x):
  return x**2
def mesh_function(f,t):
  values= np.zeros(len(t))
  values = f(t) 
  return values

def realfunc(x):
  return 2*x

figure()
title('Difference between numerical and analytical function x^2')
plot(tid,differentiate(func(tid),0.1), '*r')
plot(tid,realfunc(tid))
xlabel('time')
ylabel('x^2')

show()