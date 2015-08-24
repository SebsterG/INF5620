from matplotlib.pylab import *


def mesh_function(f,t):
  values= zeros(len(t))
  values = f(t) 
  return values, t
    
  
def function1(t):
  return exp(-t)
def function2(t):
  return exp(-3*t)
n1 = 3.
dt =0.1
t1=linspace(0,n1,(n1/dt)+1)
n2=4.
t2 = linspace(3,n2,(n2/dt)+1)


v1,time1 = mesh_function(function1, t1)
v2, time2 = mesh_function(function2,t2)

figure(1)
plot(time1,v1)
plot(time2,v2)
show()