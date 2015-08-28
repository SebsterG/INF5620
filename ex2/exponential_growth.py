from numpy import *
from matplotlib.pyplot import *

def exact(p):
	return exp(p)
def theta_calc(theta,p):
	return (1+(1-theta)*p)/(1+theta*p)


p = linspace(1,15,100)
"""plot(p,exact(p))"""
plot(p,theta_calc(0.5,p))
show()
"""
theta_values = [0,0.5,1]
for i in theta_values :
	plot(p,theta_calc(i,p))
	
show()
	"""