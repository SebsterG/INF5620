from numpy import logspace, exp
from matplotlib.pyplot import *
import time
p = logspace(-6, -0.5, 101)

def Diff_forward(p):
	y = (1-exp(-p))/p
	return semilogx(p, y)
	

def Diff_backward(p):
	y = (exp(p)-1)/p
	return semilogx(p,y)

def Diff_centered(p):
	y = (exp(-p/2)-exp(p/2))/-p
	return semilogx(p,y),y 

Diff_forward(p)
Diff_forward(p)
Diff_backward(p)
Diff_centered(p)
show()

from sympy import *
s=Symbol("s")
y1 = (1-exp(-s))/s
y2 = (exp(s)-1)/s
y3 = (exp(-s/2)-exp(s/2))/-s


print y1.series(s, 0, 6)
print y2.series(s, 0, 6)
print y3.series(s, 0, 6)











