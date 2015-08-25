from numpy import logspace, exp
from matplotlib.pyplot import *
p = logspace(-6, -0.5, 101)
y = (1-exp(-p))/p
semilogx(p, y)


def Diff_forward(u,dt):
	n = len(u)
	d_forward = zeros(n)
	for i in range(n-1):
		d_forward[n] = (u[n+1] - u[n])/dt
	return d_forward

def Diff_backward(u,dt):
	n =len(u)
	d_backward = zeros(n)
	for i in range(n-1):
		d_backward[u] = (u[n] - u[n-1])/dt
	return d_backward

def Diff_centered(u,dt):
	n=len(u)
	d_centered = zeros(n)
	for i in range(n-1):
		d_centered[u] = (u[n+0.5] - u[n-0.5])/dt
	return d_centered

