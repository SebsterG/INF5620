import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as math



def Solver(I, V, f_inn, q_inn, b, Lx, Ly, Nx, Ny, dt, T, u_real, user_Action):

	# Making space
	Nt = int(T/dt)
	x = np.linspace(0,Lx,Nx+1)
	y = np.linspace(0,Ly,Ny+1)
	t = np.linspace(0, T, Nt+1) 
	dx = x[1]-x[0]
	dy = y[1]-y[0]

	xv = x[:,np.newaxis]
	yv = y[np.newaxis,:]

	u = np.zeros((Nx+3,Ny+3)) 
	u_1 = np.zeros((Nx+3,Ny+3))
	u_2 = np.zeros((Nx+3,Ny+3))

	Cx2 = (dt/dx)**2
	Cy2 = (dt/dy)**2
	K = (dt*b/2)	

	c = 1
	f = np.zeros((Nx+3,Ny+3))
	q = np.zeros((Nx+3,Ny+3))+c

	Ix = range(1, u.shape[0]-1)
	Iy = range(1, u.shape[1]-1)
	It = range(0, t.shape[0])

	stability_limit = (1/c)*(np.sqrt((1/dx**2)+(1/dx**2)))**-1
	if dt >= stability_limit:
		print "unstable"

	# Fetching the initial condition
	for j in Iy:
		for i in Ix:
			u_1[i,j]=I(x[i-1],y[j-1])

	if user_Action == True:
		n = 0
		user_action(u_1,x,y,n)
	#user_Action = False
	
	step1 = True
	u = advance_vec(u, u_1, u_2, f, q, Cx2, Cy2, dt, V, b, Nx, Ny, step1)
	
	if user_Action == True:
		n = 1 
		user_action(u,x,y,n)

	step1 = False
	u_2, u_1, u = u_1, u, u_2

	error_count = True
	error_list = np.zeros(Nt)
	counter = 0
	counter2 = 1
	for n in range(1,Nt):
		u = advance_vec(u, u_1, u_2, f, q, Cx2, Cy2, dt, V, b, Nx, Ny, step1)
		u_2, u_1, u = u_1, u, u_2

		counter = counter + 1
		d = counter % 10
		if d == 1:
			counter2 += 1
			if user_Action == True:
				user_action(u,x,y,counter2)
				user
		if error_count == True:
			u_e = u_real(xv,yv,t[n])
			error = abs(u_e-u[1:-1,1:-1])
			error = np.amax(error)
			error_list[n] = error
	if error_count == True:
		max_error = max(error_list)
		return max_error 

	
def advance_vec(u, u_1, u_2, f, q, Cx2, Cy2, dt, V, b, Nx, Ny,step1):
	dt2 = dt**2
	K = (dt*b/2)
	# Hardcoding for simplicity (temp)
	V = 0	
	if step1:
		
		u_1[0,1:-1] = u_1[2,1:-1]
		u_1[Nx+2,1:-1] = u_1[Nx,1:-1]

		u_1[1:-1,0] = u_1[1:-1,2]
		u_1[1:-1,Ny+2] = u_1[1:-1,Ny] 


		# All interior space at first time step
		u[1:-1,1:-1] = (1+K)*dt*V + u_1[1:-1,1:-1] + \
			(Cx2/4)*((q[2:,1:-1]+q[1:-1,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
				(q[1:-1,1:-1]+q[0:-2,1:-1])*(u_1[1:-1,1:-1]-u_1[0:-2,1:-1])) + \
					(Cy2/4)*((q[1:-1:,2:]+q[1:-1,1:-1])*(u_1[1:-1,2:]-u_1[1:-1,1:-1]) - \
						(q[1:-1:,1:-1]+q[1:-1,0:-2])*(u_1[1:-1,1:-1]-u_1[1:-1,0:-2])) + \
							(dt2/2)*f[1:-1,1:-1]

			
	else:

		u_1[0,1:-1] = u_1[2,1:-1]
		u_1[Nx+2,1:-1] = u_1[Nx,1:-1]

		u_1[1:-1,0] = u_1[1:-1,2]
		u_1[1:-1,Ny+2] = u_1[1:-1,Ny] 
	

		u[1:-1,1:-1] = ((K-1)/(K+1))*u_2[1:-1,1:-1] + (2/(1+K))*u_1[1:-1,1:-1] + \
			(0.5*Cx2/(K+1))*((q[2:,1:-1]+q[1:-1,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
				(q[1:-1,1:-1]+q[0:-2,1:-1])*(u_1[1:-1,1:-1]-u_1[0:-2,1:-1])) + \
					(0.5*Cy2/(K+1)) * ((q[1:-1:,2:]+q[1:-1,1:-1])*(u_1[1:-1,2:]-u_1[1:-1,1:-1]) - \
						(q[1:-1:,1:-1]+q[1:-1,0:-2])*(u_1[1:-1,1:-1]-u_1[1:-1,0:-2])) + \
							(dt2/(K+1))*f[1:-1,1:-1]

	
	return u


def user_action(u,x,y,n): 

		Path = "/home/fenics/Desktop/Inf5620/Wave_Project/Plots"
		plt.ioff()
		x,y = np.meshgrid(x,y)
		fig = plt.figure(n)
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(x, y, u[1:-1,1:-1])
		ax.set_zlim(-2.0, 2.0)
		plt.savefig('%s/waveplot_%04d.png' % (Path,n))


def pulse(direction):

	b = 1
	V = lambda x,y: 0
	q_inn = lambda x,y: 1
	f_inn = lambda x,y,t: 0
	dt = 0.001
	T = 1
	Lx = 1
	Ly = 1

	if direction == 'y':
		Nx = 100
		Ny = 100
		def I(x,y):
			sigma = 0.025
			if math.fabs(x-0.5*Lx) > sigma :
				out = 0
			else:
				out = 1
			return out

	if direction == 'x':
		Nx = 100
		Ny = 100
		def I(x,y):
			sigma = 0.025
			if math.fabs(y-0.5*Ly) > sigma :
				out = 0
			else:
				out = 1
			return out

	if direction == 'corner':
		Nx = 100
		Ny = 100
		def I(x,y):
			sigma = 0.05
			if abs(x-0.2*Lx) > sigma or abs(y-0.6*Ly) > sigma:
				out = 0
			else:
				out = 2
			return out


	Solver(I, V, f_inn, q_inn, b, Lx, Ly, Nx, Ny, dt, T,u_real, user_Action=False)


def Standing_conv():

	mx = 1
	my = 1
	b = 0
	V = lambda x,y: 0
	q_inn = lambda x,y: 1
	f_inn = lambda x,y,t: 0

	dt = 0.01
	T = 0.5
	Nt = int(T/dt)
	Lx = 1
	Ly = 1
	Nx = 10
	Ny = 10
	t = np.linspace(0, T, Nt+1) 
	x_t = np.linspace(0,Lx,Nx+1)
	y_t = np.linspace(0,Ly,Ny+1)
	xv_t = x_t[:,np.newaxis]
	yv_t = y_t[np.newaxis,:]

	def I(x,y):
		A = 1
		return A*np.cos(mx*np.pi*x/Lx)*np.cos(mx*np.pi*y/Ly)

	def u_real(x,y,t):
		w = np.pi
		A = 1
		return A*np.cos(mx*np.pi*x/Lx)*np.cos(mx*np.pi*y/Ly)*np.cos(np.sqrt(2)*w*t)
	Fx = 10	
	Fy = 10
	err = []
	h_list = [0.1,0.01,0.001,0.0005]
	for h in h_list:
		dt = 1*h
		Nx = int(round(Lx/(Fx*h)))
		Ny = int(round(Ly/(Fy*h)))

		error = Solver(I, V, f_inn, q_inn, b, Lx, Ly, Nx, Ny, dt, T, u_real, user_Action=False)
		err.append(error)
	print "The error is:" 
	print err

	r = [np.log(err[i-1]/err[i])/ \
		np.log(h_list[i-1]/h_list[i]) \
			for i in range(1, len(h_list), 1)]
	print "... with a convergence rate of:" 
	print r

Standing_conv()