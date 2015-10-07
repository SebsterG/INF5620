import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
def solver(I,dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V):
  
  
  Nt = int(round(T/float(dt)))
 
  
  
  u = np.zeros((Nx+3,Ny+3))
  u_1 = np.zeros((Nx+3,Ny+3))
  u_2 = np.zeros((Nx+3,Ny+3))
  error_list = np.zeros(Nt+1)
  q = np.zeros((Nx+3,Ny+3))+qc  # mesh for q

  
  x = np.linspace(0,Lx,Nx+1)
  y = np.linspace(0,Ly,Ny+1)
  t = np.linspace(0, T, Nt+1)
  
  dx = float(x[1]-x[0])
  dy = float(y[1]-y[0])
  print dx ,dy
  
  
  
  #if dt > beta * dx*np.sqrt(qc):
    #sys.exit("dt is high af")
  #if dt > np.sqrt(1/qc)*np.sqrt((1/dx)**2+(1/dy)**2):
   # sys.exit("dt is high af")  
      
  
  K = dt*b/2.0
  C_x = float(dt**2/dx**2)
  C_y = float(dt**2/dy**2)
  
  xv = x[:,np.newaxis]
  yv = y[np.newaxis,:]
  
  #I = np.vectorize(I)
  Ix = range(1, u.shape[0]-1)
  Iy = range(1, u.shape[1]-1)
  It = range(0, t.shape[0])


  # Fetching the initial condition
  for j in Iy:
    for i in Ix:
      u_2[i,j]=I(x[i-1],y[j-1])

  
  u_2[-1,1:-1] = u_2[-3,1:-1]   # update ghost points
  u_2[1:-1,0] = u_2[1:-1,2]
  u_2[1:-1,-1] = u_2[1:-1,-3]
  u_2[0,1:-1] = u_2[2,1:-1]
  
  #plots 
  
  
  fig = plt.figure()
  plt.ioff()
  x,y = np.meshgrid(x,y)
  ax = Axes3D(fig)
  ax.plot_surface(x, y, u_2[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
  ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
  plt.savefig('waveplot_%04d.png' % (0))


  
  
  
  
  #first step
  #inner points
  u_1[1:-1,1:-1] = u_2[1:-1,1:-1]+ dt*V*(2 +dt*b) + \
    C_x*0.25*( (q[2:,1:-1]+q[1:-1,1:-1])*(u_2[2:,1:-1]-u_2[1:-1,1:-1]) -\
      (q[1:-1,1:-1]+q[:-2,1:-1])*(u_2[1:-1,1:-1]- u_2[:-2,1:-1])  )\
	+ C_y*0.25*( (q[1:-1,2:] + q[1:-1,1:-1])*(u_2[1:-1,2:]-u_2[1:-1,1:-1]) - \
	    (q[1:-1,1:-1]+q[1:-1,:-2])*(u_2[1:-1,1:-1]-u_2[1:-1,:-2])   )
  
  

  u_1[-1,1:-1] = u_1[-3,1:-1]   # update ghost points in first time step
  u_1[1:-1,0] = u_1[1:-1,2]
  u_1[1:-1,-1] = u_1[1:-1,-3]
  u_1[0,1:-1] = u_1[2,1:-1]
  
 
  
  fig = plt.figure()
  plt.ioff()
  ax = Axes3D(fig)
  ax.plot_surface(x, y, u_1[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')  
  ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
  plt.savefig('waveplot_%04d.png' % (1))
  #plt.show()
  
  
  
  # all time steps
  for n in range(Nt):
    # inner points
    
    u[1:-1,1:-1] = ((K-1)/(K+1))*u_2[1:-1,1:-1] + (2/(1+K))*u_1[1:-1,1:-1] + \
	    (0.5*C_x/(K+1))*((q[2:,1:-1]+q[1:-1,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
		    (q[1:-1,1:-1]+q[0:-2,1:-1])*(u_1[1:-1,1:-1]-u_1[0:-2,1:-1])) + \
			    (0.5*C_y/(K+1)) * ((q[1:-1:,2:]+q[1:-1,1:-1])*(u_1[1:-1,2:]-u_1[1:-1,1:-1]) - \
				    (q[1:-1:,1:-1]+q[1:-1,0:-2])*(u_1[1:-1,1:-1]-u_1[1:-1,0:-2])) 
					    
	
    
    u[-1,1:-1] = u[-3,1:-1]   # update ghost points in first time step
    u[1:-1,0] = u[1:-1,2]
    u[1:-1,-1] = u[1:-1,-3]
    u[0,1:-1] = u[2,1:-1]
    
    
    #plot and make figure
    
    num_frames = 100
    skip_frame = int(Nt/float(num_frames))
    if n % skip_frame == 0 or n == Nt-1:
      fig = plt.figure()
      plt.ioff()
      ax = Axes3D(fig)
      ax.plot_surface(x, y, u[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
      ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
      plt.savefig('waveplot_%04d.png' % (n+3))
    #plt.show()
      
    
    
    
    #update time steps
    #u_2[:,:] = u_1[:,:] ; u_1[:,:] = u[:,:]
    u_2, u_1, u = u_1, u, u_2  #faster 
    error = abs(u_exact(xv,yv,t[n])- u[1:-1,1:-1])
    error_list[n] = np.amax(error)
    
  error_mean = np.amax(error_list)
    
    
  return #error_mean 

def constant_solution():
  dt=0.01
  Lx = 10
  Ly = 10
  Nx = 50
  Ny = 50
  T=2
  mx = 2
  my = 2
  kx = mx*np.pi/Lx
  ky = my*np.pi/Ly
  w = np.sqrt(9.81*kx)
  A=1
  V= 0
  beta = 0.9
  b = 0
  qc = 1
  f = 0
  u_exact = lambda x, y, t: 2 
  error=solver(dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)







def plug_wave(direction):
  b = 1
  Lx = 1.0
  Ly = 1.0
  Nx = 150
  Ny = 150
  T = 4
  qc = 1
  u_exact = 0
  beta = 0.8
  f = 0
  V  = 0 
  #dt = (Lx/Nx)/(np.sqrt(qc))
  dt = np.sqrt(1/qc)*1/np.sqrt((1/(Lx/Nx))**2+(1/(Ly/Ny))**2)
  print "dt is:... ", dt
  if direction == "x":
    I = lambda x,y: 0 if abs(x-0.5*Lx) > 0.1 else 2 
  if direction == "y":
    I = lambda x,y: 0 if abs(y-0.5*Ly) > 0.1 else 2
  if direction == "middle":
    I = lambda x,y: 0 if abs(x-0.5*Lx) > 0.05 or abs(y-0.5*Ly) > 0.05 else 2

  solver(I,dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)

def standing_wave():
  
  Lx = 1.0
  Ly = 1.0
  Nx = 100.0
  Ny = 100.0
  T=2.0
  mx = 2
  my = 2
  kx = mx*np.pi/Lx
  ky = my*np.pi/Ly
  w = np.sqrt(9.81*kx)
  A=1
  V= 0
  beta = 0.9
  b = 0
  qc = 1
  f = 0
  dt = np.sqrt(1/qc)*1/np.sqrt((1/(Lx/Nx))**2+(1/(Ly/Ny))**2)
  u_exact = lambda x, y, t: A*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)*np.cos(w*t)  
  I = lambda x,y: A*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)
  print solver(I,dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)  


if __name__ == "__main__":
  #plug_wave(raw_input("x,y, or middle"))
  standing_wave()