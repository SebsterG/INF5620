import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
def solver(dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V):
  
  
  Nt = int(round(T/float(dt)))
 
  
  
  u = np.zeros((Nx+3,Ny+3))
  u_1 = np.zeros((Nx+3,Ny+3))
  u_2 = np.zeros((Nx+3,Ny+3))
  error_list = np.zeros(Nt+1)
  q = np.zeros((Nx+3,Ny+3))+qc  # mesh for q

  
  x = np.linspace(0,Lx,Nx+1)
  y = np.linspace(0,Ly,Ny+1)
  t = np.linspace(0, T, Nt+1)
  
  dx = x[1]-x[0]
  dy = y[1]-y[0]
  
  
  #if dt > (1/qc)*(np.sqrt((1/dx**2)+(1/dx**2)))**-1: # Stability quits if not reached
      #sys.exit("dt is too high")
  if dt > beta * np.sqrt(qc):
    sys.exit("dt is high af")
      
  
  K = dt*b/2.0
  C_x = float(dt**2/dx**2)
  C_y = float(dt**2/dy**2)
  
  xv = x[:,np.newaxis]
  yv = y[np.newaxis,:]
  #xv = x.reshape((x.size,1))
  #yv = y.reshape((1,y.size))
  #Init wave
  u_2[1:-1,1:-1] = u_exact(xv,yv,0) # inner points
  
  u_2[-1,1:-1] = u_2[-3,1:-1]   # update ghost points
  u_2[1:-1,0] = u_2[1:-1,2]
  u_2[1:-1,-1] = u_2[1:-1,-3]
  u_2[0,1:-1] = u_2[2,1:-1]
  
  #plots 
  
  """fig = plt.figure()
  ax = Axes3D(fig)
  X, Y = np.meshgrid(xv, yv)
  ax.plot_surface(X, Y, u_2[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
  ax.auto_scale_xyz([0, Lx], [0, Ly], [-1, 1])
  plt.show()"""
  fig = plt.figure()
  plt.ioff()
  x,y = np.meshgrid(x,y)
  ax = Axes3D(fig)
  ax.plot_surface(x, y, u_2[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
  ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
  plt.savefig('waveplot_%04d.png' % (0))
  #plt.show()

  
  
  
  
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
  
  """
  fig = plt.figure()
  ax = Axes3D(fig)
  X, Y = np.meshgrid(xv, yv)
  ax.plot_surface(X, Y, u_1[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
  ax.auto_scale_xyz([0, Lx], [0, Ly], [-1, 1])
  plt.show()"""
  
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
    """
    fig = plt.figure()
    ax = Axes3D(fig)  
    X, Y = np.meshgrid(xv, yv)
    ax.plot_surface(X, Y, u[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
    #ax.plot_surface(X, Y, u_exact(xv,yv,t[n]), rstride=1, cstride=1, cmap='hot')
    ax.auto_scale_xyz([0, Lx], [0, Ly], [-1, 1])
    #plt.savefig('frame_%04d.png'%int(n))
    plt.show()
    """
    """
    fig = plt.figure()
    ax = Axes3D(fig)  
    X, Y = np.meshgrid(xv, yv)
    ax.plot_surface(X, Y, u_exact(X,Y,t[n]), rstride=1, cstride=1, cmap='hot')
    ax.auto_scale_xyz([0, Lx], [0, Ly], [-1, 1])
    plt.savefig('frame_%04d.png'%int(n))
    #plt.show()"""
    fig = plt.figure()
    plt.ioff()
    ax = Axes3D(fig)
    ax.plot_surface(x, y, u[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
    ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
    plt.savefig('waveplot_%04d.png' % (n))
    #plt.show()
      
    
    
    
    #update time steps
    #u_2[:,:] = u_1[:,:] ; u_1[:,:] = u[:,:]
    u_2, u_1, u = u_1, u, u_2  #faster 
    error = abs(u_exact(xv,yv,t[n])- u[1:-1,1:-1])
    error_list[n] = np.amax(error)
    
  error_mean = np.amax(error_list)
    
    
  return error_mean 



def standing_wave():
  
  dt=0.01
  Lx = 1
  Ly = 1
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
  u_exact = lambda x, y, t: A*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)*np.cos(w*t)  
  print solver(dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)
  
  """
  Fx = 10	
  Fy = 10
  error = []
  h_list = [0.1,0.01,0.001,0.0005]
  for h in h_list:
	  dt = h
	  Nx = int(round(Lx/(Fx*h)))
	  Ny = int(round(Ly/(Fy*h)))

	  error.append(solver(dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V))
  print "The error is:" 
  print error

  r = [np.log(error[i-1]/error[i])/ \
	  np.log(h_list[i-1]/h_list[i]) \
		  for i in range(1, len(h_list), 1)]
  print "... with a convergence rate of:" 
  print r
  """
    
  
  
  


if __name__ == "__main__":
  standing_wave()