import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sympy as sym
import sys
def f_sympy(u,qval):

  x,y,t,q,b,kx,ky,A,B,w = sym.symbols('x, y, t, q, b, kx, ky, A, B, w')
  u = u
  q = qval
  u_tt = sym.diff(u,t,t)
  u_t = sym.diff(u,t)
  u_x = sym.diff(u,x)
  u_y = sym.diff(u,y)
  q_x = sym.diff(q,x) 
  q_y = sym.diff(q,y)
  u_xx = sym.diff(u,x,x)
  u_yy = sym.diff(u,y,y)
  f =u_tt + b*u_t -q_x*u_x - q*u_xx - q_y*u_y - q*u_yy 

  return sym.simplify(f)


def solver(I,dt,dx,dy,Lx,Ly,T,beta,b,qc,f,u_exact,V):
  
  Nt = int(round(T/float(dt)))
  Nx = Lx/dx
  Ny = Ly/dy
  
  u = np.zeros((Nx+3,Ny+3))
  u_1 = np.zeros((Nx+3,Ny+3))
  u_2 = np.zeros((Nx+3,Ny+3))
  f = np.zeros((Nx+1,Ny+1))
  error_list = np.zeros(Nt+1)

  q = np.zeros((Nx+3,Ny+3))+qc[x]  # mesh for q

  
  x = np.linspace(0,Lx,Nx+1)
  y = np.linspace(0,Ly,Ny+1)
  t = np.linspace(0, T, Nt+1)
  
  #dx = float(x[1]-x[0])
  #dy = float(y[1]-y[0])
  #dx = Lx/Nx
  #dy = Ly/Ny
  #rint "dt is %.5f.. dx is %.5f .. dy is %.5f .." %(dt,dx,dy)
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
  """Ix = range(1, u.shape[0]-1)
  Iy = range(1, u.shape[1]-1)
  It = range(0, t.shape[0])    
        # Fetching the initial condition
  for j in Iy:
    for i in Ix:
      u_2[i,j]=I(x[i-1],y[j-1])"""
  #I = np.vectorize(I(xv,yv))

  u_2[1:-1,1:-1] = I(xv,yv)

  
  u_2[-1,1:-1] = u_2[-3,1:-1]   # update ghost points
  u_2[1:-1,0] = u_2[1:-1,2]
  u_2[1:-1,-1] = u_2[1:-1,-3]
  u_2[0,1:-1] = u_2[2,1:-1]
  
  #plots 
  error = abs(u_exact(xv,yv,t[0])- u_2[1:-1,1:-1])
  error_list[0] = np.amax(error)
  #print "first error should be zero:...", error

  fig = plt.figure()
  plt.ioff()
  x,y = np.meshgrid(x,y)
  ax = Axes3D(fig)
  ax.plot_surface(xv, yv, u_2[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')
  ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
  plt.savefig('waveplot_%04d.png' % (0))


  
  
  
  
  #first step
  #inner points
  u_1[1:-1,1:-1] = u_2[1:-1,1:-1]+ dt*V*(2 +dt*b) + \
    C_x*0.25*( (q[2:,1:-1]+q[1:-1,1:-1])*(u_2[2:,1:-1]-u_2[1:-1,1:-1]) -\
      (q[1:-1,1:-1]+q[:-2,1:-1])*(u_2[1:-1,1:-1]- u_2[:-2,1:-1])  )\
	+ C_y*0.25*( (q[1:-1,2:] + q[1:-1,1:-1])*(u_2[1:-1,2:]-u_2[1:-1,1:-1]) - \
	    (q[1:-1,1:-1]+q[1:-1,:-2])*(u_2[1:-1,1:-1]-u_2[1:-1,:-2])   ) +0.5*f(x[:],y[:],t[0])
  
  

  u_1[-1,1:-1] = u_1[-3,1:-1]   # update ghost points in first time step
  u_1[1:-1,0] = u_1[1:-1,2]
  u_1[1:-1,-1] = u_1[1:-1,-3]
  u_1[0,1:-1] = u_1[2,1:-1]

  error = abs(u_exact(x,y,t[1])- u_1[1:-1,1:-1])
  error_list[1] = np.amax(error)
  fig = plt.figure()
  plt.ioff()
  ax = Axes3D(fig)
  ax.plot_surface(x, y, u_1[1:-1,1:-1], rstride=1, cstride=1, cmap='hot')  
  ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
  plt.savefig('waveplot_%04d.png' % (1))
  #plt.show()
  
  
  
  # all time steps
  for n in range(2,Nt):
    # inner points
    
    u[1:-1,1:-1] = ((K-1)/(K+1))*u_2[1:-1,1:-1] + (2/(1+K))*u_1[1:-1,1:-1] + \
	    (0.5*C_x/(K+1))*((q[2:,1:-1]+q[1:-1,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
		    (q[1:-1,1:-1]+q[0:-2,1:-1])*(u_1[1:-1,1:-1]-u_1[0:-2,1:-1])) + \
			    (0.5*C_y/(K+1)) * ((q[1:-1:,2:]+q[1:-1,1:-1])*(u_1[1:-1,2:]-u_1[1:-1,1:-1]) - \
				    (q[1:-1:,1:-1]+q[1:-1,0:-2])*(u_1[1:-1,1:-1]-u_1[1:-1,0:-2])) + f(x[:],y[:],t[n])
					    
	
    
    u[-1,1:-1] = u[-3,1:-1]   # update ghost points in first time step
    u[1:-1,0] = u[1:-1,2]
    u[1:-1,-1] = u[1:-1,-3]
    u[0,1:-1] = u[2,1:-1]
    
    
    #plot and make figure
    
    """num_frames = 100
                skip_frame = int(Nt/float(num_frames))
                if n % skip_frame == 0 or n == Nt-1:
                  fig = plt.figure()
                  plt.ioff() 
                  ax = fig.add_subplot(111, projection='3d') 
                  ax.plot_wireframe(xv,yv,u[1:-1,1:-1],rstride=1, cstride=1, cmap='hot') 
                  ax.plot_wireframe(xv,yv,u_exact(xv,yv,t[n]), rstride=1, cstride=1, cmap='hot')
                  ax.auto_scale_xyz([0, Lx], [0, Ly], [-2, 2])
                  plt.savefig('waveplot_%04d.png' % (n))"""
                #plt.show()
      
    
    
    
    #update time steps
    u_2[:,:] = u_1[:,:] ; u_1[:,:] = u[:,:]
    #u_2, u_1, u = u_1, u, u_2  #faster 
    error = abs(u_exact(x[:],y[:],t[n])- u[1:-1,1:-1])
    error_list[n] = np.amax(error)
  error_max = np.max(error_list)
    
    
  return error_max 



def manufactured():

  def make_f(b_val,kx_val,ky_val,A_val,B_val,w_val):
      x,y,t,q,b,kx,ky,A,B,w = sym.symbols('x ,y, t, q, b, kx, ky, A, B, w')
      q = 1 + x
      u = (A*sym.cos(w*t)+B*sym.sin(w*t))*sym.exp(-sym.sqrt(q)*t)*sym.cos(kx*x)*sym.cos(ky*y)
      f = f_sympy(u,q)
      f = f.subs([(b,b_val),(kx,kx_val),(ky,ky_val),(A,A_val),(B,B_val),(w,w_val)])
      print f
      return sym.lambdify((x,y,t),f)

  dx =0.05
  dy =0.05
  q = 1
  b = 0
  Lx = 1
  Ly = 1
  A = 1
  B = 1
  beta = 0.9
  mx = 1.0
  my = 1.0
  kx = mx*np.pi/Lx
  ky = my*np.pi/Ly
  qc = lambda x,y : 1+x
  w = np.sqrt(qc(1,1)*(kx**2+ky**2))
  f = make_f(b,kx,ky,A,B,w)
  print f(1,1,1)
  dt = np.sqrt(1/qc(1,1))*1/np.sqrt((1/(dx))**2+(1/(dy))**2)

  u_exact = np.vectorize(lambda x, y, t: (A*np.cos(w*t)+B*np.sin(w*t))*np.exp(-np.sqrt(q)*t)*np.cos(kx*x)*np.cos(ky*y))  
  I = np.vectorize(lambda x,y: A**np.exp(-np.sqrt(q)*t)*np.cos(kx*x)*np.cos(ky*y))


def constant_solution():
  dt=0.01
  Lx = 100
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
  u_exact = lambda x, y, t: 2 
  error=solver(dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)







def plug_wave(direction):
  b = 1
  Lx = 1.0
  Ly = 1.0
  Nx = 100
  Ny = 100
  T = 10
  qc = 1
  u_exact = 0
  beta = 0.8
  f = 0
  V  = 0 
  dt = np.sqrt(1/qc)*1/np.sqrt((1/(Lx/Nx))**2+(1/(Ly/Ny))**2)

  #dt = (Lx/Nx)/(np.sqrt(qc))
  print "dt is:... ", dt
  if direction == "x":
    I = np.vectorize(lambda x,y: 0 if abs(x-0.5*Lx) > 0.1 else 2) 
  if direction == "y":
    I = np.vectorize(lambda x,y: 0 if abs(y-0.5*Ly) > 0.1 else 2)
  if direction == "middle":
    I = np.vectorize(lambda x,y: 0 if abs(x-0.1*Lx) > 0.0 or abs(y-0.9*Ly) > 0.05 else 2)
  solver(I,dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)

def standing_wave():
  
  Lx = 1.0
  Ly = 1.0
  dx = 0.5
  dy = 0.5
  T=2
  mx = 3
  my = 3
  kx = mx*np.pi/Lx
  ky = my*np.pi/Ly
  #w = np.sqrt(9.81*kx)
  A=1
  V = 0
  beta = 0.9
  b = 0
  qc = 2
  f = 0
  w = np.sqrt(qc*(kx**2+ky**2))


  u_exact = np.vectorize(lambda x, y, t: A*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)*np.cos(w*t))  
  I = np.vectorize(lambda x,y: A*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly))
  
  dt = np.sqrt(1.0/qc)*1/np.sqrt((1.0/(dx))**2+(1.0/(dy))**2)
  #solver(I,dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)
  
  for i in range(1,5):
    print "Error is  %.25f with dt: %f, dx:%.4f , dy:%.4f " %(solver(I,dt,dx,dy,Lx,Ly,T,beta,b,qc,f,u_exact,V), dt,dx,dy)
    #dt = dt/2.0
    #dx = np.sqrt(qc)*np.sqrt(2)*dt
    #dy = dx
    dx = dx/2.0
    dy = dy/2.0
    dt = np.sqrt(1.0/qc)*1/np.sqrt((1.0/(dx))**2+(1.0/(dy))**2)

    

  """
  h = [0.1,0.01,0.001]
  for i in h:
    error = solver(I,dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)  
  """  


if __name__ == "__main__":
  #plug_wave(raw_input("x,y, or middle"))
  #standing_wave()
  manufactured()