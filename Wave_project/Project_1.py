import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sympy as sym
import sys
def f_sympy(u,qval,w):

  x,y,t,q,b,kx,ky,A,B,w = sym.symbols('x, y, t, q, b, kx, ky, A, B, w')
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

  return f


def solver(I,dt,dx,dy,Lx,Ly,T,beta,b,qc,f,u_exact,V):
  
  Nt = int(round(T/float(dt)))
  Nx = Lx/dx
  Ny = Ly/dy
  
  u = np.zeros((Nx+3,Ny+3))
  u_1 = np.zeros((Nx+3,Ny+3))
  u_2 = np.zeros((Nx+3,Ny+3))
  #f = np.zeros((Nx+1,Ny+1))
  error_list = np.zeros(Nt+1)


  
  x = np.linspace(0,Lx,Nx+1)
  y = np.linspace(0,Ly,Ny+1)
  t = np.linspace(0, T, Nt+1)

  q = np.zeros((Nx+3,Ny+3)) # mesh for q
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

  q[1:-1,1:-1] = qc(x,y)
  q[-1,1:-1] = q[-3,1:-1]   # update ghost points
  q[1:-1,0] = q[1:-1,2]
  q[1:-1,-1] = q[1:-1,-3]
  q[0,1:-1] = q[2,1:-1] 


  u_2[1:-1,1:-1] = I(xv,yv) # load initial condition

  
  u_2[-1,1:-1] = u_2[-3,1:-1]   # update ghost points
  u_2[1:-1,0] = u_2[1:-1,2]
  u_2[1:-1,-1] = u_2[1:-1,-3]
  u_2[0,1:-1] = u_2[2,1:-1]
  
  #plots 
  #error_list[0] = np.amax(error)
  #print "first error should be zero:...", error

  fig = plt.figure()
  plt.ioff()
  X,Y = np.meshgrid(x,y)
  ax = Axes3D(fig)
  ax.plot_surface(X, Y, u_2[1:-1,1:-1], rstride=1, cstride=1, cmap='hot',color='blue')
  ax.plot_wireframe(X, Y, u_exact(x,y), rstride=1, cstride=1, cmap='hot',color='red')
  ax.view_init(25,100)
  ax.auto_scale_xyz([0, Lx], [0, Ly], [0, 3])
  plt.title("dt is: %.3f, dx = dy = %.3f"%(dt,dx)) 
  plt.savefig('waveplot_%04d.png' % (0))
  #plt.show()
  
  
  
  #first step
  #inner points
  u_1[1:-1,1:-1] = u_2[1:-1,1:-1]+ dt*V*(2 +dt*b) + \
    C_x*0.25*( (q[2:,1:-1]+q[1:-1,1:-1])*(u_2[2:,1:-1]-u_2[1:-1,1:-1]) -\
      (q[1:-1,1:-1]+q[:-2,1:-1])*(u_2[1:-1,1:-1]- u_2[:-2,1:-1])  )\
	+ C_y*0.25*( (q[1:-1,2:] + q[1:-1,1:-1])*(u_2[1:-1,2:]-u_2[1:-1,1:-1]) - \
	    (q[1:-1,1:-1]+q[1:-1,:-2])*(u_2[1:-1,1:-1]-u_2[1:-1,:-2])   ) 
  #+0.5*dt**2*f(xv[:],yv[:],t[1])
  
  

  u_1[-1,1:-1] = u_1[-3,1:-1]   # update ghost points in first time step
  u_1[1:-1,0] = u_1[1:-1,2]
  u_1[1:-1,-1] = u_1[1:-1,-3]
  u_1[0,1:-1] = u_1[2,1:-1]

  #error = abs(u_exact(x,y,t[1])- u_1[1:-1,1:-1])
  #error_list[1] = np.amax(error)
  fig = plt.figure()
  plt.ioff()
  X,Y = np.meshgrid(x,y)
  ax = Axes3D(fig)
  ax.plot_wireframe(X, Y, u_1[1:-1,1:-1], rstride=1, cstride=1, cmap='hot',color='blue')
  ax.plot_wireframe(X, Y, u_exact(x,y), rstride=1, cstride=1, cmap='hot',color='red')
  ax.view_init(25,100)
  ax.auto_scale_xyz([0, Lx], [0, Ly], [0, 3])
  plt.title("dt is: %.3f, dx = dy = %.3f"%(dt,dx)) 
  plt.savefig('waveplot_%04d.png' % (0))
  #plt.show()


  
  
  
  # all time steps
  for n in range(2,Nt):
    # inner points
    
    u[1:-1,1:-1] = ((K-1)/(K+1))*u_2[1:-1,1:-1] + (2/(1+K))*u_1[1:-1,1:-1] + \
	    (0.5*C_x/(K+1))*((q[2:,1:-1]+q[1:-1,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
		    (q[1:-1,1:-1]+q[0:-2,1:-1])*(u_1[1:-1,1:-1]-u_1[0:-2,1:-1])) + \
			    (0.5*C_y/(K+1)) * ((q[1:-1:,2:]+q[1:-1,1:-1])*(u_1[1:-1,2:]-u_1[1:-1,1:-1]) - \
				    (q[1:-1:,1:-1]+q[1:-1,0:-2])*(u_1[1:-1,1:-1]-u_1[1:-1,0:-2])) 
          #dt**2*f(xv[:],yv[:],t[n])
					    
	
    
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
      X,Y = np.meshgrid(x,y)
      ax = fig.add_subplot(111, projection='3d') 
      ax.plot_wireframe(X,Y,u[1:-1,1:-1],rstride=1, cstride=1, cmap='hot',color='blue') 
      ax.plot_wireframe(X,Y,u_exact(x,y), rstride=1, cstride=1, cmap='hot',color='red')
      ax.auto_scale_xyz([0, Lx], [0, Ly], [0, 3])
      ax.view_init(25,100)
      plt.title("dt is: %.3f, dx = dy = %.3f"%(dt,dx)) 
      plt.savefig('waveplot_%04d.png' % (n))
                #plt.show()
      
    
    
    
    #update time steps
    u_2[:,:] = u_1[:,:] ; u_1[:,:] = u[:,:]
    #u_2, u_1, u = u_1, u, u_2  #faster 
    """error = abs(u_exact(x[:],y[:],t[n])- u[1:-1,1:-1])
                error_list[n] = np.amax(error)
              error_max = np.max(error_list)
                """
    
  return None



def manufactured():

  def make_f(b_val,kx_val,ky_val,A_val,B_val,w_val):
      x,y,t,q,b,kx,ky,A,B,w = sym.symbols('x ,y, t, q, b, kx, ky, A, B, w')
      q = x+1
      w = sym.sqrt(q*(kx**2+ky**2))
      u = (A*sym.cos(w*t)+B*sym.sin(w*t))*sym.exp(-sym.sqrt(q)*t)*sym.cos(kx*x)*sym.cos(ky*y)
      f = f_sympy(u,q,w)
      f = f.subs([(b,b_val),(kx,kx_val),(ky,ky_val),(A,A_val),(B,B_val),(w,w_val)])
      return np.vectorize(sym.lambdify((x,y,t),f))

  dx =0.03
  dy =0.03
  #q = 1
  b = 1
  Lx = 1
  Ly = 1
  A = 1
  B = 1
  beta = 0.9
  mx = 1.0
  my = 1.0
  kx = mx*np.pi/Lx
  ky = my*np.pi/Ly
  T = 6.0
  V = 0.0
  qc = lambda x,y : x+1 + y
  w = np.sqrt(qc(1,1)*(kx**2+ky**2))
  #print f()
  dt = float(np.sqrt(1.0/qc(1,1))*1/np.sqrt((1.0/(dx))**2+(1.0/(dy))**2))

  #from math import sqrt, cos, sin, exp
  f = make_f(b,kx,ky,A,B,w)

  u_exact = np.vectorize(lambda x, y, t: (A*np.cos(w*t)+B*np.sin(w*t))*np.exp(-np.sqrt(qc(x,y))*t)*np.cos(kx*x)*np.cos(ky*y))  
  I = np.vectorize(lambda x,y: A*np.cos(kx*x)*np.cos(ky*y))
  
  
  print solver(I,dt,dx,dy,Lx,Ly,T,beta,b,qc,f,u_exact,V)

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
  Lx = 3.0
  Ly = 3.0
  dx = 0.05
  dy = 0.05
  T = 5
  qc = lambda x,y: 2
  u_exact = lambda x,y: 0
  beta = 0.8
  f = 0
  V  = 0 
  dt = np.sqrt(1.0/qc(1,1))*1/np.sqrt((1.0/dx)**2+(1.0/dy)**2)
  I = 0

  #dt = (Lx/Nx)/(np.sqrt(qc))
  print "dt is:... ", dt, dx ,dy
  if direction == "x":
    I = np.vectorize(lambda x,y: 0 if abs(x-0.5*Lx) > 0.1 else 2) 
  if direction == "y":
    I = np.vectorize(lambda x,y: 0 if abs(y-0.5*Ly) > 0.1 else 2)
  if direction == "middle":
    I = np.vectorize(lambda x,y: 0 if abs(x-0.5*Lx) > 0.05 or abs(y-0.5*Ly) > 0.05 else 2)

  solver(I,dt,dx,dy,Lx,Ly,T,beta,b,qc,f,u_exact,V)
def standing_wave():
  
  Lx = 1.0
  Ly = 1.0
  #dx = 0.05
  #dy = 0.05
  T=2
  mx = 1
  my = 1
  kx = mx*np.pi/Lx
  ky = my*np.pi/Ly
  #w = np.sqrt(9.81*kx)
  A=1
  V = 0
  beta = 0.9
  b = 0
  qc = lambda x,y: 2
  f = lambda x,y,t : 0
  w = np.sqrt(qc(1,1)*(kx**2+ky**2))

  u_exact = np.vectorize(lambda x, y, t: A*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly)*np.cos(w*t))  
  I = np.vectorize(lambda x,y: A*np.cos(mx*np.pi*x/Lx)*np.cos(my*np.pi*y/Ly))
  #dx = Lx/(1/h)
  #dy = Lx/(1/h)
  #dt = h/10
  #print dx ,dy ,dt
  #dt = np.sqrt(1.0/qc(1,1))*1/np.sqrt((1.0/(dx))**2+(1.0/(dy))**2)
  rounds = 8
  h = 0.5
  h_values = np.zeros(rounds+1)
  err = np.zeros(rounds+1)
  for i in range(rounds):
    dt = h/10.0
    Nx = 1.0/h
    Ny = 1.0/h
    dx = Lx/(Nx)
    dy = Ly/(Ny)
    h_values[i] = h
    err[i] = solver(I,dt,dx,dy,Lx,Ly,T,beta,b,qc,f,u_exact,V)
    print err[i]
    #print "Error is  %.25f with dt: %f, dx:%.4f , dy:%.4f " %(solver(I,dt,dx,dy,Lx,Ly,T,beta,b,qc,f,u_exact,V), dt,dx,dy)
    h = h/2.0
  r = [np.log(err[i-1]/err[i])/ \
    np.log(h_values[i-1]/h_values[i]) \
      for i in range(1, len(h_values), 1)]
  print r
  """
  h = [0.1,0.01,0.001]
  for i in h:
    error = solver(I,dt,Nx,Ny,Lx,Ly,T,beta,b,qc,f,u_exact,V)  
  """  

def Hill(choice):


  b = 0
  Lx = 3.0
  Ly = 3.0
  dx = 0.04
  dy = 0.04
  T = 4
  #qc = 1
  u_exact = 0
  beta = 0.8
  f = 0
  V  = 0 

  x = np.linspace(0,Lx,(Lx/dx)+1)
  y = np.linspace(0,Ly,(Ly/dy)+1)
  xv = x[:,np.newaxis]
  yv = y[np.newaxis,:]

  #dt = (Lx/Nx)/(np.sqrt(qc))
  if choice == "b":

    I_0 = 1.5; I_a = 1.7 ; I_m = 0; I_s = 0.5
    B_0 = 0.0; B_a = 1.5; B_mx = Lx/2.0 ; B_my = Ly/2.0 ; B_s = 1 ; b_scale = 0.5
    I = np.vectorize(lambda x,y: I_0 + I_a*np.exp(-((y-I_m)/I_s)**2))
    B = np.vectorize(lambda x,y: B_0 + B_a*np.exp(-((x-B_mx)/B_s)**2 - ((y-B_my)/(b_scale*B_s))**2))
  if choice == "h":

    I_0 = 1.0; I_a = 1.0 ; I_m = Lx ; I_s = 0.5
    B_0 = 1; B_a = 1.0; B_mx =Lx/2.  ; B_my = Lx/2. ; B_s = 0.5
    I = np.vectorize(lambda x,y: I_0 + I_a*np.exp(-((x-I_m)/I_s)**2))
    B = np.vectorize(lambda x,y: B_0 + B_a*np.cos(np.pi*(x-B_mx)/(2*B_s))*np.cos(np.pi*(y-B_my)/(2*B_s))\
               if 0 <= np.sqrt((x-Lx/2.0)**2+(y-Ly/2.0)**2) <= B_s else B_0 )

  q = np.vectorize(lambda x,y:  9.81 * (I(x,y)-B(x,y)) )
  dt = np.sqrt(1.0/max(abs(q(x,y))))*1/np.sqrt((1.0/(dx))**2+(1.0/(dy))**2)
  fig = plt.figure()
  plt.ioff()
  X,Y = np.meshgrid(x,y)
  ax = Axes3D(fig)
  ax.plot_wireframe(X, Y, I(x,y), rstride=1, cstride=1, cmap='hot',color='green')
  ax.plot_wireframe(X, Y, B(x,y), rstride=1, cstride=1, cmap='hot',color='red')
  ax.view_init(25,100)
  ax.auto_scale_xyz([0, Lx], [0, Ly], [0, 2])
  plt.title("dt is: %.3f, dx = dy = %.3f"%(dt,dx)) 
  #plt.savefig('waveplot_%04d.png' % (0))
  plt.show()
  
  solver(I,dt,dx,dy,Lx,Ly,T,beta,b,q,f,B,V)
  



if __name__ == "__main__":
  #plug_wave(raw_input("x,y, or middle"))
  #standing_wave()
  #manufactured()
  Hill(raw_input("h for hat or b for beach........"))