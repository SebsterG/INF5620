import sympy as sym
V, t, I, w, dt, b,c = sym.symbols('V t I w dt b c')  # global symbols 
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = DtDt(u,dt) + w**2*u(t)   - ode_source_term(u)
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    f_0 =  w**2*I + DtDt(u,dt).subs(t,0)
    exact = (ode_source_term(u).subs(t,dt) - sym.diff(u(t),t,t).subs(t,dt))/w**2
    R = exact - 0.5*(2*dt*V + dt**2*f_0 - dt**2*w**2*I + 2*I)
    return sym.simplify(R)

def DtDt(u, dt): 
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) - 2*u(t) + u(t-dt))/dt**2

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u
    print "Initial conditions u(0)=%s, u'(0)=%s:" % \
          (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

def linear():
    main(lambda t: V*t + I)
def quadratic():
    main(lambda t : b*t**2+V*t+I)

def poly():
    main(lambda t : c*t**3 + b*t**2+V*t+I)

from numpy import *
from matplotlib.pyplot import *
def solver(I1,w,stopp,deltat):
    I1 = float(I1)
    n = int(round(stopp/deltat))
    u = zeros(n+1)
    tim = linspace(0,5,n+1)
    u[0] = I1
    u[1] = u[0]+deltat
    for i in range(1,len(u)-1):
        #f_last = w**2*u[i]*deltat**2
        f_last = 0
        u[i+1] = deltat**2*f_last - w**2*deltat**2*u[i] + 2*u[i] - u[i-1]
    return tim, u  

def test (I1,w,stopp,deltat,V1):
    tim, u = solver(I1,w,stopp,deltat)   
    figure(1)
    title("myplot")
    plot (tim,u,"*")
    plot(tim, I1*cos(w*tim)+(V1/w)*sin(w*tim))
    show()



    
if __name__ == '__main__':
    """  linear()
          quadratic()"""
    test(1,1,5,0.1,1)
    """
                print "linear"
                linear()
                print "quadratic"
                quadratic()
                print "poly"
                poly() """
