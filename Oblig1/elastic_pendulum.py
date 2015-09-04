from matplotlib.pylab import *
import math


def simulate(
    beta=0.9,                 # dimensionless parameter
    Theta=30,                 # initial angle in degrees
    epsilon=0,                # initial stretch of wire
    num_periods=6,            # simulate for num_periods
    time_steps_per_period=60, # time step resolution
    Plot=True,                # make plots or not
    Vertical = True 
    ):
	P = 2*pi
	T = P * num_periods
	dt = P / time_steps_per_period 
	n = num_periods*time_steps_per_period
	x = zeros(n)
	y = zeros(n)
	time = linspace(0,T,n)
	theta = zeros(n)
	theta[0] = Theta
	if Vertical == True: # test used only for y-values, and x values are zero
		x[0] = 0 
	else:
		x[0] = (1-epsilon)*sin(theta[0])

	y[0] = 1-(1+epsilon)*cos(theta[0])
	L = sqrt(x[0]**2+(y[0]-1)**2)
	if Vertical == True:
		x[1] = 0 
	else:
		x[1] = 0.5*(2*x[0] - (beta/(1-beta))*(1-(beta/L))*x[0]*dt**2)

	y[1] = 0.5*(2*y[0] - dt**2*(beta/(1-beta))*(1-(beta/L))*(y[0]-1)-beta*dt**2)
	theta[1] = math.degrees(math.atan(x[1]/(y[1]-1)))
	for i in range(1,n-1):
		L = sqrt(x[i]**2+(y[i]-1)**2)
		if Vertical == True:
			x[i+1] = 0
		else:
			x[i+1] = 2*x[i] - x[i-1] + dt**2*(-beta/(1-beta))*(1-(beta/L))*x[i]
		y[i+1] = 2*y[i] -y[i-1] + dt**2*((-beta/(1-beta))  *  (1-(beta/L)) * (y[i]-1)-beta)
		theta[i+1] = math.degrees(math.atan(x[i+1]/(y[i+1]-1)))
	if Plot == True:
		if Vertical == False:
			figure()
			plot(x,y)
			gca().set_aspect('equal')
			show()
		
		if Theta < 10:
		
			figure()
			plot(time, Theta*cos(time),time,theta)
			show()	
		else:
			figure()
			plot(time,theta)
			title("Theta/time")
			show()

	return x,y,theta,time
def test_1():
	simulate(epsilon = 0, Theta = 0)
def test_Vertical():
	epsilon = 0.1
	beta = 0.9
	x,y,theta,time = simulate( Vertical = True, Theta= 0, epsilon = 0.1, time_steps_per_period = 600)
	y_exact = -epsilon*cos((sqrt(beta/(1-beta)))*time)
	error = abs(y_exact - y)  # computes the error term. 
	if max(error) < (epsilon/20):  # takes the highest value of all the error terms
		print "This is nice"
	else:
		print "This is shit"			
	figure()
	plot(time,error)
	title("plots the error between the calculated and exact vibration")
	show()

def demo(beta, Theta):
	simulate(beta , Theta , Vertical = False,time_steps_per_period=600,num_periods=6 )


if __name__ == '__main__':
	 print " Vertical Test "
	 test_Vertical()
	 print "Epsilon zero test"
	 test_1()
	 print "demo with delta = 0.9 and theta = 30 "
	 demo(0.9,30)
