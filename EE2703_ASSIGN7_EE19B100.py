from pylab import * #imports everything from the module pylab
import scipy.signal as sp 
import numpy as np

"""Q1"""
alpha=0.5
#Defining numerator and denominator using polyadd and polymul
num=np.poly1d([1,alpha],[0])
den=np.polymul([1,0,2.25],[1,2*alpha,2.25+alpha**2])

#Plotting the time response of 0.5 decay using system.impulse
figure(1)
t=linspace(0,50,2001)
t,y=sp.impulse(sp.lti(num,den),None,t)
plot(t,y,'r')
xlabel(r'$t$',size=20)
ylabel(r'$x(t)$',size=20)
title("Time response with decay of 0.5")
grid(True)
show()

"""Q2"""
alpha=0.05
num=np.poly1d([1,alpha])
den=np.polymul([1,0,2.25],[1,2*alpha,2.25+alpha**2])

#Plotting the time response of 0.05 decay using system.impulse
figure(2)
t=linspace(0,50,2001)
t,y=sp.impulse(sp.lti(num,den),None,t)
plot(t,y,'r')
xlabel(r'$t$',size=20)
ylabel(r'$x(t)$',size=20)
title("Time response with decay of 0.05")
grid(True)
show()

"""Q3"""
#Dividing the frequency from 1.4 to 1.6 in the step of 0.05
freq=linspace(1.4,1.6,5)

#Plotting the time respones of different frequency
figure(3)
title("Time response with different frequency")
xlabel('$t$',size=20)
ylabel('$x(t)$',size=20)
l=[]
#Using for loop & signal.lsim to calculate the responses of each frequency
for f in freq:
    H = sp.lti([1],[1,0,2.25])#Defining the transfer function
    t=linspace(0,70,1001)
    f_=np.cos(f*t)*np.exp(-0.05*t)*(t>0)
    t,x,svec=sp.lsim(H,f_,t)
    plot(t,x,linewidth=0.7)
    l.append("Freq ="+ str(f))
legend(l)
show()

"""Q4"""
t=linspace(0,20,1001)

#Defining x(t) & y(t)
t,x=sp.impulse(sp.lti([1,0,2],[1,0,3,0]),None,t)
t,y=sp.impulse(sp.lti([2],[1,0,3,0]),None,t)

#Plotting both x(t) & y(t) on the same figure
figure(4)
title("Responses of coupled system")
xlabel('$t$',size=20)
ylabel("Time response",size=20)
plot(t,x)
plot(t,y)
legend(["$x(t)$","$y(t)$"])
show()

"""Q5"""
#Defining the transfer function
H=sp.lti([1],[1e-12,1e-4,1])

#Bode plot of H(s)
w,S,phi=H.bode()

#Plotting magnitude of H(s) vs frequency in rad(log)
semilogx()
plot(w,S)
title("Magnitude plot")
xlabel("Frequency in rad/sec (log)",size=20)
ylabel("Magnitude in dB",size=20)
grid(True)
show()

#Plotting phase of H(s) vs frequency in rad(log)
semilogx()
plot(w,phi)
title("Phase plot")
xlabel("Frequency in rad/sec (log)",size=20)
ylabel("Phase in degrees",size=20)
grid(True)
show()

"""Q6"""
#Defining the input 
def input_2(t,w1=1e3,w2=1e6):
    """Two cosines of different frequencies"""
    u_t = 1*(t>0)
    return (np.cos(w1*t)-np.cos(w2*t)) * u_t

# early response
t1=np.linspace(0,float(30e-6),1001)
t1,y1,svec = sp.lsim(H,input_2(t1),t1)

# steady state response
t2=np.linspace(0,float(10e-3),1001)
t2,y2,svec = sp.lsim(H,input_2(t2),t2)

#Plotting the early response
figure(5)
title(r"Response for 30 $\mu s$")
xlabel("$t$ (sec)",size=20)
ylabel("$v_0(t)$",size=20)
plot(t1,y1)
grid(True)
show()

#Plotting the steady state response
figure(6)
title(r"Response for 10 msec")
xlabel("$t$ (sec)",size=20)
ylabel("$v_0(t)$",size=20)
plot(t2,y2)
grid(True)
show()
