#Importing the required modules
from pylab import * 
import numpy as np 


#initializing
N = 100 #Taking N=100, sections in loop antenna
a = 10 #Radius of the loop

#Defining phi dl and r'
phi = np.linspace(0,2*np.pi,N+1)[:-1] 
dl = np.c_[(-(a*2*np.pi)/N)*np.sin(phi),((a*2*np.pi)/N)*np.cos(phi),np.zeros(N)] #x,y,z components of the vector dl
r_ = np.c_[a*np.cos(phi),a*np.sin(phi),np.zeros(N)] #(x,y,z) coordinates of points on loop antenna and z=0 always

#plotting the points on the loop antenna
figure(1)
scatter(r_[:,0],r_[:,1]) # r_[:,0] is x coordinate and r_[:,1]is y coordinate
xlabel("x-axis")
ylabel("y-axis")
title("Points on loop antenna")
axis("square")
grid(True)
show()

#Defining I
I =  np.c_[-a*np.cos(phi)*np.sin(phi),a*np.cos(phi)*np.cos(phi),np.zeros(N)] #x,y,z components of current

#Plotting the current directions on loop antenna
figure(2)
quiver(r_[:,0],r_[:,1],I[:,0],I[:,1],scale=400,color='red') # to plot current with directions
xlabel("x-axis")
ylabel("y-axis")
title("Current elements at the centre points, on loop antenna")
axis("square")
grid(True)
xlim(-12,12)
ylim(-12,12)
show()


#Breaking the volume into 3 by 3 by 1000 mesh, with mesh points separated by 1cm
x = np.arange(-1,2,1)
z = np.arange(1,1001,1)
X,Y,Z = np.meshgrid(x,x,z)#meshgrid returns all points into space

#Defining the vector r
r = np.zeros((3,3,1000,3)) #Initializing array of size 3 by 3 by 1000 by 3
r[:,:,:,0]=X
r[:,:,:,1]=Y
r[:,:,:,2]=Z

#Function to calculate the vector potential due to lth element in the loop antenna
def calc(l):
    R = norm(r-r_[l],axis=-1)
    
    A_x = (np.cos(phi[l]))*np.exp(-0.1j*R)*dl[l,0]/R
    A_y = (np.cos(phi[l]))*np.exp(-0.1j*R)*dl[l,1]/R
    return A_x,A_y

#Initializing two arrays of size 3 by 3 by 1000 and the data type as complex
A_xx , A_yy = np.zeros((3,3,1000),dtype="complex"), np.zeros((3,3,1000),dtype="complex")

#for loop to add up the potentials due to every element 
for l in range(0,N):
    A_X,A_Y = calc(l) #calling the calc function
    A_xx += A_X       #Updating the x component of potential at every point 
    A_yy += A_Y       #Updating the y component of potential at every point
    

Bz = (A_yy[1,0,:]-A_yy[-1,0,:]-(A_xx[0,1,:]-A_xx[0,-1,:]))/4 # z component of the Magnetic field
z_=np.arange(1,1001,1) #creating z_ array of length 1000

#Plotting magnetic field along z-axis
figure(3)
loglog(z_,abs(Bz),label='|Magnetic field|')
legend()
xlabel("z-axis")
ylabel("Magnitude of Bz")
title("Magnitude of Bz in loglog plot")
grid(True)
show()

#Fitting Data   
K=0
y = np.log(abs(Bz))[K:]
x = np.log(z[K:])
A = np.c_[x,np.ones(1000)[K:]]
m,b = lstsq(A,y)[0] #To obtain best fit values

#Plotting both the field and fitting model in the same graph
figure(4)
plot(x,y,label="|Magnetic Field|")
print(m,b)
plot(x,m*x+b,label="Fitting model")
xlabel("log(z)")
ylabel("log(|Bz|)")
legend()
grid(True)
show()