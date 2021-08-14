#importing the required modules
from pylab import *
import numpy as np


def dft(start,stop,steps,f,titl,y_label1,y_label2,x_label,x_lim,name,ro=False,go=False,phase_shift=False):
    #after phase wrapping
    if(phase_shift):
        #finding FFT
        sampling_rate = steps/(stop-start)
        x=linspace(start,stop,steps+1)[:-1]
        y = f(x)
        Y=fftshift(fft(y))/float(steps)
        w=sampling_rate*(linspace(-np.pi,np.pi,steps+1)[:-1])
        #plotting
        figure()
        subplot(2,1,1)
        plot(w,abs(Y),lw=2)
        xlim([-x_lim,x_lim])
        ylabel(y_label1,size=16)
        title(titl)
        grid(True)
        subplot(2,1,2)
        if (ro):
           plot(w,np.angle(Y),'ro',lw=2)
        if(go):
           ii=np.where(abs(Y)>1e-3)
           plot(w[ii],np.angle(Y[ii]),'go',lw=2)
        xlim([-x_lim,x_lim])
        ylabel(y_label2,size=16)
        xlabel(x_label,size=16)
        grid(True)
        savefig(name)
        show()
        return
    #without phase shift
    else:
        #findind FFT
        x=linspace(start,stop,steps)
        y=f(x)
        Y=fft(y)
        #plotting
        figure()
        subplot(2,1,1)
        plot(abs(Y),lw=2)
        xlim([0,x_lim])
        ylabel(y_label1,size=16)
        title(titl)
        grid(True)
        subplot(2,1,2)
        plot(np.unwrap(np.angle(Y)),lw=2)
        xlim([0,x_lim])
        ylabel(y_label2,size=16)
        xlabel(x_label,size=16)
        grid(True)
        savefig(name)
        show()
        return


#defining inputs

def f1(x):
    return np.sin(5*x)

def f2(x):
    return (1+0.1*np.cos(x))*np.cos(10*x)

def f3(x):
    return (np.sin(x))**3

def f4(x):
    return (np.cos(x))**3

def f5(x):
    return np.cos(20*x + 5*np.cos(x))


#Question1

#example1
x=rand(100)
X=fft(x)
y=ifft(X)
np.c_[x,y]
print("Absolute maximum error = ",abs(x-y).max())

#example2
#Without phase shifting
dft(0,2*np.pi,128,f1,r"Spectrum of $\sin(5t)$ without phase shift",r"$|Y|$",r"Phase of $Y$",r"$k$",140,"fig9.1.png",ro=False,go=False,phase_shift=False)
#After Phase shifting
dft(0,2*np.pi,128,f1,r"Spectrum of $\sin(5t)$ after phase shift",r"$|Y|$",r"Phase of $Y$",r"$k$",10,"fig9.2.png",ro=True,go=True,phase_shift=True)


#example3
dft(0,2*np.pi,128,f2,r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$",r"$|Y|$",r"Phase of $Y$",r"$\omega$",15,"fig9.3.png",ro=True,go=True,phase_shift=True)

dft(-4*np.pi,4*np.pi,512,f2,r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$",r"$|Y|$",r"Phase of $Y$",r"$\omega$",15,"fig9.4.png",ro=True,go=True,phase_shift=True)


#Question2

dft(-4*np.pi,4*np.pi,512,f3,r"Spectrum of $sin^3(t)$",r"$|Y|$",r"Phase of $Y$",r"$\omega$",15,"fig9.5.png",ro=False,go=True,phase_shift=True)

dft(-4*np.pi,4*np.pi,512,f4,r"Spectrum of $cos^3(t)$",r"$|Y|$",r"Phase of $Y$",r"$\omega$",15,"fig9.6.png",ro=False,go=True,phase_shift=True)

#Question3

dft(-4*np.pi,4*np.pi,512,f5,r"Spectrum of $cos(20t + 5cos(t))$",r"$|Y|$",r"Phase of $Y$",r"$\omega$",30,"fig9.5.png",ro=False,go=True,phase_shift=True)

#Question4
#defining gaussian and its expected CTFT 
def gauss(x):
    return np.exp(-0.5*x**2)

def expectedgauss(w):
    return 1/np.sqrt(2*np.pi) * np.exp(-w**2/2)

def estdft(tolerance=1e-6,samples=128,func = gauss,expectedfn = expectedgauss,wlim = 5):
    T = 8*np.pi
    N = samples
    Yold=0
    err=tolerance+1
    iters = 0
    #iterative loop to find window size
    while err>tolerance:  
        x=linspace(-T/2,T/2,N+1)[:-1]
        w = linspace(-N*np.pi/T,N*np.pi/T,N+1)[:-1]
        y = func(x)
        Y=fftshift(fft(ifftshift(y)))*T/(2*np.pi*N)
        err = sum(abs(Y[::2]-Yold))
        Yold = Y
        iters+=1
        T*=2
        N*=2
        

    #calculating error
    true_error = sum(abs(Y-expectedfn(w)))
    print("True error: ",true_error)
    print("samples = "+str(N)+" time period = pi*"+str(T/np.pi))

    mag = abs(Y)
    phi = np.angle(Y)
    phi[np.where(mag<tolerance)]=0
    
    # plot estimate
    figure()
    subplot(2,1,1)
    plot(w,abs(Y),lw=2)
    xlim([-wlim,wlim])
    ylabel('Magnitude',size=16)
    title("Estimate FFT of gaussian")
    grid(True)
    subplot(2,1,2)
    plot(w,np.angle(Y),'ro',lw=2)
    ii=np.where(abs(Y)>1e-3)
    plot(w[ii],np.angle(Y[ii]),'go',lw=2)
    xlim([-wlim,wlim])
    ylabel("Phase",size=16)
    xlabel("w",size=16)
    grid(True)
    show()

    #plotting expected output    
    Y_ = expectedfn(w)
    
    mag = abs(Y_)
    phi = np.angle(Y_)
    phi[np.where(mag<tolerance)]=0
    
    figure()
    subplot(2,1,1)
    plot(w,abs(Y),lw=2)
    xlim([-wlim,wlim])
    ylabel('Magnitude',size=16)
    title("True FFT of gaussian")
    grid(True)
    subplot(2,1,2)
    plot(w,np.angle(Y),'ro',lw=2)
    ii=np.where(abs(Y)>1e-3)
    plot(w[ii],np.angle(Y[ii]),'go',lw=2)
    xlim([-wlim,wlim])
    ylabel("Phase",size=16)
    xlabel("w",size=16)
    grid(True)
    show()

    return


estdft()