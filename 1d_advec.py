import numpy as np
import matplotlib.pyplot as plt
N=100
#N includes ghost cells
ghost=1
u=1. #speed of prop.
C=.5 #CFL number
xlow=0.
xup=2*np.pi
dx = (xup-xlow)/(N)

dt = C*dx/u #timestep set by CFL condition


grid=np.linspace(xlow,xup,N)

a=np.zeros(np.shape(grid))
newa = np.zeros(np.shape(grid))

#init loop
for i in range(ghost,N-ghost):
    if (N/3. < i and i < 2*N/3.):
        a[i]=1
 #ghost zones are initialized properly in this setup



#zeroth and last index in a are just our ghost zones

myt = 0
stride =2
count=0
while myt<10:
    myt+=dt
    count+=1

    for i in range(ghost,N-ghost):
        update=-C*dt*u*((a[i]-a[i-1])/dx) #updwind spatial deriv

        newa[i]=a[i]+update    
   
    a=newa+0 #just plus 0 to dereference a from newa
    #update ghost zones
    a[0]=a[-(1+ghost)]
    a[-1]=a[ghost]

    if count%10==0:
        plt.plot(a)
        plt.savefig("advec_out/{}.png".format(myt))

