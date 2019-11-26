import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

#taking care of MPI stuff


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

comm = MPI.COMM_WORLD
size = comm.Get_size() #number of cores we call in mpiexec
rank = comm.Get_rank()

#need to divide data in array over different processors, for now let's assume that our number of cells is (integer) divisible by the number of cores we use

cells_per_core = int(N/size)
#print(cells_per_core)


lowarr = int(cells_per_core*rank) #actually want them to overlap so we have ghost cells though...
uparr = int(cells_per_core*(rank+1))

#now account for ghost zones in each processor
lowarr -= ghost
uparr += ghost
if rank==0:
    lowarr += ghost #just undo the ghost zone on the leftmost edge since the index can't go below 0
if rank==1:
    uparr -= ghost #undoing ghost zone on rightmost edge so we don't call an index outside of array


#divide up array
subarr = a[lowarr:uparr]
newsubarr = newa[lowarr:uparr]
print(rank, subarr)

plt.plot(subarr)
plt.savefig('advec_out/rank{}_initconds.png'.format(rank))


while count<5:
    myt+=dt
   
    
    count +=1 

    for i in range(ghost,cells_per_core-ghost):
        update=-C*dt*u*((subarr[i]-subarr[i-1])/dx) #updwind spatial deriv
        newsubarr[i]=subarr[i]+update    
    
    subarr = newsubarr+0 #assign updated values to "old" values used in next time iteration
    plt.plot(subarr)
    plt.savefig('advec_out/rank{}_time{}.png'.format(rank,count))

    print('shape of subarr is {} at rank {}'.format(np.shape(subarr), rank))
    #now deal with boundary conditions and communicating info btwn prcesoors
    #first let's deal w/ left and right boundaries
    #also, for now these boundaries only work w/ ghost cell of 1 (i.e., no higher order stencil)
    if rank==0:
        #assign leftmost edge to rightmost edge
        mysend = comm.isend(subarr[ghost],dest=size-1) #send left edge of rank 0 processor to rightmost processor
        myrec = comm.irecv(source=size-1) #recieve data from right processor
        mysend.wait()
        rightedge = myrec.wait() #now this will give us the info at right edge of box
        subarr[0]=rightedge #now just set left ghost cell to be equal to value sent from right edge
        print('rank {} recieved {} from right edge of rank {}'.format(rank, rightedge,size-1))


        print('sending {} to rank 1'.format(subarr[-(1+ghost)]))
        rightsend= comm.isend(subarr[-(1+ghost)],dest=1) #send rightmost physical cell to left edge of next processor
        rightsend.wait()
        rightrec= comm.irecv(source=1)
        rightrec_data = rightrec.wait()
        subarr[-1]=rightrec_data

        print('rank {} recieved {} from left edge of rank {}'.format(rank, rightrec_data, size-1))

    if rank==size-1:
        mysend = comm.isend(subarr[-1],dest=0)
        myrec = comm.irecv(source=0)
        mysend.wait()
        leftedge = myrec.wait()
        subarr[-1]=leftedge
        print('rank {} recieved {} from left edge of rank {}'.format(rank, leftedge, 0))

        leftsend = comm.isend(subarr[ghost],dest=0) #send leftmost physical cell to right edge of previous processor
        leftsend.wait()
        leftrec = comm.irecv(source=0)
        leftrec_data = leftrec.wait()
        subarr[0]=leftrec_data
        print ('rank {} recieved {} from right edge of rank {}'.format(rank, leftrec_data, 0))

    #for tmprank in range(1,size-1): #do the rest of the processors
        #now deal with communicating values at interface of neighboring zones
        
        #if rank==tmprank:
            
        
        
        
        #assign leftmost edge to rightmost edg

    '''
    a=newa+0 #just plus 0 to dereference a from newa
    #update ghost zones
    a[0]=a[-(1+ghost)]
    a[-1]=a[ghost]

    if count%10==0:
        plt.plot(a)
        #plt.savefig("advec_out/{}.png".format(myt))
        plt.close()
    '''