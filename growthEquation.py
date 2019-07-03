import numpy as np
import scipy as sp
import scipy.linalg as la
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from scipy import sparse
import time
counter = 0 #a constant to print out the plots
v = 1 #difusion constant
scale = 100
#this part aims to creat a tridiagonal matrix
coeffMat = np.zeros((scale,scale)) #making the coefficent matrix
diameter0 = np.ones(scale)*(1+2*v) #diagonal matrix
diameter1 = np.ones(scale)*(-1*v) #diagonal matrix
diagonalMat0 = np.diag(diameter0).reshape((scale,scale))
diagonalMat = np.diag(diameter1).reshape((scale,scale))
diagonalMat1 = np.roll(diagonalMat,1)
diagonalMat2 = np.roll(diagonalMat, -1)
coeffMat = diagonalMat0 + diagonalMat1 + diagonalMat2
coeffMat[0][0] = coeffMat[scale-1][scale-1] = 1+2*v
coeffMat[scale-1][0] = coeffMat[0][scale-1] = -v
coeffMat
sparsCoeffMat = sparse.coo_matrix(coeffMat)
##############################################
#initial matrix valuse of the line (h)
knownMat = np.zeros(scale) 
unknownMat = np.zeros(scale)
noise = np.random.normal(0,0.9,scale) #adding noise due to thermal effects
meanValue = [] #to store the mean value of h each itteration
itterations = 20000 #number of times the linear algebra system would solve
w = []
for i in range(0,itterations):
    unknownMat = spla.spsolve(sparsCoeffMat, knownMat) #solving the system by sparse methode to save times
    #unknownMat[scale-1] = unknownMat[0]
    noise = np.random.normal(0,0.9,scale) #adding noise due to thermal effects
    unknownMat = unknownMat + noise
    knownMat[1:scale-1] = unknownMat[1:scale-1]
    meanValue.append(np.mean(unknownMat))
    temporalConst = 0.0
    if ((i*100)%itterations==0):
        print('{} %'.format(i/(itterations/100)))
    if (i%20==0):
        for j in range(0,scale-1):
            temporalConst += (abs(unknownMat[j] - np.mean(unknownMat)))**2

        w.append(np.sqrt(temporalConst*1/scale)) #store the w value

    if (i%200==0): #to reduce the number of plots for every 100 loops
        counter += 1
        #to create an animation
        #plt.cla()
        #plt.figure(figsize=(10,15))
        plt.title('1D Growth Equation')
        plt.xlabel('Scale')
        plt.ylabel('Growth')
        plt.ylim(-14,14)
        plt.plot(unknownMat)
        #plt.legend(loc='upper right')
        #plt.savefig('/gif/{}.png'.format(counter)) #in case of creating gif
        plt.pause(0.01)
        time.sleep(0.01)


#plt.show()
#len(meanValue)
plt.cla()
plt.plot(range(0,itterations),meanValue)
plt.legend(loc='upper right')
plt.ylabel('h')
plt.xlabel('time')
plt.savefig('averageOfH.png')
############################
#calculate w^2
plt.cla()
plt.legend(loc='upper right')
plt.plot(w,label='w^2 during the time');
plt.ylabel('w')
plt.xlabel('time')
plt.savefig('w.png')
#w1 = np.sqrt(w2)
plt.cla()
plt.legend(loc='upper right')
plt.plot(np.log(w),label='w^2 during the time');
plt.ylabel('log(w)')
plt.xlabel('time')
plt.savefig('Logw.png')

print('Done!')
