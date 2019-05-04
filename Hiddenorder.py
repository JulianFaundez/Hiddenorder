import numpy as np
import matplotlib.pyplot as plt
import math

import time 



PASO=0.1
TD=0.3 #termo de hopping entre bandas

J=0.128 #termino de exchange entre bandas
KB=8.6e2#constante de boltzman
e1=0
e3=0
e2=0
NX=55
NY=55
NZ=55
for t in np.arange(0.01,12,1):
    tol=1.0
    q0=0.1
    while tol>0.01:
          e3=0 
          for kx in np.arange(-math.pi, math.pi, PASO):
              e2=0
              for ky in np.arange(-math.pi, math.pi, PASO):
                  e1=0     
                  for kz in np.arange(-math.pi, math.pi, PASO): 

                      
                      Eak= -2*TD*(np.cos(kx)+np.cos(ky)+np.cos(kz) )
            
                      EakQ= -2*TD*(np.cos(kx+math.pi)+np.cos(ky+math.pi)+np.cos(kz+math.pi) )
         
                      Ebk= -2*TD*(np.cos(kx)+np.cos(ky)+np.cos(kz) )

                      EbkQ= -2*TD*(np.cos(kx+math.pi)+np.cos(ky+math.pi)+np.cos(kz+math.pi) )
 
                      eaux= J*J*q0*q0

                      #Emas= ( (EakQ+Ebk)/2)+ math.sqrt( pow( ((EakQ-Ebk)/2),2) + eaux)
 
                      #Emenos= ( (EakQ+Ebk)/2)- math.sqrt( pow( ((EakQ-Ebk)/2),2) + eaux)
                      p=np.roots([1, -((EakQ)+(Ebk)),  ((EakQ)*(Ebk)-eaux)])  
                      EmasN=p[0]
                      EmenosN=p[1]

                      EmasmenosN=(EmasN-EmenosN)
 
                      #Emasmenos= (Emas-Emenos)
                       
                      beta=(1/(KB*t))

                      betamasN= (beta*EmasN)

                      betamenosN=(beta*EmenosN)
                     
                      fermimas=(1/(1+np.exp(betamasN)))
  
                      fermimenos=(1/(1+np.exp(betamenosN)))

                      e1=e1+((fermimenos-fermimas)/(EmasmenosN)) 
           #fin for kz
                  e2=e2+e1 
        #finforky    
              e3=e3+e2

          e=(J*q0*e3)
          tol=abs(e-q0)
          q0=e
          
    z=e/(NX*NY*NZ)
    print("%.2f" % t,"%0.4f" %  z)
    plt.plot(z, t,'*')
#plt.xlabel('time [s]')
#plt.ylabel('signal')
plt.show()

