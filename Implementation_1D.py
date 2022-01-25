import math as math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation



#Constante gravitationnelle
g=1

N = 32
h = 1/(N+1)

#Temps Max
T = 2

#Nombre de sauvegardes : on sauvegarde tous les T/nSauvegarde (pour pas enregistrer à toutes les itérations)
nSauvegarde=30
j=0
tSauvegarde=[False for i in range(nSauvegarde)]



#Tableau avec 2 lignes et N colonnes : discrétisation du vecteur (h q) en N parties
#On commence avec h = q = 1 sur tout l'ensemble étudié
U = np.ones((2,N+2))

for i in range(int((N+2)/2),N+2):
        U[1,i] = .5

#Pour faire les calculs dans la boucle
Uprime = np.zeros((2,N+2))


t=0 #Temps dans la simulation
n=0 #Nombre d'itérations



X = np.linspace(0,1,N+2)





def fonctionF(U,i,tau):
    F=[]
    F.append( U[1,i+1] - U[1,i-1] )
    F.append( ( U[1,i+1]**2 ) / U[0,i+1] + (g*U[0,i+1]**2) / 2 - ( U[1,i-1]**2 ) / U[0,i-1] - (g*U[0,i-1]**2) / 2)
    F[0] = (tau/(2*h)) * F[0]
    F[1] = (tau/(2*h)) * F[1]
    return(F)
    
    
    
    

def affiche_U():
    
    X = [ h*i for i in range(N+2)]
    fig, axs = plt.subplots(1,2)
    
    axs[0].plot(X,U[0,:])
    axs[0].set_title("h")
    axs[0].set_ylim([0,1.5])
    
    axs[1].plot(X,U[1,:])
    axs[1].set_title("q")
    axs[1].set_ylim([0,1.5])
    
    plt.show()





def enregistre_U(n,t):

    t=round(t,4)
    
    fig, axs = plt.subplots(1,2)
    fig.suptitle("t="+str(t),fontsize=16)

    axs[0].plot(X,U[0,:])
    axs[0].set_title("h")
    axs[0].set_ylim([0,1.5])
    
    axs[1].plot(X,U[1,:])
    axs[1].set_title("q")
    axs[1].set_ylim([0,1.5])
    
    plt.savefig(r"C:\Users\R0MAIN\Documents\MASTER MATHS\Mémoire Recherche\Implementation1\sq_"+str(n)+".png")
    
    plt.close()


maxListe=[0 for i in range(N+2)]

enregistre_U(n,t)

while(t<T):

    print("t="+str(round(t,4)))
    
    
    #Recherche condition de stabilité
    
    for i in range(1,N+1):
        
        maxListe[i] = abs( U[1,i]/U[0,i] + math.sqrt(g * U[0,i] ) )
        
        if ( abs( U[1,i]/U[0,i] - math.sqrt(g * U[0,i] ) ) > maxListe[i] ):
            maxListe[i] = abs( U[1,i]/U[0,i] - math.sqrt(g * U[0,i] ) )
        
        if ( abs( U[1,i+1]/U[0,i+1] + math.sqrt(g * U[0,i+1] ) ) > maxListe[i] ):
            maxListe[i] = abs( U[1,i+1]/U[0,i+1] + math.sqrt(g * U[0,i+1] ) )
    
        if ( abs( U[1,i+1]/U[0,i+1] - math.sqrt(g * U[0,i+1] ) ) > maxListe[i] ):
            maxListe[i] = abs( U[1,i+1]/U[0,i+1] - math.sqrt(g * U[0,i+1] ) )

    M = max( maxListe )
  
    #Pour assurer la stabilité, tau doit être inférieur à h/(2* max(vp))
    tau = h/(2*M)
    
    
    
    for i in range(1,N+1):
        
        F = fonctionF( U , i , tau )
       
        m1 = maxListe[i]
        m2 = maxListe[i-1]
        
        Uprime[0,i] = U[0,i] - F[0] + (tau*m1/(2*h))*(U[0,i+1]-U[0,i]) - (tau*m2/(2*h))*(U[0,i]-U[0,i-1])
        Uprime[1,i] = U[1,i] - F[1] + (tau*m1/(2*h))*(U[1,i+1]-U[1,i]) - (tau*m2/(2*h))*(U[1,i]-U[1,i-1])

    
    for i in range(1,N+1):
        U[0,i] = Uprime[0,i]
        U[1,i] = Uprime[1,i]
    
        Uprime[0,i] = 0
        Uprime[1,i] = 0
    
    
    #Conditions aux bords
    U[0,0] = U[0,1]
    U[1,0] = U[1,1]
    U[0,N+1] = U[0,N]
    U[1,N+1] = U[1,N]

    t+=tau
    n+=1
    
    if(t > T/nSauvegarde*j and (j<nSauvegarde) and (not tSauvegarde[j])):
        enregistre_U(n,t)
        tSauvegarde[j]=True
        j+=1




print("Nombre d'itérations : " + str(n))














