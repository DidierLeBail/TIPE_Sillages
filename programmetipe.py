import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
def K(k1,k2,M):
    return np.linspace(k1,k2,M)                                #choix des vecteurs d'onde qui interféront. J'ai supposé que le mobile émettait à chaque instant 
                                                               #des ondes concentriques dont la longueur d'onde varirait de 0 à l'infini, et que ces vagues étaient
def instants(tf,N):
    return np.linspace(0,tf,N)                                 #générées à amplitudes égales. Les points d'émission des vagues sont tous les points en contact 

g=9.81                                                         #avec le fluide, mais finalement, parce que le temps de calcul était trop long, j'ai considéré
                                                               #seulement les vagues émises par les points à l'intersection des bords du mobile et de                                                                                                                                
def onde(k,h0,x0,y0,x,y,instant):                              #la surface de l'eau.
    d=np.sqrt((x-x0)**2+(y-y0)**2)          
    w=np.sqrt(k*g*np.tanh(k*h0))                               
    c=w/k 
    if d>c*instant:                                            
    	return 0
    else:
    	return (1/np.sqrt(d))*np.cos(k*d-w*instant)            
def elevationvecteuronde(lP,LP,vitesse,N,k,h0,x,y,instants,i):
	if i==0:
		return 0
	else:                                                      #lP désigne la largeur de la perturbation et LP sa longueur (mon mobile est rectangulaire    
		Y=np.linspace(-lP/2,lP/2,N)                            #et se déplace selon les x croissants). 
		X=np.linspace(0,LP,N)
		E=onde(k,h0,X[0],Y[0],x,y,instants[i])
		for l in range(1,N):                                                              #ici on calcule l'élévation de la surface de l'eau en (x,y) due à un vecteur
			E+=onde(k,h0,X[0],Y[l],x,y,instants[i])+onde(k,h0,X[-1],Y[l],x,y,instants[i]) #d'onde émis par tous les points du périmètre mouillé du mobile.
        for m in range(1,N):                                                              #(pour aller plus vite, les calculs sont faits pour deux côtés opposés
        	E+=onde(k,h0,X[m],Y[0],x,y,instants[i])+onde(k,h0,X[m],Y[-1],x,y,instants[i]) #en même temps)
    	for j in range(1,i):                                                          
            X=np.linspace(vitesse*instants[j],vitesse*instants[j]+LP,N)
            S=onde(k,h0,X[0],Y[0],x,y,instants[i-j])
            for l in range(1,N):
            	S+=onde(k,h0,X[0],Y[l],x,y,instants[i-j])+onde(k,h0,X[-1],Y[l],x,y,instants[i-j])
            for m in range(1,N):
            	S+=onde(k,h0,X[m],Y[0],x,y,instants[i-j])+onde(k,h0,X[m],Y[-1],x,y,instants[i-j])
            E+=S
    	return E


def elevationsurface(x,y,K,lP,LP,vitesse,instants,i,N,h0):                               #là on fait varier k 
	S=elevationvecteuronde(lP,LP,vitesse,N,K[0],h0,x,y,instants,i)
	for j in range(1,len(K)):
		S+=elevationvecteuronde(lP,LP,vitesse,N,K[j],h0,x,y,instants,i)
	return S

def representationelevation(K,lP,LP,vitesse,N,instants,i,lS,LS,h0):                      #lS : largeur de la piscine (LS:longueur)
    Y=np.linspace(-lS/2,lS/2,N) 
    X=np.linspace(0,LS,N)                                                                #et là on fait varier (x,y)
    Z=np.zeros((N,N),float)

    for j in range(1+int(((N-1)/float(LS))*vitesse*instants[i])):                        #les calculs ne sont pas faits pour les points du mobile
        for l in range(1+int(((N-1)/float(2*lS))*(lS-lP))):                              
            Z[i,l]=elevationsurface(X[i],Y[l],K,lP,LP,vitesse,instants,i,N,h0)
            

    for m in range(int(((N-1)/float(LS))*(LP+vitesse*instants[i])),N):
        for k in range(1+int(((N-1)/float(2*lS))*(lS-lP))):
            Z[m,k]=elevationsurface(X[m],Y[k],K,lP,LP,vitesse,instants,i,N,h0)
            
        for p in range(int(((N-1)/float(2*lS))*(lS+lP)),N):                              #j'utilise la symétrie du problème par rapport à l'axe Ox 
            Z[m,p]=Z[m,-p+N-1]
            

    for o in range(1+int(((N-1)/float(LS))*vitesse*instants[i])):
        for n in range(int(((N-1)/float(2*lS))*(lS+lP)),N):
            Z[o,n]=Z[o,-n+N-1]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X,Y = np.meshgrid(X,Y)
    surf = ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False) #ici j'ai recopié les mêmes instructions trouvées 
    fig.colorbar(surf, shrink=0.5, aspect=5)                                                         #sur une aide pour dessiner des surfaces.
    plt.xlabel('axe des x')
    plt.ylabel('axe des y')
    plt.title("hauteur de la surface de l'eau a un instant donne") 
    plt.axis('equal')
    plt.savefig('graphe simulation')
    plt.show()
    return "c'est fini"                                                           #là c'est juste pour mettre un return. Mais ça ne marche pas parce que le 
                                                                                  #programme ne s'arrête que quand on ferme la fenêtre graphique.

representationelevation(K(9.81,2*9.81,10),0.08,0.15,1,30,instants(10,10),3,6,12,2)

        

