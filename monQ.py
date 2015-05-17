import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
def K(k1,k2,M):
    return np.linspace(k1,k2,M)

def instants(tf,N):
    return np.linspace(0,tf,N)

g=9.81

def onde(k,h0,x0,y0,x,y,instants):
    d=np.sqrt((x-x0)**2+(y-y0)**2)
    w=np.sqrt(k*g*np.tanh(k*h0))
    c=w/k ; T=[0]*len(instants)
    for j in range(len(instants)):
        if d>c*instants[j]:
            T[j]=0
        else:
            T[j]=(1/d)*np.cos(k*d-w*instants[j])
    return T

def elevationvecteuronde(lP,LP,vitesse,N,k,h0,x,y,instants):
    E=[0]*len(instants)
    for i in range(len(instants)):
            X=np.linspace(vitesse*instants[i],vitesse*instants[i]+LP,N)
            Y=np.linspace(-lP/2,lP/2,N)
            S=onde(k,h0,X[0],Y[0],x,y,instants)[i]
            for l in range(1,N):
                for j in range(1,N):
                    S+=onde(k,h0,X[l],Y[j],x,y,instants)[i]
            E[i]=S
    return E

def elevationsurface(x,y,K,lP,LP,vitesse,instants,N,h0):
    L=[0]*len(instants)
    for i in range(len(instants)):
            S=elevationvecteuronde(lP,LP,vitesse,N,K[0],h0,x,y,instants)[i]
            for j in range(1,len(K)):
                S+=elevationvecteuronde(lP,LP,vitesse,N,K[j],h0,x,y,instants)[i]
            L[i]=S
    return L

def representationelevation(K,lP,LP,vitesse,N,instants,lS,LS,h0):
    Y=np.linspace(-lS/2,lS/2,N)
    X=np.linspace(0,LS,N)
    Z=np.zeros((N,N),float)
    for i in range(1+int(((N-1)/LS)*vitesse*instants[0])):
        for l in range(1+int(((N-1)/(2*lS))*(lS-lP))):
            Z[i,l]=elevationsurface(X[i],Y[l],K,lP,LP,vitesse,instants,N,h0)[0]
    for i in range(int(((N-1)/LS)*(LP+vitesse*instants[0])),N):
        for l in range(1+int(((N-1)/(2*lS))*(lS-lP))):
            Z[i,l]=elevationsurface(X[i],Y[l],K,lP,LP,vitesse,instants,N,h0)[0]
        for l in range(int(((N-1)/(2*lS))*(lS+lP)),N):
            Z[i,l]=Z[i,-l+N-1]
    for i in range(1+int(((N-1)/LS)*vitesse*instants[0])):
        for l in range(int(((N-1)/(2*lS))*(lS+lP)),N):
            Z[i,l]=Z[i,-l+N-1]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X,Y = np.meshgrid(X,Y)
    surf = ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.xlabel('axe des x')
    plt.ylabel('axe des y')
    plt.title("hauteur de la surface de l'eau a un instant donne") 
    plt.axis('equal')
    plt.savefig('graphe simulation')
    plt.show()
    return "c'est fini"



        

