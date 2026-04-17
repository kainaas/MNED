import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

path = os.path.join("..", "images")


def u_analitica(x,y): #The analytical solution to the PDE 
    return np.exp(-y)*np.sin(2*np.pi*x)

#h Values which the numerical calculation will be done
h_values = [1/10, 1/25, 1/50, 1/75, 1/100] 

#error vector (one entry for each h in h_value). Will be considered the infinity norm
erro = np.zeros(len(h_values))

#u_xx + 3u_yy -cu_y = 0
c = 4* np.pi**2 -3

#Values of the rectangle that defines the domain
x_0 = 0
x_n = 1
y_0 = -1
y_n = 1

#Dirichlet conditions 
u_dir = 0
u_esq = 0

# Robin condition on top: u_x + u = 0

#Neumann condition on bottom: u_y = k
k_baixo = lambda x: -np.e * np.sin(2 * np.pi * x)

#RHS function (does not depend on U)
f = lambda x,y: 0

#iterates for every h value
for l, h in enumerate(h_values):
    print(f"h = {h}")
    n_pontos_x = int((x_n - x_0)/h) - 1 #discard the point on Dirichlet BC (right and left)
    n_pontos_y = int((y_n - y_0)/h) + 1 
    n_pontos = n_pontos_x * n_pontos_y

    #Files where the plots are going to be saved
    plot_calculado = os.path.join(path, f"plot_calculado_h={h:.3f}.pdf")
    plot_analitico = os.path.join(path, f"plot_analitico_h={h:.3f}.pdf")


    #########
    #Creating the system
    #TODO: fazer isso com produto de kronecker
    #########
    A = np.diag([-8/(h**2)for i in range(n_pontos)]) #The diagonal of the matriz
    F = np.zeros(n_pontos)
    for i in range(n_pontos): #Calculates the RHS without considering BC
        x_i = x_0 + h * (i % n_pontos_x + 1)
        y_i = y_0 + h * (i // n_pontos_x) # using // for integer clean division
        F[i] = f(x_i, y_i)


    peso_cima = 3/(h**2) - c/(2*h)
    peso_baixo = 3/(h**2) + c/(2*h)

    peso_dir = 1/(h**2)
    peso_esq = 1/(h**2)

    #Iterates for each line of the matrix (we're using the natural ordering)
    for i in range(n_pontos): 
        #Defines the value of the node on top of the node with number i in the i-th line of matrix
        if i+n_pontos_x < n_pontos: #if i-th point is not a Robin BC
            A[i, i+n_pontos_x] += peso_cima
        else:
            A[i,i] -= 2*h*peso_cima
            A[i, i-n_pontos_x] += peso_cima

        #Point below
        if i-n_pontos_x >= 0: #if i-th point is not a Neumann condition
            A[i, i-n_pontos_x] += peso_baixo
        else: 
            x_i = x_0 + h * (i % n_pontos_x + 1)
            F[i] += 2*h*peso_baixo*k_baixo(x_i)
            A[i,i+n_pontos_x] += peso_baixo

        #if the point is totally to the right (the next point is a Dirichlet BC)
        if (i+1)%n_pontos_x == 0:
            F[i] -= u_dir*peso_dir
            A[i, i-1] += peso_esq
        
        #if the point is totally to the left (the point before is a Dirichlet BC)
        elif (i+1)%n_pontos_x == 1:
            F[i] -= u_esq*peso_esq
            A[i, i+1] += peso_dir
        else:
            A[i,i+1] += peso_dir
            A[i, i-1] += peso_esq



    U = np.linalg.solve(A, F)


    X = np.linspace(x_0, x_n, n_pontos_x+2)
    Y = np. linspace(y_0, y_n, n_pontos_y)
    xx, yy = np.meshgrid(X, Y)

    U_analitico = u_analitica(xx, yy)

    U_matriz = np.array([U[i:i+n_pontos_x] for i in [k*n_pontos_x for k in range(n_pontos_y)]])

    U_plot = np.zeros((n_pontos_y, n_pontos_x+2))
    U_plot[:,0] = np.array([u_esq for i in range(n_pontos_y)])
    U_plot[:,1:-1] = U_matriz
    U_plot[:,-1] = np.array([u_dir for i in range(n_pontos_y)])


    dif = np.abs(U_analitico - U_plot)
    erro[l] = np.max(dif)

    fig_calculado = plt.figure()
    ax = fig_calculado.add_subplot(111, projection='3d')
    ax.plot_surface(xx,yy,U_plot)
    fig_calculado.savefig(plot_calculado)

    fig_analitico = plt.figure()
    ax2 = fig_analitico.add_subplot(111, projection='3d')
    ax2.plot_surface(xx,yy,U_analitico)
    fig_analitico.savefig(plot_analitico)
    #plt.show()
    print("done")


#exports erro and h_values to a pickle file. Makes life easier to plot the loglog
with open('erro.pickle', 'wb') as f:
    pickle.dump(erro,f)

with open('h.pickle', 'wb') as f:
    pickle.dump(h_values, f)