import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

save_figs = True
refresh_pickles = False

path = os.path.join("..", "images")


def u_analitica(x,y): #The analytical solution to the PDE 
    return np.exp(-y)*np.sin(2*np.pi*x)

#h Values which the numerical calculation will be done
h_values = [1/10, 1/25, 1/50, 1/75, 1/100] 

#error vector (one entry for each h in h_value). Will be considered the infinity norm
erro = np.zeros(len(h_values))

#error vector in 1-norm for grid functions
erro_n1 = np.zeros(len(h_values))

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
    n = int((x_n - x_0)/h) - 1 #discard the point on Dirichlet BC (right and left)
    m = int((y_n - y_0)/h) - 1

    N = n*(m+2)

    

    #Files where the plots are going to be saved
    calculado2D = os.path.join(path, f"calculado2D/h={h:.3f}.pdf")
    calculado3D = os.path.join(path, f"calculado3D/h={h:.3f}.pdf")
    analitico2D = os.path.join(path, f"analitico2D/h={h:.3f}.pdf")
    analitico3D = os.path.join(path, f"analitico3D/h={h:.3f}.pdf")
    plot_dif = os.path.join(path, f"dif/h={h:.3f}.pdf")


    #===================
    #Creating the system
    #===================
    F = np.zeros(N).T

    M = np.zeros((N, N))
    In = np.eye(n)
    Im2 = np.eye(m+2)

    alpha = 3/h**2 + (4*np.pi**2 - 3)/(2*h)
    beta  = 3/h**2 - (4*np.pi**2 - 3)/(2*h)

    F[0:n] = np.array([-(6/h + 4*np.pi**2 - 3)*np.e*np.sin(2*np.pi*i*h) for i in range(1,n+1)])

    # D- + D0 + D+
    Tn = np.zeros((n, n))
    i = np.arange(0, n)
    Tn[i[1:], i[1:]-1] = 1/h**2 #D-
    Tn[i, i] = -4/h**2 #D0/2
    Tn[i[:-1], i[:-1]+1] = 1/h**2 #D+

    # A-n + B+n
    Tm2 = np.zeros((m+2, m+2))
    i = np.arange(0, m+2)
    Tm2[i[1:], i[1:]-1] = alpha #A-n
    Tm2[i, i] = -4/h**2 #D0/2
    Tm2[i[:-1], i[:-1]+1] = beta #B+n

    M = np.kron(Im2, Tn) + np.kron(Tm2, In)

    #P-n
    i = np.arange(0, N)
    j = i[n*(m+1):]
    M[j, j-n] += beta

    #G0
    M[j, j] += -2*h*beta

    #Q+n
    j = i[:n]
    M[j, j+n] += alpha


    U = np.linalg.solve(M, F)


    X = np.linspace(x_0, x_n, n+2)
    Y = np. linspace(y_0, y_n, m+2)
    xx, yy = np.meshgrid(X, Y)

    U_analitico = u_analitica(xx, yy)

    U_matriz = np.array([U[i:i+n] for i in [k*n for k in range(m+2)]])

    U_plot = np.zeros((m+2, n+2))
    U_plot[:,0] = np.array([u_esq for i in range(m+2)])
    U_plot[:,1:-1] = U_matriz
    U_plot[:,-1] = np.array([u_dir for i in range(m+2)])


    dif = np.abs(U_analitico - U_plot)
    erro[l] = np.max(dif)
    erro_n1[l] = dif.sum().sum()*h**2

    print(f"Erro na norma infinito: {float(erro[l])}")
    print(f"Erro na norma 1: {float(erro_n1[l])}")

    if save_figs:
        fig_calculado = plt.figure()
        ax = fig_calculado.add_subplot(111)
        plt.colorbar(ax.contourf(xx,yy,U_plot,levels=100,cmap='coolwarm'),label=r'$U_{i,j}$')
        fig_calculado.savefig(calculado2D)

        fig_analitico = plt.figure()
        ax2 = fig_analitico.add_subplot(111)
        plt.colorbar(ax2.contourf(xx,yy,U_analitico,levels=100,cmap='coolwarm'),label=r'$u(x_i,y_j)$')
        fig_analitico.savefig(analitico2D)

        fig_calculado = plt.figure()
        ax = fig_calculado.add_subplot(111, projection='3d')
        plt.colorbar(ax.plot_surface(xx,yy,U_plot,cmap='coolwarm'),label=r'$U_{i,j}$')
        fig_calculado.savefig(calculado3D)

        fig_analitico = plt.figure()
        ax2 = fig_analitico.add_subplot(111, projection='3d')
        plt.colorbar(ax2.plot_surface(xx,yy,U_analitico,cmap='coolwarm'),label=r'$u(x_i,y_j)$')
        fig_analitico.savefig(analitico3D)

        fig_dif = plt.figure()
        ax3 = fig_dif.add_subplot(111)
        c = ax3.pcolormesh(xx, yy, dif, cmap='inferno', shading='auto')  #creates colormap
        fig_dif.colorbar(c, ax=ax3)    # pass the mappable
        fig_dif.savefig(plot_dif)

        plt.close('all')
        #plt.show()

    print("done")


if refresh_pickles:
    #exports erro and h_values to a pickle file. Makes life easier to plot the loglog
    with open('erro.pickle', 'wb') as f:
        pickle.dump(erro,f)

    with open('h.pickle', 'wb') as f:
        pickle.dump(h_values, f)