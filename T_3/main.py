import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import lil_matrix

# Consts
eps = 0.1
P = [100, 500, 1000]
N = 430 # Number of points in x ('till 140 takes less than 3s)
Dx = 15/(N-1)
Nx = int(np.floor(15/Dx)) + 1

animate = True
show = True
save = False
feedback = True

# Cases (I), (II) e (III) for first discretization
for k,Pe in enumerate(P):
    Dt = np.min([ Dx**2 / (2/Pe + Dx), Dx ])*(1-eps)
    Nt = int(np.floor(12/Dt)) + 1
    #print(f'Dx={Dx:.4e}, Dt={Dt:.4e}')

    xj = np.arange(0,Nx,1)*Dx
    tn = np.arange(0,Nt,1)*Dt
    U0 = np.exp(-20*(xj-2)**2) + np.exp(-(xj-5)**2)

    # Initial condition
    if show or save:
        if k == 0:
            fig = plt.figure(figsize=(9,5))
            initial = fig.add_subplot(111)

            x = np.linspace(0,15,1000)
            u0 = np.exp(-20*(x-2)**2) + np.exp(-(x-5)**2)
            initial.plot(x, u0, '-k', label=r'$u(x,0)$')
            initial.plot(xj, U0, '.r', zorder=2, label=r'$U(x_j,0)$')

            initial.set_xticks(np.arange(0,16))
            initial.set_xlabel('x')
            initial.set_ylabel('u')
            initial.grid(True)
            initial.legend()

            plt.tight_layout()
            if show: plt.show()
            if save: plt.savefig('graphs/initial_condition.pdf', bbox_inches='tight', format='pdf')

    # First discretization
    m = Nx
    A = lil_matrix((m,m))

    a = Dt / (Pe*Dx**2)
    b = Dt / (2*Dx)

    # Inner lines
    for j in range(1, m-1):
        A[j, j-1] = a+b
        A[j, j] = 1-2*a
        A[j, j+1] = a-b

    # Periodic boundary conditions
    A[0, 0] = 1-2*a
    A[0, 1] = a-b
    A[0, m-1] = a+b

    A[m-1, 0] = a-b
    A[m-1, m-2] = a+b
    A[m-1, m-1] = 1-2*a

    A = A.tocsr()

    # Plots
    X, T = np.meshgrid(xj, tn)

    U = np.zeros((Nt, Nx))
    U[0, :] = U0

    if feedback: print(f'Calculating Pe={Pe}...', end='', flush=True)
    for n in range(1, Nt):
        U[n, :] = A.dot(U[n-1, :])
    if feedback and (show or save): print(' Ploting...', end='', flush=True)

    if show or save:
        fig = plt.figure(figsize=(9,5))
        final = fig.add_subplot(111, projection='3d')

        surf = final.plot_surface(X,T,U,cmap='plasma')
        fig.colorbar(surf, ax=final, shrink=0.8, pad=0.1)

        final.set_xlabel('x')
        final.set_ylabel('t')
        final.set_zlabel(r'$U$')
        final.grid(True)

        if show: final.set_title(r'$\mathbb{P}e = '+f'{Pe}'+r'$')

        plt.tight_layout()
        if show: plt.show()
        # For Pe << 1 and Pe = 1, pdf gets too big, so used png
        if save: plt.savefig(f'graphs/surface_central_Pe={Pe}.png', bbox_inches='tight', format='png', dpi=300)
    if feedback: print(' Done!')

    if animate:
        plt.close('all')

        fig = plt.figure(figsize=(10,5))
        evol = fig.add_subplot(111)

        n = U.shape[0]
        it = 100
        for k in range(n):
            if not plt.fignum_exists(fig.number):
                plt.close('all')
                break

            if k%int(n/it) == 0 or k+1 == n:
                evol.cla()
                evol.plot(X[k,:],U[k,:])
                evol.set_title(r'$\mathbb{P}e = '+f'{Pe}'+r',\ t='+f'{k}/{n-1}'+'$')
                plt.pause(1) if k == 0 else plt.pause(1/it)
    
    plt.close('all')