import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import lil_matrix

# Consts
eps = 0.1
P = [0.1, 1, 20]
N = 140 # Number of points in x (temporary for tests, 'till 140 takes less than 3s)
Dx = 15/(N-1)
Nx = int(np.floor(15/Dx)) + 1

# Cases (I), (II) e (III) for first discretization
for k,Pe in enumerate(P):
    Dt = np.min([ Dx**2 / (2/Pe + Dx), Dx ])*(1-eps)
    Nt = int(np.floor(12/Dt)) + 1

    xj = np.arange(0,Nx,1)*Dx
    tn = np.arange(0,Nt,1)*Dt
    U0 = np.exp(-20*(xj-2)**2) + np.exp(-(xj-5)**2)

    # Initial condition
    if k == 0:
        fig = plt.figure(figsize=(10,5))
        initial = fig.add_subplot(111)

        x = np.linspace(0,15,1000)
        u0 = np.exp(-20*(x-2)**2) + np.exp(-(x-5)**2)
        initial.plot(x, u0, '-k')
        initial.plot(xj, U0, '.r', zorder=2)

        initial.set_xticks(np.arange(0,16))

        plt.tight_layout()
        plt.show()

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

    for n in range(1, Nt):
        U[n, :] = A.dot(U[n-1, :])

    fig = plt.figure(figsize=(10,5))
    final = fig.add_subplot(111, projection='3d')
    final.plot_surface(X,T,U,cmap='plasma')

    final.set_xlabel('x')
    final.set_ylabel('t')
    final.set_zlabel(r'$U$')

    final.set_title(r'$\mathbb{P}e = '+f'{Pe}'+r'$')

    plt.tight_layout()
    plt.show()