#Implementa a equação de Helmholtz num domínio bidimensional retangular,
#com condições de Robin nas bordas esquerda e direitas e cond. de Dirichlet em cima e em baixo


import numpy as np
import matplotlib.pyplot as plt
import math


h = 1/10

def u_analitica(x,y):
    return np.exp(-y)*np.sin(2*np.pi*x)

#u_xx + 3u_yy -cu_y = 0
c = 4* np.pi**2 -3


x_0 = 0
x_n = 1
y_0 = -1
y_n = 1

#Dirichlet
u_dir = 0
u_esq = 0

#Condição de Robin no topo da forma u_x + u = 0

#Condição de Neumann em baixo da forma u_y = k
k_baixo = lambda x: -np.e * np.sin(2 * np.pi * x)

#Função na direita
f = lambda x,y: 0

n_pontos_x = int((x_n - x_0)/h) - 1 #descarta os pontos nas condições de Dirichlet (direita e esquerda)
n_pontos_y = int((y_n - y_0)/h) + 1 
n_pontos = n_pontos_x * n_pontos_y


A = np.diag([-8/(h**2)for i in range(n_pontos)])
F = np.zeros(n_pontos)
for i in range(n_pontos):
    x_i = x_0 + h * (i % n_pontos_x + 1)
    y_i = y_0 + h * (i // n_pontos_x) # Usando // para divisão inteira limpa
    F[i] = f(x_i, y_i)


peso_cima = 3/(h**2) - c/(2*h)
peso_baixo = 3/(h**2) + c/(2*h)

peso_dir = 1/(h**2)
peso_esq = 1/(h**2)

for i in range(n_pontos): #TODO: fazer isso com produto de kronecker
    #Toma o ponto acima do ponto que estamos e coloca-o na linha correspondente da matriz
    if i+n_pontos_x < n_pontos:
        A[i, i+n_pontos_x] = peso_cima
    else:
        A[i,i] -= 2*h*peso_cima
        A[i, i-n_pontos_x] += peso_cima

    #Ponto abaixo
    if i-n_pontos_x >= 0:
        A[i, i-n_pontos_x] = peso_baixo
    else:
        x_i = x_0 + h * (i % n_pontos_x + 1)
        F[i] += 2*h*peso_baixo*k_baixo(x_i)
        A[i,i+n_pontos_x] += peso_baixo

    #Se é um ponto totalmente a direita
    if (i+1)%n_pontos_x == 0:
        F[i] -= u_dir*peso_dir
        A[i, i-1] = peso_esq
    
    #Se é um ponto totalemnte à esquerda
    elif (i+1)%n_pontos_x == 1:
        F[i] -= u_esq*peso_esq
        A[i, i+1] = peso_dir
    else:
        A[i,i+1] = peso_dir
        A[i, i-1] = peso_esq



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

fig = plt.figure()
ax = fig.add_subplot(221, projection='3d')
ax2 = fig.add_subplot(222, projection='3d')
ax3 = fig.add_subplot(223, projection='3d')
ax.plot_surface(xx,yy,U_plot)
ax2.plot_surface(xx,yy,U_analitico)
ax3.plot_surface(xx,yy,dif)
plt.show()