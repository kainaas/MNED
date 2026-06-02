#Implementa a equação de Helmholtz num domínio bidimensional retangular,
#com condições de Robin nas bordas esquerda e direitas e cond. de Dirichlet em cima e em baixo


import numpy as np
import matplotlib.pyplot as plt
import math


h = 1/10

#u_xx + u_yy +ku = 0
k = -10


x_0 = -1
x_n = 1
y_0 = -1
y_n = 1

#Dirichlet
u_cima = 0
u_baixo = 1

#Condições de Robin da forma u_x = alpha*u
alpha_dir = -1/2
alpha_esq = 1/2

#Função na direita
f = lambda x,y: 0

n_pontos_x = (int) ((x_n - x_0)/h + 1)
n_pontos_y = (int) ((y_n - y_0)/h - 1) #descarta os pontos nas condições de Dirichlet
n_pontos = n_pontos_x * n_pontos_y


A = np.diag([-4 + k*h**2 for i in range(n_pontos)])
F = np.zeros(n_pontos)
for i in range(n_pontos):
    F[i] = f(x_0+h*(i%n_pontos_x), 
             y_0+h*(math.floor(i/n_pontos_x)+1))


for i in range(n_pontos):
    #Toma o ponto acima do ponto que estamos e coloca-o na linha correspondente da matriz
    if i+n_pontos_x < n_pontos:
        A[i, i+n_pontos_x] = 1
    else:
        F[i] += -u_cima/h**2

    #Ponto abaixo
    if i-n_pontos_x >= 0:
        A[i, i-n_pontos_x] = 1
    else:
        F[i] += -u_baixo/h**2

    #Se é um ponto totalmente a direita
    if (i+1)%n_pontos_x == 0:
        A[i, i-1] = 2
        A[i,i] += 2*h*alpha_dir
    
    #Se é um ponto totalemnte à esquerda
    elif (i+1)%n_pontos_x == 1:
        A[i, i+1] = 2
        A[i, i] += -2*h*alpha_esq
    else:
        A[i,i+1] = 1
        A[i, i-1] = 1


A *= 1/h**2

U = np.linalg.solve(A, F)


X = np.linspace(x_0, x_n, n_pontos_x)
Y = np. linspace(y_0, y_n, n_pontos_y+2)
xx, yy = np.meshgrid(X, Y)

U_matriz = np.array([U[i:i+n_pontos_x] for i in [k*n_pontos_x for k in range(n_pontos_y)]])

U_plot = np.zeros((n_pontos_y+2, n_pontos_x))
U_plot[0,:] = np.array([u_baixo for i in range(n_pontos_x)])
U_plot[1:-1,:] = U_matriz
U_plot[-1,:] = np.array([u_cima for i in range(n_pontos_x)])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xx,yy,U_plot)
plt.show()