import numpy as np
import matplotlib.pyplot as plt

m = 20 #número máximo de iterações de Newton

g = lambda U0, U1, U2, h, mi: (1/h**2)*(U0 - 2*U1  + U2 -(mi*h/2)*(U1**2*U2 - U1**2*U0 - U2 + U0)) + U1

g_f = lambda U0, U1, h, mi, sigma: (1/(h**2))*(2*U0 - 2*U1 + 2*h*sigma) + (sigma*mi)*(1- U1**2) + U1

h_vals = [1/5]

Erro = np.zeros(3)

a = 0
b = 2

alpha = 0
sigma = 1
mi = 1/2

for index, h in enumerate(h_vals):
    n=int((b-a)/h + 1) #número de pontos
    x = np.array(np.linspace(a, b, num=n)) 
    U = x*np.ones(n)

    U[0] = alpha
    

    
    for j in range(m):

        J = np.zeros((n,n))

        J[0,0] = 1
        J[n-1,n-1] = -2/h**2 - (2*sigma*mi) * U[n-1] + 1
        J[n-1, n-2] = 2/(h**2)
        
        G = np.zeros(n)
        G[0] =  U[0]- alpha

        for i in range(1,n-1):

            G[i] = g(U[i-1],U[i],U[i+1], h, mi)

            J[i,i]   = 1 -2/h**2 + (mi/h)*U[i]*(U[i-1] - U[i+1])
            J[i,i+1] = 1/h**2 + mi/(2*h)*(U[i]**2 - 1)
            J[i,i-1] = 1/h**2 + mi/(2*h)*(-U[i]**2 + 1)
        
        G[n-1] = g_f(U[n-2], U[n-1], h, mi, sigma)


        delta = np.linalg.solve(J, -G)
        U = U + delta
        erro = np.max(delta)
        print(np.linalg.norm(delta, np.inf))
        if np.linalg.norm(delta, np.inf) < 1e-5:
            print(index)
            break



# h_vals.reverse()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x, U, color = 'red', s =  5)
# ax.plot(x, U_exato)
# ax.set_title("Comparação da solução exata e da calculada")
# ax2 = fig.add_subplot(122)
# ax2.loglog(h_vals, Erro, marker = 'o',)
# ax2.set_title("Erro(h)")
plt.show()