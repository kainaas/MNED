import numpy as np
import matplotlib.pyplot as plt

m = 2 #número de iterações de Newton

f = lambda U0, U1, U2, h, x: 4*U0 + U0**2 +4*U2 + U2**2 -8*U1 + 4*(h**2) * U1 - 2*U0*U2 - 4*(h**2)*np.log(x)

a=1 #início do intervalo
b=2 #fim do intevalo
h_vals = [1/2, 1/4, 1/8, 1/16, 1/32, 1/64]

Erro = np.zeros(6)

for index, h in enumerate(h_vals):
    n=int((b-a)/h + 1) #número de pontos
    x = np.array(np.linspace(1.0, 2.0, num=n)) 
    U = np.zeros(n)

    U[0] = 0
    U[n-1] = np.log(2)


    J = np.zeros((n,n))
    J[0,0] = 1
    J[n-1,n-1] = 1

    D = np.diag((-8+4*h**2) *np.ones(n))
    D[0,0] = 0
    D[n-1,n-1] = 0
    F = np.zeros(n)
    F[0] = U[0]
    F[n-1] = U[n-1]
    for j in range(m):

        J = np.zeros((n,n))
        F = np.zeros(n)

        J[0,0] = 1
        J[n-1,n-1] = 1
        F[0] = U[0]
        F[n-1] = U[n-1] - np.log(2)

        for i in range(1,n-1):

            F[i] = f(U[i-1],U[i],U[i+1],h,x[i])

            J[i,i]   = -8 + 4*h**2
            J[i,i+1] = 4 + 2*U[i+1] - 2*U[i-1]
            J[i,i-1] = 4 + 2*U[i-1] - 2*U[i+1]

        U_sol = np.linalg.solve(J, J@U - F)
        U = U_sol

        #A solução real é u(x)= ln(x)

        U_exato = np.log(x)
        E = U - U_exato
        Erro[index] = max(np.abs(E))


h_vals.reverse()

fig = plt.figure()
ax = fig.add_subplot(121)
ax.scatter(x, U, color = 'red', s =  5)
ax.plot(x, U_exato)
ax.set_title("Comparação da solução exata e da calculada")
ax2 = fig.add_subplot(122)
ax2.loglog(h_vals, Erro, marker = 'o',)
ax2.set_title("Erro(h)")
plt.show()