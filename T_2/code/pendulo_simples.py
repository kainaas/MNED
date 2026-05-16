import numpy as np
import metodos
from RefSolution import RefSolution
import matplotlib.pyplot as plt
import pickle

#q'' = -sin(q)
#u' = f(u,t)
#Logo, fazemos u = (q, q')
#u0 = (posição inicial, velocidade inicial)
#f(u,t)=(q', -sin(q))
#      =(u2, -sin(u1))

h = 1
pos_ini = np.pi / 4
vel_ini = 0

t0 = 0
t_final = 12

u0 = np.array([pos_ini, vel_ini])
f = lambda u,t: np.array([u[1], -np.sin(u[0])])

Jf = lambda u,t: np.array([[0, 1], [-np.cos(u[0]), 0]])

result_explicito = metodos.euler_explicito(t0, t_final, h, u0, f)
result_implicito = metodos.euler_implicito(t0, t_final, h, u0, f, Jf)
result_RK4 = metodos.RK4(t0, t_final, h, u0, f)
result_ref = RefSolution(t0, t_final, h, pos_ini)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(result_ref.t, result_ref.u[0,:], '.-', label="Solução de referência")
ax.plot(result_RK4.t, result_RK4.u[0,:], '.--', label="RK4")
ax.plot(result_explicito.t, result_explicito.u[0,:], '.--', label="Euler explícito")
ax.plot(result_implicito.t, result_implicito.u[0,:], '.--', label="Euler implícito")
ax.grid()
ax.legend()
plt.show()