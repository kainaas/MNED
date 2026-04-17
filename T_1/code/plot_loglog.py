import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

path = os.path.join("..", "images")

plot_loglog_erro = os.path.join(path, "plot_loglog_erro.pdf")

with open('erro.pickle', 'rb') as f:
    erro = pickle.load(f)

with open('h.pickle', 'rb') as f:
    h_values = pickle.load(f)

x= np.array([5e-2, 8e-2])
y = x**2



fig_loglog = plt.figure()
ax = fig_loglog.add_subplot(111)
ax.loglog(np.array(h_values), erro, color='red', marker='v', label='erro calculado')
ax.loglog(x,y, color='black', ls='--', label=r'$h^2$')
ax.set_title("Plot log-log do erro na norma infinito")
ax.set_xlabel("h")
ax.set_ylabel("Erro na norma infinito")
ax.legend()

fig_loglog.savefig(plot_loglog_erro)
plt.show()
