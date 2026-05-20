import pickle
import numpy as np
import matplotlib.pyplot as plt
from os import path
from metodos import Resultado

results_dir = "resultados_numericos"
images_dir = path.join("..", "images")
show = True

with open(path.join(results_dir, "e_explicito"), 'rb') as f:
    result_explicito = pickle.load(f)
with open(path.join(results_dir, "e_implicito"), 'rb') as f:
    result_implicito = pickle.load(f)
with open(path.join(results_dir, "RK4"), 'rb') as f:
    result_RK4 = pickle.load(f)
with open(path.join(results_dir, "ref"), 'rb') as f:
    result_ref = pickle.load(f)

with open(path.join(results_dir, "erro_explicito"), 'rb') as f:
    erro_explicito = pickle.load(f)    
with open(path.join(results_dir, "erro_implicito"), 'rb') as f:
    erro_implicito = pickle.load(f)
with open(path.join(results_dir, "erro_RK4"), 'rb') as f:
    erro_RK4 = pickle.load(f)


fig = plt.figure()
ax = fig.add_subplot(121)

ax.plot(result_ref.t, result_ref.u[0,:], '.-', label="Solução de referência")
ax.plot(result_RK4.t, result_RK4.u[0,:], '.--', label="RK4")
ax.plot(result_explicito.t, result_explicito.u[0,:], '.--', label="Euler explícito")
ax.plot(result_implicito.t, result_implicito.u[0,:], '.--', label="Euler implícito")
ax.grid()
ax.legend()

detalhe = fig.add_subplot(122)
detalhe.plot(result_ref.t, result_ref.u[0,:], '.-', label="Solução de referência")
detalhe.plot(result_RK4.t, result_RK4.u[0,:], '.--', label="RK4")
detalhe.plot(result_explicito.t, result_explicito.u[0,:], '.--', label="Euler explícito")
detalhe.plot(result_implicito.t, result_implicito.u[0,:], '.--', label="Euler implícito")
detalhe.grid()
detalhe.set_xbound(result_ref.t[-1] - 2, result_ref.t[-1])
detalhe.set_ybound(0, 2)

fig.set_figheight(6)
fig.set_figwidth(11)
fig.savefig(path.join(images_dir, "Comparacao_p_simples.pdf"))
if show:
    plt.show()