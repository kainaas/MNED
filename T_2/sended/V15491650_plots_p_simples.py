import pickle
import numpy as np
import matplotlib.pyplot as plt
from os import path, makedirs
from V15491650_metodos import Resultado

results_dir = "resultados_numericos"
images_dir = "images"
makedirs(results_dir, exist_ok=True)
makedirs(images_dir, exist_ok=True)
show = False

#Abre todos os pickles de resultados
with open(path.join(results_dir, "e_explicito"), 'rb') as f:
    result_explicito = pickle.load(f)
with open(path.join(results_dir, "e_implicito"), 'rb') as f:
    result_implicito = pickle.load(f)
with open(path.join(results_dir, "RK4"), 'rb') as f:
    result_RK4 = pickle.load(f)
with open(path.join(results_dir, "ref"), 'rb') as f:
    result_ref = pickle.load(f)

#Abre todos os pickles dos erros
with open(path.join(results_dir, "erro_explicito"), 'rb') as f:
    erro_explicito = pickle.load(f)    
with open(path.join(results_dir, "erro_implicito"), 'rb') as f:
    erro_implicito = pickle.load(f)
with open(path.join(results_dir, "erro_RK4"), 'rb') as f:
    erro_RK4 = pickle.load(f)


#Função que define o plot dos gráficos de posição x tempo. Usado para fazer um grid de plots 2x2
fig_posxtempo = plt.figure()
fig_posxtempo.set_size_inches((10,10))
#fig_posxtempo.suptitle("Posição angular em função do tempo")
def plot_pos_angular(result: Resultado, title, pos_subplot, savefig: bool, save_name, color_in):
    ax = fig_posxtempo.add_subplot(pos_subplot)
    ax.plot(result.t, result.u[0,:], '.--', color=color_in)
    ax.grid()
    ax.set_xlabel("Tempo (s)")
    ax.set_ylabel("Posição angular (rad)")
    ax.set_title(title)
    if savefig:
        fig_posxtempo.savefig(path.join(images_dir, save_name),bbox_inches='tight')
    if show:
        plt.show()


#Função que define o plot dos gráficos de fase. Usado para fazer um grid de plots 2x2
fig_fase = plt.figure()
fig_fase.set_size_inches((10,10))
#fig_fase.suptitle("Retrato de fase")
def plot_fase(result: Resultado, title, pos_subplot, savefig: bool, save_name, color_in):
    ax = fig_fase.add_subplot(pos_subplot)
    ax.plot(result.u[0,:], result.u[1,:], color=color_in)
    ax.grid()
    ax.set_xlabel("Posição angular (rad)")
    ax.set_ylabel("Velocidade angular (rad/s)")
    ax.set_title(title)
    if savefig:
        fig_fase.savefig(path.join(images_dir, save_name),bbox_inches='tight')
    if show:
        plt.show()


#=======================================================================
#Gráfico que plota todos os valores no mesmo lugar. Bom para comparações
#=======================================================================
fig_comparacao = plt.figure()
#fig_comparacao.suptitle("Comparação dos resultados - Tempo x posição angular")
ax = fig_comparacao.add_subplot(121)

ax.plot(result_ref.t, result_ref.u[0,:], '.-', label="Solução de referência")
ax.plot(result_RK4.t, result_RK4.u[0,:], '.', label="RK4")
ax.plot(result_explicito.t, result_explicito.u[0,:], '.--', label="Euler explícito")
ax.plot(result_implicito.t, result_implicito.u[0,:], '.--', label="Euler implícito")
ax.set_xlabel("Tempo (s)")
ax.set_ylabel("Posição angular (rad)")
ax.grid()
ax.legend()

detalhe = fig_comparacao.add_subplot(122)
detalhe.plot(result_ref.t, result_ref.u[0,:], '.-', label="Solução de referência")
detalhe.plot(result_RK4.t, result_RK4.u[0,:], '.', label="RK4")
detalhe.plot(result_explicito.t, result_explicito.u[0,:], '.--', label="Euler explícito")
detalhe.plot(result_implicito.t, result_implicito.u[0,:], '.--', label="Euler implícito")
detalhe.grid()
detalhe.set_xlabel("Tempo (s)")
detalhe.set_ylabel("Posição angular (rad)")
detalhe.set_xbound(result_ref.t[-1] - 2, result_ref.t[-1]-1)
detalhe.set_ybound(0, 1.5)

fig_comparacao.set_figheight(6)
fig_comparacao.set_figwidth(11)
fig_comparacao.savefig(path.join(images_dir, "Comparacao_p_simples.pdf"),bbox_inches='tight')
if show:
    plt.show()


#============================
#Executa as funções dos plots
#============================
plot_pos_angular(result_explicito, "Resultado numérico - Euler explícito", 221, False, "p_simples_explicito_pos_angxtempo.pdf", 'green')
plot_pos_angular(result_implicito, "Resultado numérico - Euler implícito", 222, False, "p_simples_implicito_pos_angxtempo.pdf", 'red')
plot_pos_angular(result_RK4, "Resultado numérico - Runge-Kutta clássico", 223, False, "p_simples_RK4_pos_angxtempo.pdf", 'orange')
plot_pos_angular(result_ref, "Solução de referência", 224, True,"p_posxtempo.pdf", 'blue')

plot_fase(result_explicito, "Euler explícito", 221, False, "", 'green')
plot_fase(result_implicito, "Euler implícito", 222, False, "", 'red')
plot_fase(result_RK4, "Runge-Kutta clássico", 223, False, "",'orange')
plot_fase(result_ref, "Solução de referência", 224, True, "retrato_fase.pdf", 'green')


#=====================
#plots log-log do erro
#(norma infinito)
#(Posição angular)
#=====================
def plot_loglog(erro, title, savefig: bool, save_name, color_in, ordem_esperada: int):
    h_interval = erro[4,-1] - erro[4,0]
    x = np.array([1e-2, 2e-2])
    y = 10*x**ordem_esperada
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(erro[4,:], erro[0,:], color=color_in, marker='v', label='erro calculado')
    label_line = 0
    if ordem_esperada == 1:
        label_line = r'$O(h)$'
    else:
        label_line = r'$O(h^{%d})$' % ordem_esperada
    ax.loglog(x,y, color='black', ls='--', label=label_line)
    ax.set_xlabel("h")
    ax.set_ylabel("Erro - norma infinito")
    #ax.set_title(title)
    ax.grid(visible=True,which='minor')
    ax.legend()
    if savefig:
        fig.savefig(path.join(images_dir, save_name),bbox_inches='tight')
    if show:
        plt.show()

plot_loglog(erro_explicito, "Erro - Euler explícito", True, "erro_explicito.pdf", 'red', 1)
plot_loglog(erro_implicito, "Erro - Euler implícito", True, "erro_implicito.pdf", 'red', 1)
plot_loglog(erro_RK4, "Erro - Runge-Kutta clássico", True, "erro_RK4.pdf", 'red', 4)