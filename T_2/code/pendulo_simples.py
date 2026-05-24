import numpy as np
import metodos
from RefSolution import RefSolution
import matplotlib.pyplot as plt
import pickle
from os import path

results_dir = "resultados_numericos"
refresh_pickles = True
indice_salvar_pickle = 1


#q'' = -sin(q)
#u' = f(u,t)
#Logo, fazemos u = (q, q')
#u0 = (posição inicial, velocidade inicial)
#f(u,t)=(q', -sin(q))
#      =(u2, -sin(u1))


h_values = [0.1, 0.05, 0.01, 0.005, 0.001]

#linhas 0 e 1 de "erros" são os erros de q e p na norma infinito
#linhas 2 e 3 de "erros" são os erros de q e p nas norma quadratica
#linha 4 tem os valor de h correspondente
erro_explicito = np.zeros((5,len(h_values)))
erro_implicito = np.zeros((5,len(h_values)))
erro_RK4 = np.zeros((5,len(h_values)))

pos_ini = np.pi / 4
vel_ini = 0

t0 = 0
t_final = 15

u0 = np.array([pos_ini, vel_ini])
f = lambda u,t: np.array([u[1], -np.sin(u[0])])

Jf = lambda u,t: np.array([[0, 1], [-np.cos(u[0]), 0]])

for i, h in enumerate(h_values):
    result_explicito = metodos.euler_explicito(t0, t_final, h, u0, f)
    result_implicito = metodos.euler_implicito(t0, t_final, h, u0, f, Jf)
    result_RK4 = metodos.RK4(t0, t_final, h, u0, f)
    result_ref = RefSolution(t0, t_final, h, pos_ini)
    
    
    erro_explicito[:2,i] = metodos.erro_infinito(result_ref, result_explicito)
    erro_explicito[2:4, i] = metodos.erro_norma2(result_ref, result_explicito)
    erro_explicito[4, i] = h

    erro_implicito[:2,i] = metodos.erro_infinito(result_ref, result_implicito)
    erro_implicito[2:4, i] = metodos.erro_norma2(result_ref, result_implicito)
    erro_implicito[4, i] = h

    erro_RK4[:2,i] = metodos.erro_infinito(result_ref, result_RK4)
    erro_RK4[2:4, i] = metodos.erro_norma2(result_ref, result_RK4)
    erro_RK4[4, i] = h

    if i == indice_salvar_pickle and refresh_pickles:
        with open(path.join(results_dir, "e_explicito"), 'wb') as file:
            pickle.dump(result_explicito, file)
        
        with open(path.join(results_dir, "e_implicito"), 'wb') as file:
            pickle.dump(result_implicito, file)

        with open(path.join(results_dir, "RK4"), 'wb') as file:
            pickle.dump(result_RK4, file)

        with open(path.join(results_dir, "ref"), 'wb') as file:
            pickle.dump(result_ref, file)


if refresh_pickles:
    with open(path.join(results_dir, "erro_explicito"), 'wb') as f:
        pickle.dump(erro_explicito, f)
    
    with open(path.join(results_dir, "erro_implicito"), 'wb') as f:
        pickle.dump(erro_implicito, f)

    with open(path.join(results_dir, "erro_RK4"), 'wb') as f:
        pickle.dump(erro_RK4, f)
    