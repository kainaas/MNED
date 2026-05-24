import numpy as np
import metodos
import matplotlib.pyplot as plt
import pickle
from os import path

#u = (q1, q2, p1, p2)
#A1 = (p1*p2 * sin(q1 - q2)) / (1+ sin^2(q1 - q2))
#A2 = {[p1^2 + 2p2^2 - 2p1*p2 * cos(q1-q2)] * sin(q1 - q2) * cos(q1 - q2)} / (1 + sin^2(q1 - q2))^2
'''
f(u,t) = [(p1 - p2 * cos(q1 - q2) / (2 - cos^2(q1 - q2)),
          (2p2 - p1 * cos(q1 -q2)) / (2 - cos^2(q1 - q2)),
          -A1 + A2 - 2sin(q1),
          A1 - A2 - sin(q2)]
'''
#u0 = (ang_ini1, ang_ini2, momento_ini1, momento_ini2)

'''
q1 = u[0]
q2 = u[1]
p1 = u[2]
p2 = u[3]
'''

results_dir = "resultados_numericos"

def converter_cartesiano(result: metodos.Resultado):
    x1 = np.sin(result.u[0, :])
    x2 = x1 + np.sin(result.u[1, :])
    y1 = -np.cos(result.u[0, :])
    y2 = y1 - np.cos(result.u[1, :])
    return x1, x2, y1, y2


#Experimento da cond inicial próxima
ang_ini1_1 = np.pi / 2
ang_ini2_1 = np.pi / 2

ang_ini1_2 = np.pi / 2
ang_ini2_2 = np.pi / 2 + 1e-2

#Experimento das cond iniciais perto do estado de equilíbrio
ang_ini1_eq = 0.2
ang_ini2_eq = 0.2

#Todos começam com velocidade nula
momento_ini1 = 0
momento_ini2 = 0


u0_1 = np.array([ang_ini1_1, ang_ini2_1, momento_ini1, momento_ini2])
u0_2 = np.array([ang_ini1_2, ang_ini2_2, momento_ini1, momento_ini2])
u0_eq = np.array([ang_ini1_eq, ang_ini2_eq, momento_ini1, momento_ini2])

h = 0.01
t0 = 0
t_final = 50

A1 = lambda u: (u[2]*u[3]*np.sin(u[0] - u[1])) / (1 + (np.sin(u[0] - u[1]))**2 )
A2 = lambda u: ((u[2]**2 + 2*u[3]**2 - 2*u[2]*u[3] * np.cos(u[0] - u[1]))*np.sin(u[0] - u[1]) * np.cos(u[0]- u[1])) / (1 + (np.sin(u[0] - u[1]))**2)**2

f = lambda u, t: np.array([
    (u[2] - u[3] * np.cos(u[0] - u[1])) / (2 - (np.cos(u[0] - u[1]))**2),
    (2*u[3] - u[2] * np.cos(u[0] - u[1])) / (2 - (np.cos(u[0] - u[1]))**2),
    -A1(u) + A2(u) - 2 * np.sin(u[0]),
    A1(u) - A2(u) - np.sin(u[1])
])


result_RK4_1 = metodos.RK4(t0, t_final, h, u0_1, f)
result_RK4_2 = metodos.RK4(t0, t_final, h, u0_2, f)
result_RK4_eq = metodos.RK4(t0, t_final, h, u0_eq, f)

with open(path.join(results_dir, "duplo_1"), "wb") as file:
    pickle.dump(result_RK4_1, file)
with open(path.join(results_dir, "duplo_2"), "wb") as file:
    pickle.dump(result_RK4_2, file)
with open(path.join(results_dir, "duplo_eq"), "wb") as file:
    pickle.dump(result_RK4_eq, file)


