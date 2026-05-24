import numpy as np
from V15491650_metodos import Resultado
import matplotlib.pyplot as plt
import pickle
from os import path, makedirs
from V15491650_pendulo_duplo import converter_cartesiano

results_dir = "resultados_numericos"
images_dir = "images"
makedirs(results_dir, exist_ok=True)
makedirs(images_dir, exist_ok=True)
show = False

with open(path.join(results_dir, "duplo"), 'rb') as f:
    result = pickle.load(f)

with open(path.join(results_dir, "duplo_1"), 'rb') as f:
    result_1 = pickle.load(f)

with open(path.join(results_dir, "duplo_2"), 'rb') as f:
    result_2 = pickle.load(f)

with open(path.join(results_dir, "duplo_eq"), 'rb') as f:
    result_eq = pickle.load(f)


def plot_trajetoria(result: Resultado, title, savefig: bool, save_name, color_in):
    x1, x2, y1, y2 = converter_cartesiano(result)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x2, y2, color=color_in)
    ax.grid()
    #ax.set_title(title)
    if savefig:
        fig.savefig(path.join(images_dir, save_name),bbox_inches='tight')

plot_trajetoria(result, "Trajetória do pêndulo", True, "pendulo_duplo.pdf", 'red')
plot_trajetoria(result_1, "Trajetória do pêndulo - Condição inicial 1", True, "pendulo_duplo1.pdf", 'red')
plot_trajetoria(result_2, "Trajetória do pêndulo - Condição inicial 2", True, "pendulo_duplo2.pdf", 'red')
plot_trajetoria(result_eq, "Trajetória do pêndulo - próximo ao equilíbrio", True, "pendulo_duploeq.pdf", 'red')

if show:
    plt.show()