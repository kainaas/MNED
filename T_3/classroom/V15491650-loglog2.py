import numpy as np
import matplotlib.pyplot as plt
from os import makedirs

makedirs("graphs", exist_ok=True)

def f(x, Pe):
    return x**2/((2/Pe) + x)

X = np.linspace(-5,1, 50)
X_exp = 10**X


P=[0.001, 1000000]

for i, Pe in enumerate(P):
    F = f(X_exp,Pe)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$\Delta x$')
    label_str1 = r'$\frac{\Delta x^2}{2 \mathbb{P}e^{-1} + \Delta x}; \mathbb{P}e = 0.001$'
    label_str2 = r'$\frac{\Delta x^2}{2 \mathbb{P}e^{-1} + \Delta x}; \mathbb{P}e = 1000000$'
    
    if i == 0:
        ax.loglog(X_exp, F, label=label_str1)
        ax.loglog([10**(-4), 10**(-2)], [10**(-8), 10**(-4)], '--r', label=r'$f(x) = x^2$')
    elif i == 1:
        ax.loglog(X_exp, F, label=label_str2)
        ax.loglog([10**(-4), 10**(-2)], [1.5*10**(-4), 1.5*10**(-2)], '--r', label=r'$f(x) = x$')
    ax.legend()

    fig.savefig(f'graphs/retricao_Pe={Pe}.pdf', bbox_inches='tight', dpi=300)
    plt.show()
