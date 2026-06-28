import numpy as np
import matplotlib.pyplot as plt

Dx = 1

def f(x):
    return Dx**2/((2/x) + Dx)

X = np.linspace(-5,5, 50)
X_exp = 10**X
F = f(X_exp)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$\mathbb{P}e$')
ax.loglog(X_exp, F, label=r'$\frac{\Delta x^2}{2 \mathbb{P}e^{-1} + \Delta x}$')
ax.plot([10**(-5), 10**(5)], [Dx, Dx], label=r'$\Delta x$')
ax.legend()

fig.savefig(f'graphs/retricao_Pe.pdf', bbox_inches='tight', dpi=300)
plt.show()
