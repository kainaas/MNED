import numpy as np
from tabulate import tabulate

regressiva = lambda u,x,h: (u(x) - u(x-h))/h
progressiva = lambda u,x,h: (u(x+h) - u(x))/h
central = lambda u,x,h: (u(x+h)-u(x-h))/(2*h)

u1 = lambda x: np.sin(x)
u2 = lambda x: np.exp(-x)

u1_linha = lambda x: np.cos(x)
u2_linha = lambda x: -np.exp(-x)

x1 = 2
x2 = 1

cabecalho = ["h", "regressiva", "erro regressiva", "progressiva", "erro progressiva", "central", "erro central"]

tabela1 = []
for i in range(3):
    h = 10**(-(2**i))
    reg = regressiva(u1, x1, h)
    prog = progressiva(u1, x1, h)
    cent = central(u1,x1,h)
    erro_reg = reg - u1_linha(x1)
    erro_prog = prog - u1_linha(x1)
    erro_cent = cent - u1_linha(x1)
    tabela1.append([h, reg, erro_reg, prog, erro_prog, cent, erro_cent])
print("Tabela u(x) = sen(x)")
print(tabulate(tabela1, headers = cabecalho))

tabela2 = []
for i in range(3):
    h = 10**(-(2**i))
    reg = regressiva(u2, x2, h)
    prog = progressiva(u2, x2, h)
    cent = central(u2,x2,h)
    erro_reg = reg - u1_linha(x2)
    erro_prog = prog - u2_linha(x2)
    erro_cent = cent - u2_linha(x2)
    tabela2.append([h, reg, erro_reg, prog, erro_prog, cent, erro_cent])
print("\n\nTabela u(x) = exp(-x)")
print(tabulate(tabela2, headers = cabecalho))