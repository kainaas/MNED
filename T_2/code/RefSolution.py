''' Program designed to express q as a function of t:

 Input:
 - q_0: Amplitude of the oscillation.
 - t: discretized time.'''

from scipy.special import ellipj
from scipy.special import ellipk
import numpy as np

class Resultado:
    def __init__(self, qnt_passos, h, u, t):
        self.qnt_passos = qnt_passos #Quantidade de passos executada
        self.h = h #Passo de t
        self.u = u #Vetor com a solução
        self.t = t #Vetor com os valores do tempo utilizados. Note que  t_final <= t[-1] < t_final + h




def RefSolution(t0, t_final, h, q_0):
    qnt_passos = int((t_final-t0)/h)
    t = np.linspace(t0, t_final, qnt_passos+1)

    # Auxiliary variable  
    k_0 = np.sin(0.5*q_0)  
 
    # Incomplete elliptic integral of the first kind 
    K = ellipk(k_0**2)
 
    # Jacobi elliptic functions: sn and cn  
    sn, cn, dn, ph = ellipj(K - t, k_0**2) 
 
    # Angular displacement  
    q = 2*np.arcsin(k_0*sn)
 
    # Angular velocity  
    p = -2*k_0*cn

    u = np.vstack([q,p])
    
    return Resultado(qnt_passos, h, u, t)