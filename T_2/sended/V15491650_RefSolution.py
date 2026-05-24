''' Program designed to express q as a function of t:

 Input:
 - q_0: Amplitude of the oscillation.
 - t: discretized time.'''

from scipy.special import ellipj
from scipy.special import ellipk
import numpy as np
from V15491650_metodos import Resultado


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