import numpy as np
from math import floor

EPS = 10**-7 #Epsilon pro método de Newton
MAX_IT = 20 #Número máximo de iterações para Newton

#Classe para retornar o resultado dos métodos
class Resultado:
    def __init__(self, qnt_passos, h, u, t):
        self.qnt_passos = qnt_passos #Quantidade de passos executada
        self.h = h #Passo de t
        self.u = u #Vetor com a solução
        self.t = t #Vetor com os valores do tempo utilizados. Note que  t_final <= t[-1] < t_final + h


        

def calcula_qnt_passos(t0, t_final, h):
    tam_intervalo = t_final - t0

    qnt_passos = tam_intervalo/h
    if qnt_passos != floor(qnt_passos): qnt_passos = floor(qnt_passos) + 1
    else: qnt_passos = int(qnt_passos)
    qnt_passos += 1 #Acrescentar u0 no vetor solução
    return qnt_passos

def inicializa_metodo(qnt_passos, t0, u0):
    u = np.zeros((len(u0), qnt_passos)) #Cria uma matriz em que as colunas são a solução em um determinado tempo e as linhas são a solução de uma variável no tempo todo
    t = np.zeros(qnt_passos)

    u[:, 0] = u0
    t[0] = t0
    return u, t

def euler_explicito(t0, t_final, h, u0: np.array, f):
    qnt_passos = calcula_qnt_passos(t0, t_final, h)
    u, t = inicializa_metodo(qnt_passos, t0, u0)

    for i in range(1,qnt_passos):
        t[i] = t[i-1] + h
        u[:, i] = u[:, i-1] + h*f(u[:, i-1], t[i-1])

    return Resultado(qnt_passos, h, u, t)
     

def RK4(t0, t_final, h, u0: np.array, f):
    qnt_passos = calcula_qnt_passos(t0, t_final, h)
    u, t = inicializa_metodo(qnt_passos, t0, u0)

    for i in range(1,qnt_passos):
        t[i] = t[i-1] + h
        K1 = f(u[:, i-1], t[i-1])
        K2 = f(u[:, i-1] + (h/2)*K1, t[i-1] + h/2)
        K3 = f(u[:, i-1] + (h/2)*K2, t[i-1] + h/2)
        K4 = f(u[:, i-1] + h*K3, t[i-1] + h)
        u[:,i] = u[:,i-1] + (h/6)*(K1 + 2*K2 + 2*K3 + K4)

    return Resultado(qnt_passos, h, u, t)


def euler_implicito(t0, t_final, h, u0: np.array, f, Jf): #Jf é a matriz da derivada de f em relação ao vetor U^{n+1}
    qnt_passos = calcula_qnt_passos(t0, t_final, h)
    u, t = inicializa_metodo(qnt_passos, t0, u0)

    dim = len(u0)

    for i in range(1,qnt_passos):
        t[i] = t[i-1] + h
        u_it = u[:,i-1].copy()
        for j in range(MAX_IT): #Loop do método de Newton
             # G(x)
            G = u_it - u[:, i-1] - h*f(u_it, t[i])

            # Jacobiana de G
            JG = np.eye(dim) - h*Jf(u_it, t[i])

            delta = np.linalg.solve(JG, -G)

            u_it += delta
            if np.linalg.norm(delta) < EPS:
                u[:,i] = u_it   
                break
            if j == MAX_IT - 1:
                print(delta)

    return Resultado(qnt_passos, h, u, t)
    
        
def erro_infinito(analitica: Resultado, calculada: Resultado):
    return np.max(np.abs(analitica.u-calculada.u),axis=1)

def erro_norma2(analitica: Resultado, calculada: Resultado):
    #(h * soma(erro ao quadrado))^(1/2) - Definição do LeVeque de erro em funções de malha
    return np.sqrt(analitica.h * np.sum(np.square(analitica.u - calculada.u),axis=1))


if __name__ == "__main__":

    #Testando num problema simples todos os métodos
    h = 0.1
    t0 = 0
    t_final = 1
    u0 = np.array([1])
    f = lambda u,t: u
    Jf = lambda u,t: np.array([[1]])

    resultado_euler = euler_explicito(t0, t_final, h, u0, f)
    print("Euler:")
    print(resultado_euler.u)
    print(resultado_euler.t)

    resultado_RK4 = RK4(t0, t_final, h, u0, f)
    print("RK4")
    print(resultado_RK4.u)
    print(resultado_RK4.t)

    resultado_implicito = euler_implicito(t0, t_final, h, u0, f, Jf)
    print("Euler Implícito")
    print(resultado_implicito.u)
    print(resultado_implicito.t)