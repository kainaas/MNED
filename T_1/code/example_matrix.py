import numpy as np

h = 0.1
n = 4 -1
m = 2*(n+1) -1

N = n*(m+2)
A = np.zeros((N,N))
F = np.zeros((1,N))

np.set_printoptions(linewidth=np.inf)

# Top
k = 0
for i in range(0,n): # n iterations
    if (i != 0): A[k, i-1] = 1/h**2
    A[k, i] = -8/h**2
    if (i != n-1): A[k, i+1] = 1/h**2
    A[k, i+n] = 6/h**2

    F[0, k] = -(6/h + 4*np.pi**2 - 3)*np.e*np.sin(2*np.pi*i*h)

    k += 1

# Inner
for j in range(1,m+1):
    for i in range(0,n): # n*m iterations
        if (i != 0): A[k, i-1+n*j] = 1/h**2
        A[k, i+n*j] = -8/h**2
        if (i != n-1): A[k, i+1+n*j] = 1/h**2

        A[k, i+n*(j-1)] = 3/h**2 + (4*np.pi**2 - 3)/(2*h)
        A[k, i+n*(j+1)] = 3/h**2 - (4*np.pi**2 - 3)/(2*h)

        k +=1

# Bottom
for i in range(0,n): # n iterations
    if (i != 0): A[k, i-1+n*(m+1)] = 1/h**2
    A[k, i+n*(m+1)] = 4*np.pi**2-3-8/h**2-6/h
    if (i != n-1): A[k, i+1+n*(m+1)] = 1/h**2

    A[k, i+n*m] = 6/h**2

    k += 1

### My method

M = np.zeros((N, N))
In = np.eye(n)
Im2 = np.eye(m+2)

alpha = 3/h**2 + (4*np.pi**2 - 3)/(2*h)
beta  = 3/h**2 - (4*np.pi**2 - 3)/(2*h)

# D- + D0 + D+
Tn = np.zeros((n, n))
i = np.arange(0, n)
Tn[i[1:], i[1:]-1] = 1/h**2 #D-
Tn[i, i] = -4/h**2 #D0/2
Tn[i[:-1], i[:-1]+1] = 1/h**2 #D+

# A-n + B+n
Tm2 = np.zeros((m+2, m+2))
i = np.arange(0, m+2)
Tm2[i[1:], i[1:]-1] = alpha #A-n
Tm2[i, i] = -4/h**2 #D0/2
Tm2[i[:-1], i[:-1]+1] = beta #B+n

M = np.kron(Im2, Tn) + np.kron(Tm2, In)

#P-n
i = np.arange(0, N)
j = i[n*(m+1):]
M[j, j-n] += beta

#G0
M[j, j] += -2*h*beta

#Q+n
j = i[:n]
M[j, j+n] += alpha

###

np.savetxt(
    "out.csv",
    A,
    delimiter=",",
    fmt="%.10e"
)

np.savetxt(
    "out2.csv",
    M,
    delimiter=",",
    fmt="%.10e"
)