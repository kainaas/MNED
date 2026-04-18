import numpy as np

h = 0.1
n = 3 -1
m = 2*(n+1) -1

N = n*(m+2)
A = np.zeros((N+2,N+2))
F = np.zeros((1,N+2))


# Bottom
k = 1
for i in range(1,n+1): # n iterations
    A[k, i-1] = 0 if (i == 1) else 1/h**2
    A[k, i] = -8/h**2
    A[k, i+1] = 0 if (i == n) else 1/h**2
    A[k, i+n] = 6/h**2

    F[0, k] = -(6/h + 4*np.pi**2 - 3)*np.e*np.sin(2*np.pi*i*h)

    k += 1

# Inner
for j in range(1,m+1):
    for i in range(1,n+1): # n*m iterations
        A[k, i-1+n*j] = 0 if (i == 1) else 1/h**2
        A[k, i+n*j] = -8/h**2
        A[k, i+1+n*j] = 0 if (i == n) else 1/h**2

        A[k, i+n*(j-1)] = 3/h**2 + (4*np.pi**2 - 3)/(2*h)
        A[k, i+n*(j+1)] = 3/h**2 - (4*np.pi**2 - 3)/(2*h)

        k +=1

# Top
for i in range(1,n+1): # n iterations
    A[k, i-1+n*(m+1)] = 0 if (i == 1) else 1/h**2
    A[k, i+n*(m+1)] = 4*np.pi**2-3-8/h**2-6/h
    A[k, i+1+n*(m+1)] = 0 if (i == n) else 1/h**2

    A[k, i+n*m] = 6/h**2

    k += 1

### My method

M = np.zeros((N + 2, N + 2))

#D0
i = np.arange(1, N + 1)
M[i,i] = -8/h**2

#D-
i = np.arange(2, N + 1)
j = i[i % n != 1]
M[j, j-1] = 1/h**2

#D+
i = np.arange(1, N)
j = i[i % n != 0]
M[j, j+1] = 1/h**2

#A-n
i = np.arange(n+1, N + 1)
alpha = 3/h**2 + (4*np.pi**2 - 3)/(2*h)
M[i, i-n] = alpha

#B+n
i = np.arange(1, N + 1 - n)
beta  = 3/h**2 - (4*np.pi**2 - 3)/(2*h)
M[i, i+n] = beta

#P-n
i = np.arange(n*(m+1)+1, N + 1)
M[i, i-n] += beta

#Q+n
i = np.arange(1, n+1)
M[i, i+n] += alpha

#G0
i = np.arange(n*(m+1)+1, N + 1)
M[i, i] += -2*h*beta

###

np.savetxt(
    "out.csv",
    A[1:N+1, 1:N+1],
    delimiter=",",
    fmt="%.10e"
)

np.savetxt(
    "out2.csv",
    M[1:N+1, 1:N+1],
    delimiter=",",
    fmt="%.10e"
)