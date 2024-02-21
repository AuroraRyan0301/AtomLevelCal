import numpy as np
import matplotlib.pyplot as plt
# from utils import *
from constants import *


laser_power = 1000 # W
hbar = 1.0545718e-34 # J s
dipole_moment1 = 3.461079186*ea0 # C m
electic_field = np.sqrt(2 * laser_power / (c * epsion0)) # V/m
# import ipdb; ipdb.set_trace()
print(dipole_moment1)
rabi_freq1 = dipole_moment1 * electic_field / hbar

# import ipdb; ipdb.set_trace()
oma = 1
omb = 1
dels = 0
N = 3
R = 401

del_vals = np.zeros(R)
M = np.zeros((N**2, N**2),dtype = 'complex_')
rho = np.zeros((N, N),dtype = 'complex_')
Ham = np.zeros((N, N),dtype = 'complex_')
Q = np.zeros((N, N),dtype = 'complex_')
W = np.zeros((N**2 - 1, N**2 - 1),dtype = 'complex_')
S = np.zeros((N**2 - 1, 1),dtype = 'complex_')
B = np.zeros((N**2 - 1, 1),dtype = 'complex_')
A = np.zeros((N**2, R),dtype = 'complex_')

for m in range(1,R+1):
    del_vals[m-1] = (m - (R + 1) / 2) / 10
    Ham = np.array([[del_vals[m-1]/2, 0, oma/2],
                    [0, -del_vals[m-1]/2, omb/2],
                    [oma/2, omb/2, (dels + 0.5j) * (-1)]])
    for n in range(1,N**2+1):
        for p in range(1,N**2+1):
            print("n {}, p {} ".format(n, p))
            beta = n % N if n % N != 0 else N
            alpha = (1 + (n - beta) / N)
            sigma = p % N if p % N != 0 else N
            epsilon = (1 + (p - sigma) / N)
            rho = np.zeros((N, N))
            rho[int(epsilon)-1, int(sigma)-1] = 1
            #print epsilon, sigma
            # print("epsilon {}, sigma {} ".format(epsilon, sigma))
            Q = (Ham @ rho - rho @ np.conj(Ham)) * (0 - 1j)
            # print(Q)
            Q[0, 0] += rho[2, 2] / 2
            Q[1, 1] += rho[2, 2] / 2
            M[int(n-1), int(p-1)] = Q[int(alpha - 1), int(beta - 1)]
            # print(M)
            # if n==2 and p==2:
            #     import ipdb; ipdb.set_trace()
            # print("alpha {}, beta {} ".format(alpha, beta))
    S = M[:N**2 - 1, N**2 - 1:]
    W = M[:N**2 - 1, :N**2 - 1]
    import ipdb; ipdb.set_trace()
    for d in range(1,N):
        W[:, ((d-1) * N + d-1)] -= S.flatten()

    B = np.linalg.solve(W, -S)
    rhonn = 1
    for f in range(1, N):
        rhonn -= B[((f-1) * N + f-1),0]

    A[:N**2 - 1, m-1] = B.flatten()
    A[N**2 - 1, m-1] = rhonn

plt.plot(del_vals, np.real(A[N**2 - 1, :]))
plt.show()