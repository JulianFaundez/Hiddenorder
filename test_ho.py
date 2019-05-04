import numpy as np

# Running main function for hidden order

step = 0.1
TD = 0.3    # termo de hopping entre bandas
J = 0.128   # termino de exchange entre bandas
KB = 8.6e2  # constante de boltzman
e1, e2, e3 = 0, 0, 0
NX, Ny, Ny = 55, 55, 55

k_argument = np.arange(-np.pi, np.pi, step)
q0 = 0.1

kx, ky, kz = np.meshgrid(k_argument, k_argument, k_argument)

Eak = -2 * TD * (np.cos(kx) + np.cos(ky) + np.cos(kz))
EakQ = -2 * TD * (np.cos(kx + np.pi) + np.cos(ky + np.pi) + np.cos(kz + np.pi))
Ebk = -2 * TD * (np.cos(kx) + np.cos(ky) + np.cos(kz))

eaux = J ** 2 * q0 ** 2

A = np.ones(Eak.shape)
B = - EakQ - Ebk
C = EakQ * Ebk - eaux

# C.shape = (63, 63, 63)

input_roots = np.array([A, B, C])
# (3, 63, 63, 63)
# [:, 0, 0, 0]

# trying order two
store_data = np.zeros((2, 63, 63, 63))

# dumb iteration, just to try
for i in range(63):
    for j in range(63):
        for k in range(63):
            store_data[:, i, j, k] = np.roots(input_roots[:, i, j, k])
