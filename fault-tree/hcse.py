lbd_hw = 5e-5                   # Hardware
lbd_mi = 5e-5                   # Memory Interface
lbd_m = 6e-5                    # Memory
lbd_p = 1e-4                    # Processor
lbd_ps = 1e-4                   # Processor Swap
lbd_b = 1e-6                    # Bus
N1 = 4                          # Number of replicas
k = 2                           # Mininum number of working replicas
T = gmm = 100                   # Mission time
lbd = [lbd_hw] + [lbd_m] * 5 + [lbd_mi] * 2 + [lbd_p] * 2 + [lbd_ps] + [lbd_b] * 2
EG_IDX = [0, 1, 1, 2, 1, 1, 3, 3, 4, 4, 5, 6, 6]
lbd *= N1
EG_IDX *= N1
M = len(lbd)
M1 = len(set(EG_IDX))
