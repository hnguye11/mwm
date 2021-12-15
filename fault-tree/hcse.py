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

''' MWM solutions. '''
gq_4 = [0.001266, 0.000119, 0.000119, 0.000393, 0.000196, 0.000196, 0.000004]
gq_3 = [0.003170, 0.000218, 0.000218, 0.000646, 0.000297, 0.000297, 0.000007]
gq_2 = [0.005225, 0.000297, 0.000296, 0.000841, 0.000367, 0.000367, 0.000009]
gq_1 = [0.007373, 0.000363, 0.000363, 0.001009, 0.000423, 0.000423, 0.000011]
