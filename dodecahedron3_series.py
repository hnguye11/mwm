# S-T UNRELIABILITY (UNDIRECTED GRAPH)

V = ['s', '1_a', '2_a', '3_a', '4_a', '5_a', '6_a', '7_a', '8_a', '9_a', '10_a', '11_a', '12_a', '13_a', '14_a', '15_a', '16_a', '17_a', '18_a', 't_a', '1_b', '2_b', '3_b', '4_b', '5_b', '6_b', '7_b', '8_b', '9_b', '10_b', '11_b', '12_b', '13_b', '14_b', '15_b', '16_b', '17_b', '18_b', 't_b', '1_c', '2_c', '3_c', '4_c', '5_c', '6_c', '7_c', '8_c', '9_c', '10_c', '11_c', '12_c', '13_c', '14_c', '15_c', '16_c', '17_c', '18_c', 't']

E = [('s', '1_a'), ('s', '2_a'), ('s', '3_a'), ('1_a', '5_a'), ('1_a', '6_a'), ('2_a', '7_a'), ('2_a', '8_a'), ('3_a', '4_a'), ('3_a', '9_a'), ('4_a', '5_a'), ('6_a', '7_a'), ('8_a', '9_a'), ('4_a', '10_a'), ('5_a', '11_a'), ('6_a', '12_a'), ('7_a', '13_a'), ('8_a', '14_a'), ('9_a', '15_a'), ('10_a', '15_a'), ('11_a', '12_a'), ('13_a', '14_a'), ('10_a', '16_a'), ('11_a', '16_a'), ('12_a', '17_a'), ('13_a', '17_a'), ('14_a', '18_a'), ('15_a', '18_a'), ('16_a', 't_a'), ('17_a', 't_a'), ('18_a', 't_a')] + [('t_a', '1_b'), ('t_a', '2_b'), ('t_a', '3_b'), ('1_b', '5_b'), ('1_b', '6_b'), ('2_b', '7_b'), ('2_b', '8_b'), ('3_b', '4_b'), ('3_b', '9_b'), ('4_b', '5_b'), ('6_b', '7_b'), ('8_b', '9_b'), ('4_b', '10_b'), ('5_b', '11_b'), ('6_b', '12_b'), ('7_b', '13_b'), ('8_b', '14_b'), ('9_b', '15_b'), ('10_b', '15_b'), ('11_b', '12_b'), ('13_b', '14_b'), ('10_b', '16_b'), ('11_b', '16_b'), ('12_b', '17_b'), ('13_b', '17_b'), ('14_b', '18_b'), ('15_b', '18_b'), ('16_b', 't_b'), ('17_b', 't_b'), ('18_b', 't_b')] + [('t_b', '1_c'), ('t_b', '2_c'), ('t_b', '3_c'), ('1_c', '5_c'), ('1_c', '6_c'), ('2_c', '7_c'), ('2_c', '8_c'), ('3_c', '4_c'), ('3_c', '9_c'), ('4_c', '5_c'), ('6_c', '7_c'), ('8_c', '9_c'), ('4_c', '10_c'), ('5_c', '11_c'), ('6_c', '12_c'), ('7_c', '13_c'), ('8_c', '14_c'), ('9_c', '15_c'), ('10_c', '15_c'), ('11_c', '12_c'), ('13_c', '14_c'), ('10_c', '16_c'), ('11_c', '16_c'), ('12_c', '17_c'), ('13_c', '17_c'), ('14_c', '18_c'), ('15_c', '18_c'), ('16_c', 't'), ('17_c', 't'), ('18_c', 't')]

VG = ['VG0', 'VG1', 'VG1', 'VG1', 'VG4', 'VG4', 'VG4', 'VG4', 'VG4', 'VG4', 'VG10', 'VG10', 'VG10', 'VG10', 'VG10', 'VG10', 'VG16', 'VG16', 'VG16', 'VG19', 'VG20', 'VG20', 'VG20', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG23', 'VG20', 'VG20', 'VG20', 'VG19', 'VG16', 'VG16', 'VG16', 'VG10', 'VG10', 'VG10', 'VG10', 'VG10', 'VG10', 'VG4', 'VG4', 'VG4', 'VG4', 'VG4', 'VG4', 'VG1', 'VG1', 'VG1', 'VG0']

EG_UNIQUE = ['VG16-VG19', 'VG23-VG23', 'VG0-VG1', 'VG19-VG20', 'VG4-VG4', 'VG1-VG4', 'VG10-VG16', 'VG20-VG23', 'VG10-VG10', 'VG10-VG4']

EG_IDX = [2, 2, 2, 5, 5, 5, 5, 5, 5, 4, 4, 4, 9, 9, 9, 9, 9, 9, 8, 8, 8, 6, 6, 6, 6, 6, 6, 0, 0, 0, 3, 3, 3, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7, 7, 7, 3, 3, 3, 0, 0, 0, 6, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9, 9, 9, 9, 9, 4, 4, 4, 5, 5, 5, 5, 5, 5, 2, 2, 2]

p_1 = [1 - 1e-2 for Ei in E]
p_2 = [1 - 1e-3 for Ei in E]
p_3 = [1 - 1e-4 for Ei in E]

''' MWM solutions. '''
gq_1 = [8.901659e-01, 9.659005e-01, 8.901659e-01, 8.901659e-01, 9.900000e-01, 9.659005e-01, 9.659005e-01, 9.659005e-01, 9.900000e-01, 9.659005e-01]
gq_2 = [8.547872e-01, 9.871275e-01, 8.547872e-01, 8.547872e-01, 9.990000e-01, 9.871275e-01, 9.871275e-01, 9.871275e-01, 9.990000e-01, 9.871275e-01]
gq_3 = [8.405608e-01, 9.956634e-01, 8.405608e-01, 8.405608e-01, 9.999000e-01, 9.956634e-01, 9.956634e-01, 9.956634e-01, 9.999000e-01, 9.956634e-01]
