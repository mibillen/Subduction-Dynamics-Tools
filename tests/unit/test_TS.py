#!/opt/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import python_scripts.TS as TS

def test_solve_angle_of_subduction_CD():
	'''
	implemented by Haoyuan
	test iplementation of solve_angle_of_subduction_CD
	'''
	## Arc-corner Values for a theta = pi/4 slab is derived in 6.115 and 6.116
	theta = np.pi/4.0
	tolerance = 1e-6
	Cstd = -np.pi * 2**0.5 / 2 / (2 - np.pi**2.0/4.0)
	Dstd = -2.0**0.5 * (2 - np.pi/2.0) / (2 - np.pi**2.0/4.0)
	C, D = TS.solve_angle_of_subduction_arc_CD(theta)
	assert(abs(C-Cstd)/Cstd < tolerance)
	assert(abs(D-Dstd)/Dstd < tolerance)
	print(C, D)
	## Arc-corner Values for a theta = pi/4 slab is derived in 6.122 and 6.123
	## Haoyuan: TODO

