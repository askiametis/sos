"""
Description: This script is used to comput SOS
for multiple forward anb backward run FEP simulations
"""

import numpy as np
import math

windows = 100			# number of the windows of "alchemical" transformation 
replicas = 65			# number of intermediate states (replicas) for each window  
run =  1				# number of run simulations 
K = 0.0019872041		# Boltzman constant (kcal/mol/K)
T = 298.15
beta = 1/(K*T)

#forward transformation:

# energy.frwd.dat: data file with energies (windows, replicas)
# temperature.frwd.dat: data file with respecitve temperatures (windows, replicas)


energy_frwd = np.loadtxt("energy.frwd.dat")
temp_frwd = np.loadtxt("temperature.frwd.dat")

if len(energy_frwd) == len(temp_frwd):
	pass
else:
	print("Error message: input files must have the same lenght.")

if run == 1: 

	energy = np.array(energy_frwd).reshape(windows, replicas)
	temp = np.array(temp_frwd).reshape(windows, replicas)

	xijf = ((-1)/(2*K*temp))*(energy)			# xijf 	 
	xijf2 = (xijf)**(2)							# xijf squared
	xif_av = np.mean(xijf, axis=1) 				# average over all replicas 
	xif_av2 = (xif_av)**(2) 					# average over all replicas squared	
#	print(xijf)
#	print(xijf2)
	print(xif_av)
#	print(xif_av2)
	

#backward transformation:

# energy.back.dat: data file with energies (windows, replicas)
# temperature.back.dat: data file with respecitve temperatures (windows, replicas)


energy_back = np.loadtxt("energy.back.dat")
temp_back = np.loadtxt("temperature.back.dat")

if len(energy_back) == len(temp_back):
    pass
else:
	print("Error message: input files must have the same lenght.")

if run == 1:

	energy = np.array(energy_back).reshape(windows, replicas)
	temp = np.array(temp_back).reshape(windows, replicas)

	xijb = (-1/(2*K*temp))*(energy)              # xijb
	xijb2 = (xijb)**(2)                          # xijb squared
	xib_av = np.mean(xijb, axis=1)				 # average over all replicas
	xib_av2 = (xib_av)**(2)						 # average over all replicas squared
print(xib_av)


# Average

if len(energy_frwd) == len(temp_frwd):
	xi_av = (xif_av)/(xib_av)
else:
    print("Error message: forward and backward transformations must have the same lenght.")

print(xi_av)
# Computing the SOS total free energy (DeltaA)

ln = np.log(xi_av)
deltaA = (-1/beta)*(np.sum(ln))
print(deltaA)

# Standard errors (e2)

# forward transformation 

eijf2 = (np.sum(xijf2))/(replicas**(2))
eif_av2 = (xif_av2)/(replicas)
e2f = eijf2 - eif_av2
#print(e2f)

# backward transformation

eijb2 = (np.sum(xijb2))/(replicas**(2))
eib_av2 = (xib_av2)/(replicas)
e2b = eijb2 - eib_av2
#print(e2b)

# Standard error for deltaA

ef = (e2f/xif_av)**(2) 
eb = (e2b/xib_av)**(2)
e = ef + eb
s = (np.sum(e, axis=0))
sqrt = np.sqrt(s)
#print(e)
#print(s)
#print(sqrt)
error = (1/beta)*(sqrt)
#print(deltaA)
#print(error)
