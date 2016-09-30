"""
Description: This script is used to comput SOS
for multiple forward anb backward run FEP simulations
"""

import numpy as np
import math

windows = 100			# number of the windows of "alchemical" transformation (i) 
replicas = 65			# number of intermediate states (replicas) for each window (j)  
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

	xijf = np.exp((-1/(2*K*temp))*(energy))	    # xijf 	 
	xijf2 = (xijf)**(2)							# xijf squared
	xif_av = np.mean(xijf, axis=1) 				# average over all replicas 
	xif_av2 = (xif_av)**(2) 					# average over all replicas squared	
#	print(xijf)
#	print(xijf.shape)
#	print(xijf2)
#	print(xif_av)
#	print(xif_av.shape)
	

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
#	print(energy)
#	print(energy.shape)


	xijb = np.exp((-1/(2*K*temp))*(energy))      # xijb
	xijb2 = (xijb)**(2)                          # xijb squared
	xib_av = np.mean(xijb, axis=1)				 # average over all replicas
	xib_av2 = (xib_av)**(2)						 # average over all replicas squared	
#	print(xijb)
#	print(xijb.shape)
#	print(xib_av)
#	print(xib_av.shape)


# Computing the SOS total free energy (DeltaA)

xi_av = (xif_av)/(xib_av)
#print(xi_av)
ln = np.log(xi_av)
deltaA = (-1/beta)*(np.sum(ln))
print(ln)
print (deltaA)

# Standard errors (e2)

# forward transformation 

eijf2 = (np.sum(xijf2, axis=1))/(replicas**(2))
eif_av2 = (xif_av2)/(replicas)
e2f = eijf2 - eif_av2
#print(xijf2)
#print(xijf2.shape)
#print(xif_av2)
#print(xif_av2.shape)
#print(e2f)
#print(e2f.shape)

# backward transformation

eijb2 = (np.sum(xijb2, axis=1))/(replicas**(2))
eib_av2 = (xib_av2)/(replicas)
e2b = eijb2 - eib_av2
#print(e2b)

# Standard error for deltaA

ef = (e2f/xif_av)**(2) 
eb = (e2b/xib_av)**(2)
#print(ef)
#print(eb)
e = ef + eb
s = (np.sum(e, axis=0))
sqrt = np.sqrt(s)
#print(s)
#print(sqrt)
error = (1/beta)*(sqrt)
#print(deltaA)
print(error)

with open("fep.dat", "w") as fout:
    fout.write("{:.15f}".format(deltaA)) 
