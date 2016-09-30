"""
Description: This script is used to comput SOS
for multiple forward anb backward run FEP simulations
"""

import numpy as np
import math

windows = 100			# number of the windows of "alchemical" transformation 
states =   65			# each of intermediate state for each window  
replicas =  1			# number of run replicas simulations 
K = 0.0019872041		# Boltzman constant (kcal/mol/K)
T = 298.15
beta = 1/(K*T)

#forward transformation: xijf

# energy.frwd.dat: data file with energies (windows, states)
# temperature.frwd.dat: data file with respecitve temperatures (windows, states)


energy_frwd = np.loadtxt("energy.frwd.dat")
temp_frwd = np.loadtxt("temperature.frwd.dat")

if len(energy_frwd) == len(temp_frwd):
	pass
else:
	print("Error message: input files must have the same lenght.")

if replicas == 1: 

	energy = np.array(energy_frwd).reshape(windows, states)
	temp = np.array(temp_frwd).reshape(windows, states)

	xf = (-1/(2*K*temp))*(energy)	 

	xijf = np.mean(xf, axis=1) 		# average over all intermediate states for each window
	print(xijf)
	print(xijf.shape)
	
else:

	print("Number of run forward FEP simulations: ", replicas)
	energy = np.array(energy_frwd).reshape(replicas, windows, states)
	temp = np.array(temp_frwd).reshape(replicas, windows, states)
    
	xf = (-1/(2*K*temp))*(energy)  

	xi = np.mean(xf, axis=0)		# average over all intermediate states for each window
	xijf = np.mean(xi, axis=1)		# average over all replicas for each window 
	print(xij_f)
	print(xij_f.shape)


#backward transformationi: xijb

# energy.back.dat: data file with energies (windows, states)
# temperature.back.dat: data file with respecitve temperatures (windows, states)


energy_back = np.loadtxt("energy.back.dat")
temp_back = np.loadtxt("temperature.back.dat")

if len(energy_back) == len(temp_back):
    pass
else:
	print("Error message: input files must have the same lenght.")

if replicas == 1:

    energy = np.array(energy_back).reshape(windows, states)
    temp = np.array(temp_back).reshape(windows, states)

    xb = (-1/(2*K*temp))*(energy)

    xijb = np.mean(xf, axis=1)		# average over all intermediate states for each window
    print(xijb)
    print(xijb.shape)

else:

    print("Number of run backward FEP simulations: ", replicas)
    energy = np.array(energy_back).reshape(replicas, windows, states)
    temp = np.array(temp_back).reshape(replicas, windows, states)

    xf = (-1/(2*K*temp))*(energy)

    xi = np.mean(xf, axis=0)		# average over all intermediate states for each window
    xijb = np.mean(xi, axis=1)		# average over all replicas for each window
    print(xij_b)
    print(xij_b.shape)

# Average

if len(energy_frwd) == len(temp_frwd):
	xij_av = (xijf)/(xijb)    
else:
    print("Error message: forward and backward transformations must have the same lenght.")


# Computing the SOS total free energy (DeltaG)

ln = np.log(xij_av)
deltaG = (-beta)*(np.sum(ln))
print(deltaG)

# Standard errors (eps)

# forward 

eps =  

