#Este codigo corre toyMC para diferentes parametros de entrada a y b

import os,sys

command = ""
for N in range(1,2):
	for iDC in range(1,10):
		DC=iDC*100
		for A in [3500]: 
			for B in [10]: 
				for R in range(1,100):   # Este loop determina cuantas veces se repite el MC con los mismos parametros
					command += "./toyMC.exe "+str(N)+" "+str(DC)+" "+str(A)+" "+str(B)+" "+str(R)+" && "
command = command[:-2]     
print command
os.system(command)
