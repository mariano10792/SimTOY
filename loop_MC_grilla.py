#Este codigo corre toyMC para diferentes parametros de entrada a y b

import os,sys

command = ""
for N in range(100,101):
	for iDC in range(1,2):
		DC=iDC*543
		for A in [5000]: 
			for B in [7]: 
				for R in range(1,128):   # Este loop determina cuantas veces se repite el MC con los mismos parametros
					command += "./toyMC.exe "+str(N)+" "+str(DC)+" "+str(A)+" "+str(B)+" "+str(R)+" && "
command = command[:-3]     
print command
os.system(command)
