#Este codigo corre toyMC para diferentes parametros de entrada a y b

import os,sys

command = ""
for N in range(200,201):
	for iDC in range(0,1):
		DC=iDC*10
		for A in [10000]: 
			for B in [100]: 
				for R in range(1,10):  # Este loop determina cuantas veces se repite el MC con los mismos parametros
					command += "./toyMC.exe "+str(N)+" "+str(DC)+" "+str(A)+" "+str(B)+" "+str(R)+" && "
command = command[:-2]     
print command
os.system(command)
