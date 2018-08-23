# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

import os,sys
#
command = ""
for N in range(100,101):
	for iDC in range(1,2):
		DC=iDC*543
		for A in [4000]: 
			for B in [4,7,10,13,16]: 
				for R in range(1,128):   # Este loop determina cuantas veces se repite el MC con los mismos parametros
					command += "./skExtract.exe -c extractConfig.xml MC_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".fits -o output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".root"+" && "
command = command[:-3]     
print command
os.system(command)
