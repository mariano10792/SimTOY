# Este codigo suma usando hadd todas las salidas de skExtract.exe que correspondan a iguales parametros fisicos de entrada

import os,sys



command = ""
for N in range(100,101):
	for iDC in range(1,2):
		DC=iDC*543
		for A in [5000]: 
			for B in [7]: 
				command += "hadd output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+".root output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R=*.root && "
command = command[:-3]     
print command
os.system(command)


