# Este codigo suma usando hadd todas las salidas de skExtract.exe que correspondan a iguales parametros fisicos de entrada

import os,sys



command = ""
for N in range(1,2):
	for iDC in range(1,10):
		DC=iDC*100
		for A in [3500]: 
			for B in [10]: 
				command += "hadd output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+".root output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R=*.root && "
command = command[:-2]     
print command
os.system(command)


