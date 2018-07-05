# Este codigo suma usando hadd todas las salidas de skExtract.exe que correspondan a iguales parametros fisicos de entrada

import os,sys



command = ""
for N in range(200,201):
	for iDC in range(0,1):
		DC=iDC*10
		for A in [10000]: 
			for B in [100]: 
				command += "hadd output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+".root output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R=*.root && "
command = command[:-2]     
print command
os.system(command)


