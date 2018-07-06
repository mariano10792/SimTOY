# Este codigo suma usando hadd todas las salidas de skExtract.exe que correspondan a iguales parametros fisicos de entrada

import os,sys



command = ""
for N in range(200,201):
	for iDC in range(0,1):
		DC=iDC*10
		for A in [100,250,500,750,1000,1500,2000,2500,3000,5000]: 
			for B in [1,5,10,15,20,25,30]: 
				command += "hadd output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+".root output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R=*.root && "
command = command[:-2]     
print command
os.system(command)


