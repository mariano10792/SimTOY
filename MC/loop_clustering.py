# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

import os,sys
#,,,,,,,,,5000]
command = ""
for N in range(200,201):
	for iDC in range(0,1):
		DC=iDC*10
		for A in [5000]: 
			for B in [1,5,10,15,20,25,30]: 
				for R in range(1,50): 
					command += "./skExtract.exe -c extractConfig.xml MC_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".fits -o output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".root"+" && "
command = command[:-2]     
print command
os.system(command)
