# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

import os,sys
#,
command = ""
for N in range(200,201):
	for iDC in range(1,10):
		DC=iDC*10
		for A in [3500]: 
			for B in [10]: 
				for R in range(2,50):   # Este loop determina cuantas veces se repite el MC con los mismos parametros
					#command += "rm MC_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".fits"+" && "
					#	para los pdf
					command += "rm CCD_N0"+str(N)+"DC"+str(DC)+"A"+str(A)+".000000B"+str(B)+".000000R"+str(R)+".pdf"+" && "
					# para los root
					#command += "rm output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".root"+" && "
					
					
command = command[:-2]     
print command
os.system(command)




