# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

import os,sys
#,
command = ""
for N in range(100,101):
	for iDC in range(1,2):
		DC=iDC*543
		for A in [3800,4000]: 
			for B in [4,7,10,13,16]: 
				for R in range(1,128):   # Este loop determina cuantas veces se repite el MC con los mismos parametros
					#command += "rm MC_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".fits"+" && "
					#	para los pdf
					#command += "rm CCD_N0"+str(N)+"DC"+str(DC)+"A"+str(A)+".000000B"+str(B)+".000000R"+str(R)+".pdf"+" && "
					# para los root
					command += "rm output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".root"+" && "
					
					
command = command[:-3]     
print command
os.system(command)




