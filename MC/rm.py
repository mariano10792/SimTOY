# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

import os,sys
#,
command = ""
for N in range(200,201):
	for iDC in range(0,1):
		DC=iDC*10
		for A in [4000,4250,4500,4750,5000]: 
			for B in [1,5,10,15,20,25,30]: 
				for R in range(1,50): 
					#command += "rm MC_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".fits"+" && "
					#	para los pdf
					#command += "rm CCD_N0"+str(N)+"DC"+str(DC)+"A"+str(A)+".000000B"+str(B)+".000000R"+str(R)+".pdf"+" && "
					# para los root
					command += "rm output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".root"+" && "
					
					
command = command[:-2]     
print command
os.system(command)




