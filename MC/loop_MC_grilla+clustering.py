#Este codigo corre toyMC para diferentes parametros de entrada a y b

import os,sys

command = ""

N1=200
N2=201

iDC1=0
iDC2=1

a= [1,10,100,1000]

b= [1,10,100,1000]

R1=1
R2=2

for N in range(N1,N2):
	for iDC in range(iDC1,iDC2):
		DC=iDC*10
		for A in a: 
			for B in b: 
				for R in range(R1,R2):  # Este loop determina cuantas veces se repite el MC con los mismos parametros
					command += "./toyMC.exe "+str(N)+" "+str(DC)+" "+str(A)+" "+str(B)+" "+str(R)+" && "
command = command[:-2]     
print command
os.system(command)

# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

command = ""
for N in range(N1,N2):
	for iDC in range(iDC1,iDC2):
		DC=iDC*10
		for A in a: 
			for B in b: 
				for R in range(R1,R2):
					command += "./skExtract.exe -c extractConfig.xml MC_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".fits -o output_N0="+str(N)+"_DC="+str(DC)+"_A="+str(A)+"_B="+str(B)+"_R="+str(R)+".root"+" && "
command = command[:-2]     
print command
os.system(command)
