# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

import os,sys
#
command = ""
for N in range(2,11):
#	for iDC in range(1,2):
#		DC=iDC*543
#		for A in [5000]: 
#			for B in [7]: 
	for R in range(1,128):   # Este loop determina cuantas veces se repite el MC con los mismos parametros
					command += "./skExtract.exe -c extractConfig.xml ./grupos_de_"+str(N)+"/*.fits -o ./grupos_de_"+str(N)+"/*.root"+" && "
command = command[:-3]     
print command
os.system(command)
