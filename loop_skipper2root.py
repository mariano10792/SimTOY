# Este codigo corre skExtract.exe para diferentes outputs con parametros de entrada N0, DC, A y B
# y para cada run R

import os,sys
#
command = "" 

for R in range (1,61)
	for S in range(1,5):   # Este loop determina cuantas veces se repite el MC con los mismos parametros
					command += "./skipper2root.exe /home/ccds/Soft/ltaDaemon-eth/images/BsAs/30Sep2018/skp_55Fe_Mylar_-2000Samp-4line_"+str(R)+"_"+str(S)+".fits -o skp_55fe -o /home/ccds/Soft/ltaDaemon-eth/images/BsAs/30Sep2018/rootfiles/skp_55Fe_Mylar_-2000Samp-4line_"+str(R)+"_"+str(S)+" && "
command = command[:-3]     
print command
os.system(command)
