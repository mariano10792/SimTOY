# Este codigo copia los archivos desde la carpeta de andre a una carpeta ficticia llamada caca en mi escritorio

import os,sys

command = ""
for a in range(73,86):

		command += "\n scp monsoon@estela:/home/monsoon/andre/macros/datos/tracer/newSetup/fe55/dataApril12/proc_fe552018-04-12_00"+str(a)+".fits /home/mariano/Desktop/caca/"
		#command += "./skExtract.exe -c extractConfig.xml interactions+dc_A="+str(a)+"_B="+str(b)+".fits -o output_"+str(a)+"_"+str(b)+".root"+" && "


command = command[:]

print command
os.system(command)
command = ""

for a in range(87,100):

		command += "\n scp monsoon@estela:/home/monsoon/andre/macros/datos/tracer/newSetup/fe55/dataApril12/proc_fe552018-04-13_00"+str(a)+".fits /home/mariano/Desktop/caca/"
		#command += "./skExtract.exe -c extractConfig.xml interactions+dc_A="+str(a)+"_B="+str(b)+".fits -o output_"+str(a)+"_"+str(b)+".root"+" && "


command = command[:]

print command
os.system(command)
command = ""

for a in range(100,117):

		command += "\n scp monsoon@estela:/home/monsoon/andre/macros/datos/tracer/newSetup/fe55/dataApril12/proc_fe552018-04-13_0"+str(a)+".fits /home/mariano/Desktop/caca/"
		#command += "./skExtract.exe -c extractConfig.xml interactions+dc_A="+str(a)+"_B="+str(b)+".fits -o output_"+str(a)+"_"+str(b)+".root"+" && "


command = command[:]

print command
os.system(command)
command = ""

for a in range(117,267):

		command += "\n scp monsoon@estela:/home/monsoon/andre/macros/datos/tracer/newSetup/fe55/dataApril12/proc_fe552018-04-14_0"+str(a)+".fits /home/mariano/Desktop/caca/"
		#command += "./skExtract.exe -c extractConfig.xml interactions+dc_A="+str(a)+"_B="+str(b)+".fits -o output_"+str(a)+"_"+str(b)+".root"+" && "


command = command[:]

print command
os.system(command)
