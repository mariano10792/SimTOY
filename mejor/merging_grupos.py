from astropy.io import fits
import numpy as np
import os,sys


#this script merges fits file in groupos or 2,3,4 or whatever. The only condition is that they are named as they are numbered: "1.fits, 2.fits, etc"

for h in range(1,26):
	a=h #amount or files i wanna merge
	B=127/a # number of group of files merged
	#b=round(B)

	command = "mkdir grupos_de_"+str(a)+""
	os.system(command)
	for j in range(1,B+1):
		tt = np.zeros((50,500))	
		for i in range (1,a+1):
			k1=a*(j-1)+i
			hdu1 = fits.open('/home/mariano/Desktop/SimTOY/mejor/'+str(k1)+'.fits')
			t=hdu1[3].data
			tt=tt+t
		
		#T=np.sum(tt)
		#print (T)
		hdu1 = fits.PrimaryHDU(tt)
		hdu1l = fits.HDUList([hdu1])
		hdu1l.writeto('./grupos_de_'+str(a)+'/merge_'+str(k1-(a-1))+'_'+str(k1)+'.fits')
		command = "./skExtract.exe -c extractConfig.xml ./grupos_de_"+str(a)+"/merge_"+str(k1-(a-1))+"_"+str(k1)+".fits -o ./grupos_de_"+str(a)+"/merge_"+str(k1-(a-1))+"_"+str(k1)+".root -q"
		os.system(command)
	#command2= "hadd ./grupos_de_"+str(a)+"/merged_"+str(a)+".root ./grupos_de_"+str(a)+"/merge_"+str(k1-(a-1))+"_"+str(k1)+".root"
	command2= "hadd ./merged_"+str(a)+".root ./grupos_de_"+str(a)+"/merge_*.root"
	os.system(command2)
	
	
print("termine")
	




