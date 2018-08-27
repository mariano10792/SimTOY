from astropy.io import fits


#arrayhdu = ()

a=5 #amount or files i wanna merge
B=120/a # number of group of files merged
b=round(B)

for j in range(1,b) # 1 a 24
	
	hduinit = fits.open('/home/mariano/Desktop/SimTOY/mejor/'+j+'.fits')
	tinit=hduinit[3].data
	hduinit = fits.PrimaryHDU(tinit)
	hduinitl = fits.HDUList([hduinit])
	hduinitl.writeto('merge1.fits')
	for i in range(a-1,a+1): # 
	
	


#archivo n
hduinit = fits.open('/home/mariano/Desktop/SimTOY/mejor/'+n+'.fits')
tinit=hduinit[3].data
hduinit = fits.PrimaryHDU(tinit)
hduinitl = fits.HDUList([hduinit])
hduinitl.writeto('merge1.fits')
count=2
#while (count < 127):

#el inicio debe ser n+1 y el final n+a dodne a es n cantidad de archivos que quiero mergear
#y n el primer archivo a mergear
#for count range(2,127)
for count range(n+1,n+a)
        coun=count-1
        Coun=str(coun)
        Count=str(count)
        hdu0=fits.open('/home/mariano/Desktop/SimTOY/mejor/merge'+Coun+'.fits')
        t0=hdu0[0].data

	hdu = fits.open('/home/mariano/Desktop/SimTOY/mejor/'+Count+'.fits')
        t=hdu[3].data

        T=t0+t

#        arrayhdu[count]=hdu
	

        hdusave = fits.PrimaryHDU(T)
        hdusavel = fits.HDUList([hdusave])
        hdusavel.writeto('merge'+Count+'.fits')
	
	
	count+=1
	
#antes de terminar el for hay que borrar el primer archivo de merge!! (el "n0")
	
print("termine")
	


#hdulist2 = fits.open('/home/mariano/Desktop/caca/proc_fe552018-04-14_0143.fits')


#t1=hdulist[3].data
#t2=hdulist2[3].data

#T=t1+t2

#hdu = fits.PrimaryHDU(T)
#hdul = fits.HDUList([hdu])
#hdul.writeto('merge142143.fits')




