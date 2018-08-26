from astropy.io import fits


#arrayhdu = ()


hduinit = fits.open('/home/mariano/Desktop/SimTOY/mejor/1.fits')
tinit=hduinit[3].data
hduinit = fits.PrimaryHDU(tinit)
hduinitl = fits.HDUList([hduinit])
hduinitl.writeto('merge1.fits')



count=2
while (count < 127):
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
	
	
print("termine")
	


#hdulist2 = fits.open('/home/mariano/Desktop/caca/proc_fe552018-04-14_0143.fits')


#t1=hdulist[3].data
#t2=hdulist2[3].data

#T=t1+t2

#hdu = fits.PrimaryHDU(T)
#hdul = fits.HDUList([hdu])
#hdul.writeto('merge142143.fits')




