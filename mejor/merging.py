from astropy.io import fits


arrayhdu = ()
count= 73
while (count < 267)

	
	array = (array, count)
	hdu = fits.open('/home/mariano/Desktop/caca/proc_fe552018-04-14_0142.fits')
	arrayhdu[count]=hdu
	
	
	
	
	
	
	
	


hdulist2 = fits.open('/home/mariano/Desktop/caca/proc_fe552018-04-14_0143.fits')


t1=hdulist[3].data
t2=hdulist2[3].data

T=t1+t2

hdu = fits.PrimaryHDU(T)
hdul = fits.HDUList([hdu])
hdul.writeto('merge142143.fits')



