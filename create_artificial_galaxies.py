"""
Program to create artificial galaxies in random fields and positions

By Sandra Nogueira dos Reis, June 2017 @ IA - Lisbon
"""

from astropy.io import fits
import numpy as np
import os

#--------------------------------paths-----------------------------------------
images_path    = "/home/sreis/candels_images/"

#---------------------------------constants------------------------------------
zp_imaging = 25.95
pix_scale  = 0.06 #[arcsec pix^-1] for H-band; for I-band is 0.03
size_of_fit = 200 #[pix]: for H-band is 200x200 and for I-band is 400x400
mag_min = 15 #[mag]
mag_max = 19 #[mag]
mag_step = 5
re_min = 0.3/pix_scale #[arcsec]/[arcsec pix^-1] = [pix] => must be in [pix] for GALFIT
re_max = 3.3/pix_scale #[pix]
re_step = 10
nn_min = 2.3
nn_max = 7.3
nn_step = 6
ar_min = 0.4
ar_max = 1.
ar_step = 7
fields_names = ['goodsn','cosmos','uds','egs']
candels_fields = {'goodsn':images_path+"goodsn_all_wfc3_ir_f160w_060mas_v1.0",'cosmos':images_path+"hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0",'uds':images_path+"hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0","egs":images_path+"egs_all_wfc3_ir_f160w_060mas_v1.0"}
#------------------------------------------------------------------------------

mag_range = np.linspace(mag_min,mag_max,mag_step) #linspace return evenly spaced numbers over a specified interval
re_range = np.linspace(re_min,re_max,re_step)
nn_range = np.linspace(nn_min,nn_max,nn_step)
ar_range = np.linspace(ar_min,ar_max,ar_step)

count = 0 #to count the iterations
for mag_ele in mag_range:
    for re_ele in re_range:
        for nn_ele in nn_range:
            for ar_ele in ar_range:
                count = count+1
                field = np.random.choice(fields_names) #choose a random field
                seg_img = fits.open(field+"_segmentation.fits")
                seg_data = seg_img[0].data #read the values from segmentation file
                seg_hdr = seg_img[0].header
                image_of_field = candels_fields[field]+"_drz.fits"
                rms_of_field = candels_fields[field]+"_rms.fits"
                #print count, field, image_of_choosen_field

                #-----------------take values from the header-------------------
                exptime = seg_hdr['EXPTIME']
                x_naxis = seg_hdr['NAXIS1']
                y_naxis = seg_hdr['NAXIS2']

                #-----------------------------read catalogs---------------------
                #it is possible to read SExtractor tables with the astropy.io.ascii module:
                #from astropy.io import ascii
                #data=ascii.read("/home/sreis/simulations/goodsn.cat")
                field_cat = field+".cat"
                number,xx,yy,mag_aper,flux_aper,fluxerr_aper,flux_auto,fluxerr_auto,ra,dec,fwhm,sm_sex,xmin,xmax,ymin,ymax,ellip,theta = np.loadtxt(field_cat, unpack = True)

                #---------------parameters to create artificial galaxy----------
                gal_zp  = zp_imaging - 2.5 * np.log10 (exptime)
                #---galaxy position---
                # x_pix = np.random.randint(200,x_naxis-199) #random.randint return random INTEGERS from low (inclusive) to high (exclusive); I use 200 for low value to avoid edges
                # y_pix = np.random.randint(200,y_naxis-199)
                x_pix = np.random.uniform(200.,x_naxis-199.)
                y_pix = np.random.uniform(200.,y_naxis-199.)
                xmin = x_pix - size_of_fit/2.
                xmax = x_pix + size_of_fit/2.
                ymin = y_pix - size_of_fit/2.
                ymax = y_pix + size_of_fit/2.
                pa_value = np.random.uniform(-90,91) #any value within the given interval is equally likely to be drawn by uniform

                #--------------------create galfit script----------------------
                if seg_data[y_pix,x_pix] == 0:
                    file = open(str(count)+"_artificial_galaxy.script", "w")
                    file.write("# IMAGE PARAMETERS\n")
                    file.write(" A) "+image_of_field+"  # Input Data image (FITS file)\n")
                    file.write(" B) "+str(count)"_artificial_galaxy.fits  # Name for the output image\n")
                    file.write(" C) "+rms_of_field+"  # Noise image name (made from data if blank or 'none')\n")
                    file.write(" D) "+images_path+"gs_deep_f160w_v0.5_psf.fits  # Input PSF image and (optional) diffusion kernel\n")
                    file.write(" E) 1  # PSF oversampling factor relative to data\n")
                    file.write(" F) none  # Pixel mask (ASCII file or FITS file with non-0 values)\n")
                    file.write(" G) none  # Parameter constraint file (ASCII)\n")
                    file.write(" H) "+str(xmin)+" "+str(xmax)+" "+str(ymin)+" "+str(ymax)+" # Image region to fit (xmin xmax ymin ymax)\n")
                    file.write(" I) "+str(size_of_fit)+" "+str(size_of_fit)+"  # Size of convolution box (x y)\n")
                    file.write(" J) "+str(gal_zp)+"  # Magnitude photometric zeropoint\n")
                    file.write(" K) "+str(pix_scale)+" "+str(pix_scale)+"  # Plate scale (dx dy)\n")
                    file.write(" O) regular  # Display type (regular, curses, both)\n")
                    file.write(" P) 1  # Create ouput only? (1=yes; 0=optimize)\n") #item P set to 1: GALFIT will create a model image based on your input parameters and immediately quit
                    file.write(" S) 0  # Modify/create objects interactively?\n")
                    file.write("\n")
                    file.write("\n")
                    file.write("# Component number: 1\n")
                    file.write(" 0) sky                    #  Component type\n")
                    file.write(" 1) 0.00      1       #  Sky background at center of fitting region [ADUs]\n")
                    file.write(" 2) 0.000     1       #  dsky/dx (sky gradient in x)     [ADUs/pix]\n")
                    file.write(" 3) 0.000     1       #  dsky/dy (sky gradient in y)     [ADUs/pix]\n")
                    file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
                    file.write("\n")
                    #-----singlefit-----
                    file.write("# Component number: 2 - galaxy\n")
                    file.write(" 0) sersic                 #  Component type\n")
                    file.write(" 1) "+str(x_pix)+" "+str(y_pix)+" 1 1  #  Position x, y\n")
                    file.write(" 3) "+str(mag_ele)+"     1          #  Integrated magnitude\n")
                    file.write(" 4) "+str(re_ele)+"      1          #  R_e (effective radius)   [pix]\n")
                    file.write(" 5) "+str(nn_ele)+"           1          #  Sersic index n (de Vaucouleurs n=4)\n")
                    file.write(" 9) "+str(ar_ele)+"      1          #  Axis ratio (b/a)  \n")
                    file.write("10) "+str(pa_value)+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
                    file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
                    file.write("\n")
                    file.close()
