"""
Program to create artificial galaxies in random fields and positions

By Sandra Nogueira dos Reis, June 2017 @ IA - Lisbon
"""

from astropy.io import fits
import numpy as np
import os

#--------------------------------paths-----------------------------------------
images_path    = "/home/sreis/candels_images/"
masks_path     = "/home/sreis/candels_images/masks/"

#---------------------------------constants------------------------------------
zp_imaging = 25.95
pix_scale  = 0.06

#-----------------------------read catalogs------------------------------------
sample       = catalogs_path+"samplecat_with_morphologies.txt"
gal_id,morphology,field,ra,dec,redshift,lmass= np.loadtxt(sample,dtype=str,unpack=True)
#------------------------------------------------------------------------------

for ii in range(len(gal_id)):
    gal = int(gal_id[ii])
    img_cut = fits.open(images_path+str(gal)+"_cut.fits")
    hdr     = img_cut[0].header #read the header of the cutted image
    exptime = hdr['EXPTIME']
    gal_zp  = zp_imaging - 2.5 * np.log10 (exptime)

    #---------take the values from the header of the model from galfit---------
##vou pegar numa posição aleatória onde o valor do segmentation map seja igual a zero e criar aí uma galáxia

    img = fits.open(doublefit_path+str(gal)+"_doublefit_fixedn1.fits")
    model = img[2]

    x_pix               = model.header['NAXIS1']
    y_pix               = model.header['NAXIS2']

    gal_x_disk_model    = model.header['2_XC']
    gal_x_disk          = float(gal_x_disk_model.split()[0])
    gal_y_disk_model    = model.header['2_YC']
    gal_y_disk          = float(gal_y_disk_model.split()[0])
    gal_mag_disk_model  = model.header['2_MAG']
    gal_mag_disk        = float(gal_mag_disk_model.split()[0])
    gal_re_disk_model   = model.header['2_RE']
    gal_re_disk         = float(gal_re_disk_model.split()[0])
    gal_n_disk_model    = model.header['2_N']
    gal_n_disk          = float(1.0000)
    gal_ar_disk_model   = model.header['2_AR']
    gal_ar_disk         = float(gal_ar_disk_model.split()[0])
    gal_pa_disk_model   = model.header['2_PA']
    gal_pa_disk         = float(gal_pa_disk_model.split()[0])

    gal_x_bulge_model   = model.header['3_XC']
    gal_x_bulge         = float(gal_x_bulge_model.split()[0])
    gal_y_bulge_model   = model.header['3_YC']
    gal_y_bulge         = float(gal_y_bulge_model.split()[0])
    gal_mag_bulge_model = model.header['3_MAG']
    gal_mag_bulge       = float(gal_mag_bulge_model.split()[0])
    gal_re_bulge_model  = model.header['3_RE']
    gal_re_bulge        = float(gal_re_bulge_model.split()[0])
    gal_n_bulge_model   = model.header['3_N']
    gal_n_bulge         = float(gal_n_bulge_model.split()[0])
    gal_ar_bulge_model  = model.header['3_AR']
    gal_ar_bulge        = float(gal_ar_bulge_model.split()[0])
    gal_pa_bulge_model  = model.header['3_PA']
    gal_pa_bulge        = float(gal_pa_bulge_model.split()[0])

    #-----------------------------write galfit script--------------------------
    file = open(str(gal)+"_model.script", "w")
    file.write("# IMAGE PARAMETERS\n")
    file.write(" A) "+images_path+str(gal)+"_cut.fits  # Input Data image (FITS file)\n")
    file.write(" B) "+str(gal)+"_model.fits  # Name for the output image\n")
    file.write(" C) "+images_path+str(gal)+"_rms.fits  # Noise image name (made from data if blank or 'none')\n")
    file.write(" D) "+images_path+"gs_deep_f160w_v0.5_psf.fits  # Input PSF image and (optional) diffusion kernel\n")
    file.write(" E) 1  # PSF oversampling factor relative to data\n")
    file.write(" F) "+masks_path+str(gal)+"_mask.fits  # Pixel mask (ASCII file or FITS file with non-0 values)\n")
    file.write(" G) "+doublefit_path+str(gal)+"_constraints.txt  # Parameter constraint file (ASCII)\n")
    file.write(" H) 1         "+str(x_pix)+ "  1          "+str(y_pix)+ "  # Image region to fit (xmin xmax ymin ymax)\n")
    file.write(" I) "+str(x_pix)+"      "+str(y_pix)+"  # Size of convolution box (x y)\n")
    file.write(" J) "+str(gal_zp)+"  # Magnitude photometric zeropoint\n")
    file.write(" K) "+str(pix_scale)+"     "+str(pix_scale)+"  # Plate scale (dx dy)\n")
    file.write(" O) regular  # Display type (regular, curses, both)\n")
    file.write(" P) 3  # Create ouput only? (1=yes; 0=optimize)\n")
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
    #-----disk-----
    file.write("# Component number: 2 - disk\n")
    file.write(" 0) sersic                 #  Component type\n")
    file.write(" 1) "+str(gal_x_disk)+" "+str(gal_y_disk)+" 1 1  #  Position x, y\n")
    file.write(" 3) "+str(gal_mag_disk)+"     1          #  Integrated magnitude\n")
    file.write(" 4) "+str(gal_re_disk)+"      1          #  R_e (effective radius)   [pix]\n")
    file.write(" 5) "+str(gal_n_disk)+"           1          #  Sersic index n (de Vaucouleurs n=4)\n")
    file.write(" 9) "+str(gal_ar_disk)+"      1          #  Axis ratio (b/a)  \n")
    file.write("10) "+str(gal_pa_disk)+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
    file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
    file.write("\n")
    #-----bulge-----
    file.write("# Component number: 3 - bulge\n")
    file.write(" 0) sersic                 #  Component type\n")
    file.write(" 1) "+str(gal_x_bulge)+" "+str(gal_y_bulge)+" 1 1  #  Position x, y\n")
    file.write(" 3) "+str(gal_mag_bulge)+"     1          #  Integrated magnitude\n")
    file.write(" 4) "+str(gal_re_bulge)+"      1          #  R_e (effective radius)   [pix]\n")
    file.write(" 5) "+str(gal_n_bulge)+"           1          #  Sersic index n (de Vaucouleurs n=4)\n")
    file.write(" 9) "+str(gal_ar_bulge)+"      1          #  Axis ratio (b/a)  \n")
    file.write("10) "+str(gal_pa_bulge)+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
    file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
    file.write("\n")
    #-----bar-----
    if gal in [13942, 37194]:
        gal_x_bar_model   = model.header['4_XC']
        gal_x_bar         = float(gal_x_bar_model.split()[0])
        gal_y_bar_model   = model.header['4_YC']
        gal_y_bar        = float(gal_y_bar_model.split()[0])
        gal_mag_bar_model = model.header['4_MAG']
        gal_mag_bar       = float(gal_mag_bar_model.split()[0])
        gal_re_bar_model  = model.header['4_RE']
        gal_re_bar        = float(gal_re_bar_model.split()[0])
        gal_n_bar_model   = model.header['4_N']
        gal_n_bar         = float(gal_n_bar_model.split()[0])
        gal_ar_bar_model  = model.header['4_AR']
        gal_ar_bar        = float(gal_ar_bar_model.split()[0])
        gal_pa_bar_model  = model.header['4_PA']
        gal_pa_bar        = float(gal_pa_bar_model.split()[0])

        file.write("# Component number: 4 - bar\n")
        file.write(" 0) sersic                 #  Component type\n")
        file.write(" 1) "+str(gal_x_bar)+" "+str(gal_y_bar)+" 1 1  #  Position x, y\n")
        file.write(" 3) "+str(gal_mag_bar)+"     1          #  Integrated magnitude\n")
        file.write(" 4) "+str(gal_re_bar)+"      1          #  R_e (effective radius)   [pix]\n")
        file.write(" 5) "+str(gal_n_bar)+"           1          #  Sersic index n (de Vaucouleurs n=4)\n")
        file.write(" 9) "+str(gal_ar_bar)+"      1          #  Axis ratio (b/a)  \n")
        file.write("10) "+str(gal_pa_bar)+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
        file.write("C0) 0.5      0          #  Diskyness(-)/Boxyness(+)\n")
        file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
        file.write("\n")

    file.close()
