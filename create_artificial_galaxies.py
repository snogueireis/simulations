"""
Program to create artificial galaxies in random fields and parameters

Input:
- Science image
- Weight rms image
- Segmentation image from SExtractor
- PSF image
Output:
- The galaxy stamp and its correspondent rms image
- A catalog containing all the input values for creating the mock images

By Sandra Nogueira dos Reis, July 2017 @ IA - Lisbon
Updated: March 2019 @ IA - Lisbon
"""

from astropy.io import fits
from astropy.table import Table
import montage_wrapper as montage
import numpy as np
import os
import pdb #for debugging purposes only

# ------------------------------ paths -----------------------------------------
path_img = "../../../images_candels/"
path_rms = "../../../images_candels/"
path_seg = "../../../images_candels/fields_segmentation/"
path_psf = "../psfs/"

# ----------------------------- variables --------------------------------------
name_img = ["hlsp_candels_hst_wfc3_gn-tot-60mas_f160w_v1.0_drz.fits",  # GOODS-N
            "hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits",       # COSMOS
            "hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits",       # UDS
            "hlsp_candels_hst_wfc3_egs-tot-60mas_f160w_v1.0_drz.fits"] # EGS

name_rms = ["hlsp_candels_hst_wfc3_gn-tot-60mas_f160w_v1.0_rms.fits",  # GOODS-N
            "hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_rms.fits",       # COSMOS
            "hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_rms.fits",       # UDS
            "hlsp_candels_hst_wfc3_egs-tot-60mas_f160w_v1.0_rms.fits"] # EGS

name_seg = ["goodsn_segmentation.fits", # GOODS-N
            "cosmos_segmentation.fits", # COSMOS
            "uds_segmentation.fits",    # UDS
            "egs_segmentation.fits"]    # EGS

name_psf = "new_psf_hband.fits"

exptimes   = {'goodsn': 284527.6, 'cosmos': 130415.1, 'uds': 135114.9, 'aegis': 159257.8}
x_naxis    = {'goodsn':    20480, 'cosmos':    14000, 'uds':    30720, 'aegis':    40800}
y_naxis    = {'goodsn':    20480, 'cosmos':    36000, 'uds':    12800, 'aegis':    12600}
zp_imaging = {'goodsn':  25.9463, 'cosmos':  25.9463, 'uds':    25.96, 'aegis':    25.96} # taken from the official CANDELS readme files

pix_scale   = 0.06            # [arcsec/pix]
size_of_fit = 300             # [pix]
side        = size_of_fit - 1 # [pix]

mag_min = 15.5          # [mag]
mag_max = 19            # [mag]
re_min  = 0.3/pix_scale # [arcsec]/[arcsec pix^-1] = [pix] => must be in [pix] for GALFIT
re_max  = 3./pix_scale  # [pix]
nn_min  = 1.
nn_max  = 6.
ar_min  = 0.3
ar_max  = 1.

number_of_simulated_gals = 4500                         # because 7*9*10*7 (steps for each parameter) = 4410, so I round it up to 4500
number_of_gals_per_field = number_of_simulated_gals / 4 # divided by 4 because I have 4 fields
# ------------------------------------------------------------------------------


jj = 0
gal_id    = []
xx_center = []
yy_center = []
mag_in    = []
re_in     = []
nn_in     = []
ar_in     = []
pa_in     = []
# ----- for each mosaic to cut -----
for ii in range(len(name_img)):
    # --- detecting the mosaic's patch of the sky ---
    if name_img[ii].find("gn")    != -1:
        field = "goodsn"
    elif name_img[ii].find("cos") != -1:
        field = "cosmos"
    elif name_img[ii].find("uds") != -1:
        field = "uds"
    elif name_img[ii].find("egs") != -1:
        field = "aegis"

    # ----- obtaining the galaxy coordinates in pixels -----
    image_of_field = path_img + name_img[ii]
    image_of_rms   = path_img + name_rms[ii]
    field_img      = fits.open (image_of_field)
    seg_img        = fits.open (path_seg + name_seg[ii])
    # print "field image ", name_img[ii], ", segmentation image ", name_seg[ii]

    # ----- parameters to create artificial galaxy -----
    for kk in range(number_of_gals_per_field):
        gal     = jj + 1                                                  # because I don't want a galaxy with name zero
        gal_zp  = zp_imaging[field] - 2.5 * np.log10 ( exptimes[field] )
        mag     = np.random.uniform (mag_min, mag_max)                    # any value within the given interval is equally likely to be drawn by uniform
        re      = np.random.uniform (re_min,   re_max)
        nn      = np.random.uniform (nn_min,   nn_max)
        ar      = np.random.uniform (ar_min,   ar_max)
        pa      = np.random.uniform (-90,          91)

        # --- galaxy position ---
        x_pix = np.random.uniform (size_of_fit, x_naxis[field] - size_of_fit); x_pix = int(x_pix)
        y_pix = np.random.uniform (size_of_fit, y_naxis[field] - size_of_fit); y_pix = int(y_pix)

        # --- galaxy stamp ---
        xmin = int(x_pix - size_of_fit/2.)
        xmax = int(x_pix + size_of_fit/2.) - 1
        ymin = int(y_pix - size_of_fit/2.)
        ymax = int(y_pix + size_of_fit/2.) - 1
        stamp = field_img[0].data[ ymin:ymax+1, xmin:xmax+1 ]

        # --- if it is a position of a sextracted object or a zero value pixel, choose another position ---
        while not stamp.all() and (seg_img[0].data[y_pix, x_pix] != 0 or field_img[0].data[y_pix, x_pix] == 0):
            x_pix = np.random.uniform (size_of_fit, x_naxis[field] - size_of_fit); x_pix = int(x_pix)
            y_pix = np.random.uniform (size_of_fit, y_naxis[field] - size_of_fit); y_pix = int(y_pix)
            xmin  = int(x_pix - size_of_fit/2.)
            xmax  = int(x_pix + size_of_fit/2.) - 1
            ymin  = int(y_pix - size_of_fit/2.)
            ymax  = int(y_pix + size_of_fit/2.) - 1
            stamp = field_img[0].data[ ymin:ymax+1, xmin:xmax+1 ]

        # ----- create stamp with original image (without model galaxy) -----
        xx = xmin + 1
        yy = ymin + 1
        montage.mSubimage_pix (image_of_field, str(gal)+"_original_stamp.fits", xx, yy, side)
        montage.mSubimage_pix (image_of_rms, str(gal)+"_rms.fits", xx, yy, side)


        # ----- write galfit script to create model galaxy -----
        file = open(str(gal)+"_artificial_galaxy.script", "w")
        file.write("# IMAGE PARAMETERS\n")
        file.write(" A) "+str(gal)+"_original_stamp.fit  # Input Data image (FITS file)\n")
        file.write(" B) "+str(gal)+"_artificial_galaxy.fits  # Name for the output image\n")
        file.write(" C) "+str(gal)+"_rms.fits  # Noise image name (made from data if blank or 'none')\n")
        file.write(" D) "+path_psf+name_psf+" # Input PSF image and (optional) diffusion kernel\n")
        file.write(" E) 1  # PSF oversampling factor relative to data\n")
        file.write(" F) none  # Pixel mask (ASCII file or FITS file with non-0 values)\n")
        file.write(" G) none  # Parameter constraint file (ASCII)\n")
        file.write(" H) "+str(xmin)+" "+str(xmax)+" "+str(ymin)+" "+str(ymax)+" # Image region to fit (xmin xmax ymin ymax)\n")
        file.write(" I) "+str(size_of_fit)+" "+str(size_of_fit)+"  # Size of convolution box (x y)\n")
        file.write(" J) "+str(gal_zp)+"  # Magnitude photometric zeropoint\n")
        file.write(" K) "+str(pix_scale)+" "+str(pix_scale)+"  # Plate scale (dx dy)\n")
        file.write(" O) regular  # Display type (regular, curses, both)\n")
        file.write(" P) 1  # Create ouput only? (1=yes; 0=optimize)\n") # item P set to 1: GALFIT will create a model image based on your input parameters and immediately quit
        file.write(" S) 0  # Modify/create objects interactively?\n")
        file.write("\n")
        file.write("\n")

        # --- singlefit ---
        file.write("# Component number: 1 - galaxy\n")
        file.write(" 0) sersic                 #  Component type\n")
        file.write(" 1) "+str(x_pix)+" "+str(y_pix)+" 1 1  #  Position x, y\n")
        file.write(" 3) "+str(mag)+"     1          #  Integrated magnitude\n")
        file.write(" 4) "+str(re)+"      1          #  R_e (effective radius)   [pix]\n")
        file.write(" 5) "+str(nn)+"           1          #  Sersic index n (de Vaucouleurs n=4)\n")
        file.write(" 9) "+str(ar)+"      1          #  Axis ratio (b/a)  \n")
        file.write("10) "+str(pa)+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
        file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
        file.write("\n")
        file.close()


        # ----- run galfit to create model galaxy -----
        os.system("galfit "+str(gal)+"_artificial_galaxy.script")

        # ----- adding the model to the real image -----
        img_model  = fits.open (str(gal)+"_artificial_galaxy.fits") # open model galaxy I just created with GALFIT
        model_data = img_model[0].data

        # --- add model galaxy to the original image ---
        stamp_header = fits.getheader (str(gal)+"_original_stamp.fits")
        final_img    = stamp + model_data
        fits.writeto (str(gal)+"_mock_galaxy.fits", final_img, header = stamp_header, overwrite = True) # save final fits stamp, with original image from the CANDELS field + model galaxy
        jj = jj + 1

        # ----- saving everything in the catalog -----
        gal_id.append(gal)
        xx_center.append(xx)
        yy_center.append(yy)
        mag_in.append(mag)
        re_in.append(re)
        nn_in.append(nn)
        ar_in.append(ar)
        pa_in.append(pa)
    new_table = Table([gal_id,xx_center,yy_center,mag_in,re_in,nn_in,ar_in,pa_in],names=['gal','x','y','mag','re','n','ar','pa'])
new_table.write("catalog_mock_galaxies_hband.cat", formats = {'mag':'%0.4f', 're':'%0.4f', 'n':'%0.4f', 'ar':'%0.4f', 'pa':'%0.4f'}, format = "ascii.commented_header", overwrite = True)
