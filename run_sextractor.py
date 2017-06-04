"""
Run SExtractor in CANDELS fields
Sandra Nogueira dos Reis, June 2017 @ Lisbon
"""

import numpy as np
import os

#-----------------------------------paths------------------------------------
images_path = "/home/sandra/candels_images/"

#--------------------------------constants-----------------------------------
conf_file = "wfc3_H.sex"
mag_zp = 25.95
goodsn = "goodsn_all_wfc3_ir_f160w_060mas_v1.0"
cosmos = "hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0"
uds = "hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0"
egs = "egs_all_wfc3_ir_f160w_060mas_v1.0"
candels_fields = [goodsn,cosmos,uds,egs]
#sample    = "galaxies.txt"
#gal,field = np.loadtxt(sample,dtype=str,unpack=True)
#----------------------------------------------------------------------------

for ii in range(len(candels_fields)):
	field_image = images_path+candels_fields[ii]+"_drz.fits"
	rms_img     = images_path+candels_fields[ii]+"_rms.fits"
	field_cat   = candels_fields[ii]+".cat"
	seg_img     = candels_fields[ii]+"_segmentation.fits"

	if candels_fields[ii] == goodsn:
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+field_cat+' -CHECKIMAGE_TYPE SEGMENTATION'+' -CHECKIMAGE_NAME '+seg_img+' -MAG_ZEROPOINT '+str(mag_zp)+" -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE "+rms_img)
	elif candels_fields[ii] == cosmos:
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+field_cat+' -CHECKIMAGE_TYPE SEGMENTATION'+' -CHECKIMAGE_NAME '+seg_img+' -MAG_ZEROPOINT '+str(mag_zp)+" -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE "+rms_img)
	elif candels_fields[ii] == uds:
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+field_cat+' -CHECKIMAGE_TYPE SEGMENTATION'+' -CHECKIMAGE_NAME '+seg_img+' -MAG_ZEROPOINT '+str(mag_zp)+" -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE "+rms_img)
	elif candels_fields[ii] == egs:
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+field_cat+' -CHECKIMAGE_TYPE SEGMENTATION'+' -CHECKIMAGE_NAME '+seg_img+' -MAG_ZEROPOINT '+str(mag_zp)+" -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE "+rms_img)
