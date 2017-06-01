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
goodsn = images_path+"goodsn_all_wfc3_ir_f160w_060mas_v1.0_drz.fits"
cosmos = images_path+"hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits"
uds = images_path+"hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits"
egs = images_path+"egs_all_wfc3_ir_f160w_060mas_v1.0_drz.fits"
candels_fields = [goodsn,cosmos,uds,egs]
#sample    = "galaxies.txt"
#gal,field = np.loadtxt(sample,dtype=str,unpack=True)
#----------------------------------------------------------------------------

for ii in range(len(candels_fields)):
	field_image = candels_fields[ii]
	gal_cat  = candels_fields[ii]+".cat"
	seg_img  = candels_fields[ii]+"_segmentation.fits"
	rms_img  = candels_fields[ii]+"_rms.fits"

	if field[ii] == "GOODS-N":
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+candels_cat+' -CHECKIMAGE_TYPE SEGMENTATION'+' -CHECKIMAGE_NAME '+seg_img+' -MAG_ZEROPOINT '+str(mag_zp))
	elif field[ii] == "COSMOS":
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+gal_cat+" -CHECKIMAGE_NAME "+seg_img+" -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE "+wht_img)
	elif field[ii] == "UDS":
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+gal_cat+" -CHECKIMAGE_NAME "+seg_img+" -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE "+wht_img)
	elif field[ii] == "AEGIS":
		os.system("sextractor "+field_image+" -c "+conf_file+" -CATALOG_NAME "+gal_cat+" -CHECKIMAGE_NAME "+seg_img+" -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE "+wht_img)
