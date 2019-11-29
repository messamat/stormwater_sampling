from arcpy.sa import *
import arcpy
import os
import numpy as np

'''
To get rpy2 running
- download whl from https://www.lfd.uci.edu/~gohlke/pythonlibs/#rpy2 into c:\python27\arcgisx6410.7\Scripts
- terminal: cd c:\python27\arcgisx6410.6=7\Scripts
- terminal: pip install rpy2-2.7.8-cp27-none-win_amd64.whl
- in R: normalizePath(R.home()) - add it to R_HOME system variable in environment variables
- in R:  Sys.info()['user'] - add it to R_USER system variable in environment variables
'''
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
base = importr('base')
stats = importr('stats')

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False
arcpy.env.cellSize = "MINOF"

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'
moddir = os.path.join(rootdir, 'results/data_modeling')
Binggdb = os.path.join(rootdir, 'results/bing/postprocess.gdb')
transitgdb = os.path.join(rootdir, 'results/transit.gdb')
Seattlegdb = os.path.join(rootdir, 'results/Seattle_sampling.gdb')
pollutgdb = os.path.join(rootdir,'results/pollution_variables.gdb')
pollutPSgdb = os.path.join(rootdir, 'results/pollution_variables_PS.gdb')
NLCDreclass = os.path.join(rootdir, 'results/LU.gdb/NLCD_reclass_final')

#Create new gdb for predictions
predgdb = os.path.join(rootdir,'results/PSpredictions.gdb')
if arcpy.Exists(predgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'PSpredictions.gdb')

#Define variables
heatsubAADTlog100 = os.path.join(pollutPSgdb, 'heatPSAADTlog100')
heatbustransitpow100_1 = os.path.join(pollutgdb, 'heatbustransitpow100_1')
heatbustransitlog200 = os.path.join(pollutgdb, 'heatbustransitlog200')
nlcd_imp_ps_mean = os.path.join(rootdir, 'results/nlcd_imp_ps_mean.tif')
predzn36 = os.path.join(predgdb, 'predzn36')
predcu19 = os.path.join(predgdb, 'predcu19')
predpi30 = os.path.join(predgdb, 'predpi30')

#Load models
mods = ro.r['readRDS'](os.path.join(moddir,'fieldXRFmodels.rds'))
modscale = ro.r['readRDS'](os.path.join(moddir, 'fieldXRFmodels_scaling.rds'))
print(base.names(mods)) #Get model names

#Predict Zn
znmod = mods.rx('logZnmod')[0]
znmodc = stats.coef(znmod)
print(znmodc)

#Model 36
predzn = Int((100*
             Exp(
                 znmodc.rx('(Intercept)')[0] +

                 (znmodc.rx('heatsubAADTlog100frt')[0] *
                  ((100*Float(Con(IsNull(Raster(heatsubAADTlog100)), 0, Raster(heatsubAADTlog100)))
                  /modscale.rx('heatsubAADTlog100')[0][0])**(0.25))) +

                 (znmodc.rx('nlcd_imp_ps_mean')[0] *
                  Float(Con(IsNull(Raster(nlcd_imp_ps_mean)), 0, Raster(nlcd_imp_ps_mean))))
             )) + 0.5
             )
print('Saving Zinc predictions...')
predzn.save(predzn36)

#Export attribute table
predzn36_land = arcpy.sa.Con(Raster(NLCDreclass) != 11, Raster(predzn36))
predzn36_land.save(predzn36 + 'land')
arcpy.CopyRows_management(predzn36 + 'land', out_table = os.path.join(predgdb, 'predzn36_tab'))

#Predict Cu
cumod = mods.rx('Cumod')[0]
cumodc = stats.coef(cumod)
print(cumodc)

predcu = Int((100*
              (
                      cumodc.rx('(Intercept)')[0] +

                      (cumodc.rx('heatbustransitlog200sqrt')[0] * \
                       ((100*Float(Con(IsNull(Raster(heatbustransitlog200)), 0, Raster(heatbustransitlog200)))
                         /modscale.rx('heatbustransitlog200')[0][0])**(0.5))) +

                      (cumodc.rx('nlcd_imp_ps_mean')[0] *
                       Float(Con(IsNull(Raster(nlcd_imp_ps_mean)), 0, Raster(nlcd_imp_ps_mean))))
              )) + 0.5
             )

print('Saving Copper predictions...')
predcu.save(predcu19)

#Predict pollution index
pimod = mods.rx('PImod')[0]
pimodc = stats.coef(pimod)
print(pimodc)

predpi = Int((100*
              (
                      pimodc.rx('(Intercept)')[0] +

                      (pimodc.rx('heatsubAADTlog100thd')[0] *
                       ((100*Float(Con(IsNull(Raster(heatsubAADTlog100)), 0, Raster(heatsubAADTlog100)))
                         /modscale.rx('heatsubAADTlog100')[0][0])**(1.0/3.0))) +

                      (pimodc.rx('nlcd_imp_ps_mean')[0] *
                       Float(Con(IsNull(Raster(nlcd_imp_ps_mean)), 0, Raster(nlcd_imp_ps_mean)))) +

                      (pimodc.rx('heatbustransitpow100_1')[0] * \
                       ((100 * Float(Con(IsNull(Raster(heatbustransitpow100_1)), 0, Raster(heatbustransitpow100_1)))
                         / modscale.rx('heatbustransitpow100_1')[0][0])))
              )) + 0.5
             )
print('Saving Pollution index predictions...')
predpi.save(predpi30)

XRFsites = os.path.join(pollutgdb, 'XRFsites_aea')
sitescheck = os.path.join(predgdb,'sitescheck')
arcpy.CopyFeatures_management(XRFsites, sitescheck)
arcpy.sa.ExtractMultiValuesToPoints(sitescheck, [predzn36, predcu19, predpi30], bilinear_interpolate_values='NONE')