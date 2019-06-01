from arcpy.sa import *
import arcpy
import os
import numpy as np

'''
To get rpy2 running
- download whl from https://www.lfd.uci.edu/~gohlke/pythonlibs/#rpy2 into c:\python27\arcgisx6410.6\Scripts
- terminal: cd c:\python27\arcgisx6410.6\Scripts
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

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'
moddir = os.path.join(rootdir, 'results/data_modeling')
PSgdb=os.path.join(rootdir,'results/PSOSM.gdb')
Binggdb = os.path.join(rootdir, 'results/bing/postprocess.gdb')
transitgdb = os.path.join(rootdir, 'results/transit.gdb')
Seattlegdb = os.path.join(rootdir, 'results/Seattle_sampling.gdb')
pollutgdb = os.path.join(rootdir,'results/pollution_variables.gdb')

#Create new gdb for predictions
predgdb = os.path.join(rootdir,'results/PSpredictions.gdb')
# if arcpy.Exists(predgdb):
#     print('Geodatabase already exists')
# else:
#     arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'PSpredictions.gdb')

#Define variables
heatsubAADTlog100 = os.path.join(pollutgdb, 'heatsubAADTlog100')
heatbustransitlog200 = os.path.join(transitgdb, 'heatbustransitlog300')
heat_binglog300 = os.path.join(rootdir, 'results/bing/binglog300')
heatbing1902log300 = os.path.join(rootdir, 'results/bing/heatbing1902log300.tif')
nlcd_imp_ps_mean = os.path.join(rootdir, 'results/nlcd_imp_ps_mean.tif')

#Load models
mods = ro.r['readRDS'](os.path.join(moddir,'fieldXRFmodels.rds'))
modscale = ro.r['readRDS'](os.path.join(moddir, 'fieldXRFmodels_scaling.rds'))
print(base.names(mods)) #Get model names

#Predict Zn
znmod = mods.rx('logZnmod')[0]
znmodc = stats.coef(znmod)
print(znmodc)

#Model 33
predzn = Int((100*
             Exp(
                 znmodc.rx('(Intercept)')[0] +
                 (znmodc.rx('heatbing1902log300')[0] * 100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0]) +
                 (znmodc.rx('heatsubAADTlog100frt')[0] * ((100*Float(Con(IsNull(Raster(heatsubAADTlog100)), 0, Raster(heatsubAADTlog100)))
                  /modscale.rx('heatsubAADTlog100')[0][0])**(1/4))) +
                 (znmodc.rx('nlcd_imp_ps_mean')[0] * Float(Con(IsNull(Raster(nlcd_imp_ps_mean)), 0, Raster(nlcd_imp_ps_mean))))+
                 (znmodc.rx('nlcd_imp_ps_mean:heatbing1902log300')[0] *
                  (100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0]) *
                  Float(Con(IsNull(Raster(nlcd_imp_ps_mean)), 0, Raster(nlcd_imp_ps_mean))))
             )) + 0.5
             )
predzn.save(os.path.join(predgdb, 'predzn33'))

heatbing = (znmodc.rx('heatbing1902log300')[0] * 100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0])
heatbing.save(os.path.join(predgdb, 'testbing'))
checkaadt = (znmodc.rx('heatsubAADTlog100frt')[0] * ((100*Float(Con(IsNull(Raster(heatsubAADTlog100)), 0, Raster(heatsubAADTlog100)))
                  /modscale.rx('heatsubAADTlog100')[0][0])**(1/4)))
checkaadt.save(os.path.join(predgdb, 'testaadt'))
checkinterac =                  (znmodc.rx('nlcd_imp_ps_mean:heatbing1902log300')[0] *
                  (100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0]) *
                  Float(Con(IsNull(Raster(nlcd_imp_ps_mean)), 0, Raster(nlcd_imp_ps_mean))))
checkinterac.save(os.path.join(predgdb, 'checkinterac'))

XRFsites_projUTM = os.path.join(PSgdb, 'XRFsites_projUTM')

sitescheck = os.path.join(predgdb,'sitescheck')
UTM10 = arcpy.SpatialReference(26910)
arcpy.Project_management(os.path.join(rootdir, 'results/PSpredictions.gdb/XRFsites_proj'), sitescheck, UTM10)
arcpy.sa.ExtractMultiValuesToPoints(sitescheck, [os.path.join(predgdb, ras) for ras in ['testbing','checkinterac']])
arcpy.CopyRows_management(sitescheck, os.path.join(rootdir, 'results/checkznint.csv'))

arcpy.CopyFeatures_management(os.path.join(rootdir, 'results/PSpredictions.gdb/XRFsites_proj'), sitescheck+'notUTM')
arcpy.sa.ExtractMultiValuesToPoints(sitescheck, [heatbing1902log300, nlcd_imp_ps_mean], bilinear_interpolate_values='NONE')
arcpy.sa.ExtractMultiValuesToPoints(sitescheck+'notUTM', [heatsubAADTlog100])
arcpy.CopyRows_management(sitescheck, os.path.join(rootdir, 'results/checkznint2.csv'))
arcpy.CopyRows_management(sitescheck+'notUTM', os.path.join(rootdir, 'results/checkznint3.csv'))

#Predict pollution index
pmod = mods.rx('pollutmod')[0]
pmodc = stats.coef(pmod)
print(pmodc)

pollutind43 = Int(pmodc.rx('(Intercept)')[0] +
                  (pmodc.rx('heatbing1902log300')[0] * 100*Float(Con(IsNull(Raster(heatbing1902log300)), 0 , Raster(heatbing1902log300)))
                  /modscale.rx('heatbing1902log300')[0][0]) +
                  (pmodc.rx('heatbustransitlog300')[0] * 100*Float(Con(IsNull(Raster(heatbustransitlog300)), 0 , Raster(heatbustransitlog300)))
                  /modscale.rx('heatbustransitlog300')[0][0]) +
                  (pmodc.rx('nlcd_imp_ps')[0] * Float(Con(IsNull(Raster(nlcd_imp_ps)), 0, Raster(nlcd_imp_ps)))) +
                  0.5
                  )
pollutind43.save(os.path.join(predgdb, 'pollutind43'))
#arcpy.Delete_management(os.path.join(predgdb, 'pollutind43'))




#If still doesn't work, make sure that scalings are right