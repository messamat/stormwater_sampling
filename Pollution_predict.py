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
rootdir = 'C:/Mathis/ICSL/stormwater'
moddir = os.path.join(rootdir, 'results/data_modeling')
PSgdb=os.path.join(rootdir,'results/PSOSM.gdb')
Binggdb = os.path.join(rootdir, 'results/bing/postprocess.gdb')
transitgdb = os.path.join(rootdir, 'results/transit.gdb')
Seattlegdb = os.path.join(rootdir, 'results/Seattle_sampling.gdb')

#Create new gdb for predictions
predgdb = os.path.join(rootdir,'results/PSpredictions.gdb')
# if arcpy.Exists(predgdb):
#     print('Geodatabase already exists')
# else:
#     arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'PSpredictions.gdb')

#Define variables
heatOSMAADTlog50 = os.path.join(Seattlegdb, 'heatOSMAADTlog50')
heatOSMAADTlog100 = os.path.join(Seattlegdb, 'heatOSMAADTlog100')
heatSPDlog300 = os.path.join(PSgdb, 'heatSPDlog300')
heatbustransitlog300 = os.path.join(transitgdb, 'heatbustransitlog300')
heat_binglog300 = os.path.join(rootdir, 'results/bing/binglog300')
heatbing1902log300 = os.path.join(rootdir, 'results/bing/heatbing1902log300.tif')
nlcd_imp_ps = os.path.join(rootdir, 'results/nlcd_imp_ps')
heatOSMgradientlog300 = os.path.join(PSgdb, 'heatOSMgradientlog300')

#Load models
mods = ro.r['readRDS'](os.path.join(moddir,'fieldXRFmodels.rds'))
modscale = ro.r['readRDS'](os.path.join(moddir, 'fieldXRFmodels_scaling.rds'))
print(base.names(mods)) #Get model names

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



# pollutind = Int(pmodc.rx('(Intercept)')[0] +
#     pmodc.rx('heat_binglog300')[0] * 100*Float(arcpy.sa.Raster(heat_binglog300))/modscale.rx('heat_binglog300')[0][0] +
#     pmodc.rx('heatbustransitlog300')[0] * 100*Float(Con(IsNull(Raster(heatbustransitlog300)), 0 , Raster(heatbustransitlog300)))
#                 /modscale.rx('heatbustransitlog300')[0][0]+
#     pmodc.rx('heatOSMAADTlog50')[0] * 100*Float(Con(IsNull(Raster(heatOSMAADTlog50)), 0, Raster(heatOSMAADTlog50)))
#                 /modscale.rx('heatOSMAADTlog50')[0][0] +
#     pmodc.rx('heat_binglog300')[0] * 100*Float(Raster(heat_binglog300))/modscale.rx('heat_binglog300')[0][0] +
#     pmodc.rx('nlcd_imp_ps')[0] * Float(Con(IsNull(Raster(nlcd_imp_ps)), 0, Raster(nlcd_imp_ps))) + 0.5)
# pollutind.save(os.path.join(predgdb, 'pollutind'))


#Predict Zn
znmod = mods.rx('logZnmod')[0]
znmodc = stats.coef(znmod)
print(znmodc)

#Model 30
# predzn = Int(100*
#              Exp(
#                  znmodc.rx('(Intercept)')[0] +
#                  (znmodc.rx('heatbustransitlog300')[0] * 100*Float(Con(IsNull(Raster(heatbustransitlog300)), 0 , Raster(heatbustransitlog300)))
#                   /modscale.rx('heatbustransitlog300')[0][0]) +
#                  (znmodc.rx('heatbing1902log300')[0] * 100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0]) +
#                  (znmodc.rx('heatOSMAADTlog100')[0] * 100*Float(Con(IsNull(Raster(heatOSMAADTlog100)), 0, Raster(heatOSMAADTlog100)))
#                   /modscale.rx('heatOSMAADTlog100')[0][0]) +
#                  (znmodc.rx('heatSPDlog300')[0] * 100 * Float(Con(IsNull(Raster(heatSPDlog300)), 0, Raster(heatSPDlog300)))
#                   /modscale.rx('heatSPDlog300')[0][0]) +
#                  (znmodc.rx('heatSPDlog300:heatOSMAADTlog100')[0] *
#                   (100*Float(Con(IsNull(Raster(heatSPDlog300)), 0, Raster(heatSPDlog300)))/modscale.rx('heatSPDlog300')[0][0]) *
#                   (100*Float(Con(IsNull(Raster(heatOSMAADTlog100)), 0, Raster(heatOSMAADTlog100)))/modscale.rx('heatOSMAADTlog100')[0][0]))
#              )
#              )
#
# #Model 11
# predzn = Int(100*
#              Exp(
#                  znmodc.rx('(Intercept)')[0] +
#                  (znmodc.rx('heatbustransitlog300')[0] * 100*Float(Con(IsNull(Raster(heatbustransitlog300)), 0 , Raster(heatbustransitlog300)))
#                   /modscale.rx('heatbustransitlog300')[0][0]) +
#                  (znmodc.rx('heatbing1902log300')[0] * 100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0]) +
#                  (znmodc.rx('heatSPDlog300')[0] * 100 * Float(Con(IsNull(Raster(heatSPDlog300)), 0, Raster(heatSPDlog300)))
#                   /modscale.rx('heatSPDlog300')[0][0]) +
#                  (znmodc.rx('heatSPDlog300:heatbing1902log300')[0] *
#                   (100*Float(Con(IsNull(Raster(heatSPDlog300)), 0, Raster(heatSPDlog300)))/modscale.rx('heatSPDlog300')[0][0]) *
#                   (100*Float(Con(IsNull(Raster(heatbing1902log300)), 0, Raster(heatbing1902log300)))/modscale.rx('heatbing1902log300')[0][0]))
#              )
#              )

# #Model 12
# predzn12 = Int(100*
#              Exp(
#                  znmodc.rx('(Intercept)')[0] +
#                  (znmodc.rx('heatbustransitlog300')[0] * 100*Float(Con(IsNull(Raster(heatbustransitlog300)), 0 , Raster(heatbustransitlog300)))
#                   /modscale.rx('heatbustransitlog300')[0][0]) +
#                  (znmodc.rx('heatbing1902log300')[0] * 100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0]) +
#                  (znmodc.rx('heatSPDlog300')[0] * 100 * Float(Con(IsNull(Raster(heatSPDlog300)), 0, Raster(heatSPDlog300)))
#                   /modscale.rx('heatSPDlog300')[0][0]) +
#                  (znmodc.rx('I(heatSPDlog300^2)')[0] *
#                   Power((100*Float(Con(IsNull(Raster(heatbing1902log300)), 0, Raster(heatbing1902log300)))
#                          /modscale.rx('heatbing1902log300')[0][0]),
#                         2)
#                   )
#              )
#              )
# predzn12.save(os.path.join(predgdb, 'predzn12'))

#Model 28
# predzn28 = Int(100*
#                Exp(
#                    znmodc.rx('(Intercept)')[0] +
#                    (znmodc.rx('heatbustransitlog300')[0] * 100*Float(Con(IsNull(Raster(heatbustransitlog300)), 0 , Raster(heatbustransitlog300)))
#                     /modscale.rx('heatbustransitlog300')[0][0]) +
#                    (znmodc.rx('heatbing1902log300')[0] * 100*Float(arcpy.sa.Raster(heatbing1902log300))/modscale.rx('heatbing1902log300')[0][0]) +
#                    (znmodc.rx('heatSPDlog300')[0] * 100 * Float(Con(IsNull(Raster(heatSPDlog300)), 0, Raster(heatSPDlog300)))
#                     /modscale.rx('heatSPDlog300')[0][0]) +
#                    (znmodc.rx('heatOSMgradientlog300')[0] * 100 * Float(
#                        Con(IsNull(Raster(heatOSMgradientlog300)), 0, Raster(heatOSMgradientlog300)))
#                     / modscale.rx('heatOSMgradientlog300')[0][0]) +
#                    (znmodc.rx('heatSPDlog300:heatOSMgradientlog300')[0] *
#                     (100*Float(Con(IsNull(Raster(heatSPDlog300)), 0, Raster(heatSPDlog300)))/modscale.rx('heatSPDlog300')[0][0]) *
#                     (100*Float(Con(IsNull(Raster(heatOSMgradientlog300)), 0, Raster(heatOSMgradientlog300)))/modscale.rx('heatOSMgradientlog300')[0][0])
#                     )
#                )
#                )
# predzn28.save(os.path.join(predgdb, 'predzn28'))



#If still doesn't work, make sure that scalings are right