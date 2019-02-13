from arcpy.sa import *
import arcpy
import os

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
heatSPDlog300 = os.path.join(PSgdb, 'heatOSMAADTlog50')
heatbustransitlog300 = os.path.join(transitgdb, 'heatbustransitlog300')
heat_binglog300 = os.path.join(rootdir, 'results/bing/binglog300')
nlcd_imp_ps = os.path.join(rootdir, 'results/nlcd_imp_ps')

#Load models
mods = ro.r['readRDS'](os.path.join(moddir,'fieldXRFmodels.rds'))
modscale = ro.r['readRDS'](os.path.join(moddir, 'fieldXRFmodels_scaling.rds'))
print(base.names(mods)) #Get model names

#Predict pollution index
pmod = mods.rx('pollutmod')[0]
pmodc = stats.coef(pmod)
print(pmodc)


pollutind = Int(pmodc.rx('(Intercept)')[0] +
    pmodc.rx('heat_binglog300')[0] * 100*Float(arcpy.sa.Raster(heat_binglog300))/modscale.rx('heat_binglog300')[0][0] +
    pmodc.rx('heatbustransitlog300')[0] * 100*Float(Con(IsNull(Raster(heatbustransitlog300)), 0 , Raster(heatbustransitlog300)))
                /modscale.rx('heatbustransitlog300')[0][0]+
    pmodc.rx('heatOSMAADTlog50')[0] * 100*Float(Con(IsNull(Raster(heatOSMAADTlog50)), 0, Raster(heatOSMAADTlog50)))
                /modscale.rx('heatOSMAADTlog50')[0][0] +
    pmodc.rx('heat_binglog300')[0] * 100*Float(Raster(heat_binglog300))/modscale.rx('heat_binglog300')[0][0] +
    pmodc.rx('nlcd_imp_ps')[0] * Float(Con(IsNull(Raster(nlcd_imp_ps)), 0, Raster(nlcd_imp_ps))) + 0.5)
pollutind.save(os.path.join(predgdb, 'pollutind'))
