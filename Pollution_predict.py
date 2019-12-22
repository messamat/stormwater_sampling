from arcpy.sa import *
import arcpy
import os
import numpy as np

'''
To get rpy2 running
- download whl from https://www.lfd.uci.edu/~gohlke/pythonlibs/#rpy2 into c:\python27\arcgisx6410.7\Scripts
- terminal: cd c:\python27\\arcgisx6410.6=7\Scripts
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
pollutfieldclean_cast = os.path.join(rootdir, 'results/data_modeling/pollutfieldclean_cast.dbf')
template_ras = os.path.join(rootdir, 'results/bing/heatbing1902log300proj.tif')
restemplate = arcpy.GetRasterProperties_management(template_ras, 'CELLSIZEX')
PSdissolve = os.path.join(rootdir, 'results/PSwtshd_dissolve.shp')

#Create new gdb for predictions
predgdb = os.path.join(rootdir,'results/PSpredictions.gdb')
if arcpy.Exists(predgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'PSpredictions.gdb')

#Define variables
PSdissolveras = os.path.join(predgdb, 'PSwtshd_dissolve')
heatsubAADTlog100 = os.path.join(pollutPSgdb, 'heatPSAADTlog100')
heatbustransitpow100_1 = os.path.join(pollutPSgdb, 'heatPSbustransitpow100_1')
heatbustransitlog200 = os.path.join(pollutPSgdb, 'heatPSbustransitlog200')
nlcd_imp_ps_mean = os.path.join(rootdir, 'results/nlcd_imp_ps_mean.tif')
predzn36 = os.path.join(predgdb, 'predzn36')
predzn36index = os.path.join(predgdb, 'predzn36index')
predcu19 = os.path.join(predgdb, 'predcu19')
predcu19index = os.path.join(predgdb, 'predcu19index')
predpi30 = os.path.join(predgdb, 'predpi30')

#Create template raster for predictions
arcpy.env.snapRaster = template_ras
if not arcpy.Exists(PSdissolveras):
    arcpy.CalculateField_management(PSdissolve, field='Id', expression = 1, expression_type='PYTHON')
    arcpy.PolygonToRaster_conversion(PSdissolve, 'Id', PSdissolveras, cellsize = restemplate)

#Define extent
arcpy.env.extent = PSdissolveras

#Load models
mods = ro.r['readRDS'](os.path.join(moddir,'fieldXRFmodels.rds'))
modscale = ro.r['readRDS'](os.path.join(moddir, 'fieldXRFmodels_scaling.rds'))
print(base.names(mods)) #Get model names

#Predict Zn
znmod = mods.rx('logZnstandmod')[0]
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

#Standardize to 0-100
minznpred = float(arcpy.GetRasterProperties_management (predzn36, "MINIMUM").getOutput (0))
maxznpred = float(arcpy.GetRasterProperties_management (predzn36, "MAXIMUM").getOutput (0))
predznindex = 100*(Float(Raster(predzn36))- minznpred)/(maxznpred - minznpred)
predznindex.save(predzn36index)

#Export attribute table
predzn36_land = arcpy.sa.Con((Raster(NLCDreclass) != 11) & (Raster(PSdissolveras) == 1), Raster(predzn36))
predzn36_land.save(predzn36 + 'land')
arcpy.CopyRows_management(predzn36 + 'land', out_table = os.path.join(predgdb, 'predzn36_tab'))

#Predict Cu
cumod = mods.rx('Custandmod')[0]
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

#Standardize to 0-100
predcu19index = os.path.join(predgdb, 'predcu19index')
mincupred = float(arcpy.GetRasterProperties_management (predcu19, "MINIMUM").getOutput (0))
maxcupred = float(arcpy.GetRasterProperties_management (predcu19, "MAXIMUM").getOutput (0))
predcuindex = 100*(Float(Raster(predcu19))- mincupred)/(maxcupred - mincupred)
predcuindex.save(predcu19index)

#Export attribute table
predcu_land = arcpy.sa.Con((Raster(NLCDreclass) != 11) & (Raster(PSdissolveras) == 1), Raster(predcu19))
predcu_land.save(predcu19 + 'land')
arcpy.CopyRows_management(predcu19 + 'land', out_table = os.path.join(predgdb, 'predcu_tab'))

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

#Export attribute table
predpi_land = arcpy.sa.Con((Raster(NLCDreclass) != 11) & (Raster(PSdissolveras) == 1), Raster(predpi30))
predpi_land.save(predpi30 + 'land')
arcpy.CopyRows_management(predpi30 + 'land', out_table = os.path.join(predgdb, 'predpi_tab'))

#Produce boolean vector showing 1 - land pixels with Cu or Zn beyond background levels and 0 - pixels at background level
predallbool_land = arcpy.sa.Con((Raster(predcu19index)>0) | (Raster(predzn36index)>0), 1, 0)
predallbool_land.save(os.path.join(predgdb, 'predcuznbool_land'))

#Extract results to check that they match R outputs
XRFsites = os.path.join(pollutgdb, "XRFsites_aeasel")
sitescheck = os.path.join(predgdb,'sitescheck')
arcpy.CopyFeatures_management(XRFsites, sitescheck)
arcpy.sa.ExtractMultiValuesToPoints(sitescheck, [predzn36, predcu19, predpi30], bilinear_interpolate_values='NONE')

#Get residuals
arcpy.MakeFeatureLayer_management(XRFsites, 'sites')
arcpy.AddJoin_management('sites', 'SiteIDPair', pollutfieldclean_cast, 'SiteIDPair')
arcpy.CopyFeatures_management('sites', os.path.join(rootdir, 'results//PSpredictions.gdb/XRFsites_aeasel_resid'))