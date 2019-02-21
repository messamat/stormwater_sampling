import arcpy
import math
import glob
import re
import time
from os import *
import arcpy
from pathos.multiprocessing import ProcessPool
from collections import defaultdict
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

#Set up paths
#rootdir = 'F:/Levin_Lab/stormwater/'
rootdir = 'C:/Mathis/ICSL/stormwater'
#res= path.join(rootdir, 'results/bing/')
res = path.join(rootdir, 'results/airdata/tiles')
arcpy.env.workspace=res

#UTM10 = arcpy.SpatialReference(26910)

#Note that 180623_02_45 was deleted because too incomplete and 180625_22_45 is missing tiles (scheduled task crashed)
gdb = path.join(rootdir,'results/airdata/postprocess.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir,'results/bing'), 'postprocess.gdb')

#Output variables
bingmean_odd = path.join(res,'bingmeanod')
bingmean_even = path.join(res,'bingmeanev')
traffic3am_ras = path.join(res, 'traffic3am')
traffic2am_ras = path.join(res, 'traffic2am')

bingclean_odd = path.join(res, 'bingcleanod')
bingclean_even = path.join(res, 'bingcleanev')
bingclean_mean = path.join(gdb, 'bingcleanmea')
bingeuc = path.join(res, 'bingeuc')


#Create list of layers for odd and even hours to remove bing logo
mlclist = arcpy.ListRasters('*tif')


#Clean out bing logo and road closure artefacts then run euclidean allocation on tiles
arcpy.ResetEnvironments()
ext_template = bingmean_odd
totalextent = arcpy.Describe(ext_template).Extent
n_2d = 20 #Number of side tiles (total number of tiles = n_2d^2)
xpl = range(int(totalextent.XMin), int(math.ceil(totalextent.XMax)), int(math.ceil((totalextent.XMax - totalextent.XMin)/n_2d))+1) + [int(math.ceil(totalextent.XMax))]
ypl = range(int(totalextent.YMin), int(math.ceil(totalextent.YMax)), int(math.ceil((totalextent.YMax - totalextent.YMin)/n_2d))+1) + [int(math.ceil(totalextent.YMax))]
x=0
for xmin, xmax in zip(xpl[:-1], xpl[1:]):
    for ymin, ymax in zip(ypl[:-1], ypl[1:]):
        tilext = arcpy.Extent(xmin, ymin, xmax, ymax)
        print(tilext)
        arcpy.env.extent = tilext
        arcpy.env.snapRaster = ext_template
        #Odd hour tile processing
        outeuc = bingeuc + str(x)

        try:
            arcpy.Delete_management(clean_mean)
        except:
            pass
        try:
            arcpy.Delete_management(euc)
        except:
            pass

        if not arcpy.Exists(outeuc):
            print(outeuc)
            clean_mean = Int(1000*
                             Con((Raster(bingmean_odd)<1.0) & IsNull(Raster(logomask_odd)) & (Raster(traffic3am_ras)<0.2) &
                             (Raster(bingmean_even)<1.0) & IsNull(Raster(logomask_even)) & (Raster(traffic2am_ras)<0.2),
                             (Raster(bingmean_odd)+Raster(bingmean_even))/2,
                             Con((Raster(bingmean_odd)<1.0) & IsNull(Raster(logomask_odd)) & (Raster(traffic3am_ras)<0.2),
                                 Raster(bingmean_odd),
                                 Con((Raster(bingmean_even)<1.0) & IsNull(Raster(logomask_even)) & (Raster(traffic2am_ras)<0.2),
                                     Raster(bingmean_even)
                                     )
                                 )
                             ))
            print('Done cleaning and averaging, running euclidean allocation')
            euc = arcpy.sa.EucAllocation(clean_mean, maximum_distance=50)
            euc.save(bingeuc)
            arcpy.Rename_management(bingeuc, outeuc)
        x+=1

#Mosaic all tiles
arcpy.ResetEnvironments()
if not arcpy.Exists(path.join(gdb, 'bingeuc')):
    arcpy.env.workspace = res
    arcpy.MosaicToNewRaster_management(arcpy.ListRasters('bingeuc*'), gdb, 'bingeuc', coordinate_system_for_the_raster= UTM10,
                                       pixel_type = '16_BIT_UNSIGNED', number_of_bands=1, mosaic_method='LAST')

#Create heatmap - focal statistics
keyw = '_pow.*(?<=2)[.]txt'
#keyw = '_log'
for kertxt in listdir(path.join(rootdir, 'results/bing')):
    if re.compile('kernel'+keyw).match(kertxt):
        outext = '_bingheat_bing{}'.format(path.splitext(kertxt)[0][7:])
        if not arcpy.Exists(path.join(gdb,outext)):
            print(outext)
            tic = time.time()
            kernel = NbrWeight(path.join(rootdir, 'results/bing', kertxt))
            heat_bing = Int(arcpy.sa.FocalStatistics(path.join(gdb, 'bingeuc'), neighborhood= kernel, statistics_type='SUM', ignore_nodata='DATA')/100+0.5)
            heat_bing.save(path.join(gdb,outext))
            toc = time.time()
            print(toc-tic)

#Project, and compute heat index out of 100
arcpy.ProjectRaster_management(path.join(gdb, 'bingeuc'), 'bingeuc_proj', UTM10, resampling_type='BILINEAR') #Project
#Compute a heat index out of 100 (standardized, but non-transformed -- for communication)
bingmax = arcpy.GetRasterProperties_management('heat_bing_int', 'MAXIMUM')
bingheatindex = 100*Float(Raster('heat_bing_int'))/float(bingmax.getOutput(0))
bingheatindex.save('heat_bing_index')