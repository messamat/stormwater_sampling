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
oddclasfiles = [path.join(res, f) for f in mlclist if (int(f[7:9]) % 2 > 0)]
evenclasfiles = [path.join(res, f) for f in mlclist  if (int(f[7:9]) % 2 == 0)]

# for tile in evenclasfiles:
#     print(tile)
#     for i in range(28,39):
#         arcpy.CopyRaster_management(tile, tile.replace('16',str(i)))

tileset = defaultdict(list)
for i in evenclasfiles+oddclasfiles:
    tileset[path.split(i)[1][13:22]].append(i)

ntiles = max([len(v) for v in tileset.values()])
def bingmean(tile, N=ntiles):
    bingmean = arcpy.sa.Float(arcpy.sa.CellStatistics(tileset[tile], statistics_type='SUM', ignore_nodata='DATA')) / N
    bingmean.save(path.join(res, 'mean{}'.format(tile)))
tilesetk = tileset.keys()


tic = time.time()
p = ProcessPool(nodes=4)
p.map(bingmean, tilesetk)
p.close
# for tile in tileset:
#     print(tile)
#     bingmean(tile, 10)
print(time.time()-tic)


#53 tiles, 22 layers - 122s - with multiprocessing: 56s (56/(53*22))=0.048027 #see parallel_test2.py
#1500 tiles, 168 layers - 1500*168*0.048027=12103/3600 = 3.5 h


#Compute mean across time for even and odd hours separately - 5h for a week of data
t0 = time.time()
arcpy.env.extent = "MAXOF"
if not arcpy.Exists(bingmean_odd):
    bingmean = CellStatistics(oddclasfiles, statistics_type='MEAN', ignore_nodata='DATA')
    bingmean.save(bingmean_odd) #Saving to .tif within /bing does not work, saving to grid within F: works
    print('Done averaging traffic for odd hours')
if not arcpy.Exists(bingmean_even):
    bingmean = CellStatistics(evenclasfiles, statistics_type='MEAN', ignore_nodata='DATA')
    bingmean.save(bingmean_even)
    print('Done averaging traffic for even hours')
t1=time.time()
print(t1-t0)

#Compute mean for all days at 3 am to identify semi-permanent road closures
if not arcpy.Exists(traffic3am_ras):
    arcpy.env.extent = "MAXOF"
    print('Produce 3 am traffic layer')
    traffic3am_raslist = arcpy.ListRasters('reclas_*_03_00')
    traffic3am_mean  = CellStatistics(traffic3am_raslist, statistics_type='MEAN', ignore_nodata='DATA')
    traffic3am_mean.save(traffic3am_ras)
if not arcpy.Exists(traffic2am_ras):
    arcpy.env.extent = "MAXOF"
    print('Produce 2 am traffic layer')
    traffic2am_raslist = arcpy.ListRasters('reclas_*_02_00')
    traffic2am_mean  = CellStatistics(traffic2am_raslist, statistics_type='MEAN', ignore_nodata='DATA')
    traffic2am_mean.save(traffic2am_ras)

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