import arcpy
import math
import glob
import re
import time
from os import *
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

#Set up paths
rootdir = 'F:/Levin_Lab/stormwater/'
res= path.join(rootdir, 'results/bing/')
arcpy.env.workspace=res
mlclist = [filenames for (dirpath, dirnames, filenames) in walk(res)][0]
regex=re.compile(".+mlc\.tif$")
clasfiles = [m.group(0) for l in mlclist for m in [regex.search(l)] if m] #Look into this

kernel_log100 = NbrWeight(path.join(rootdir, 'results/logkernel100.txt')) #UPDATE

#Bing logo masks
logomask_odd = path.join(rootdir, 'results/boolean_logo.tif')
logomask_even = path.join(rootdir, 'results/boolean_logoalt.tif')

#Note that 180623_02_45 was deleted because too incomplete and 180625_22_45 is missing tiles (scheduled task crashed)
gdb = path.join(rootdir,'results/bing/postprocess.gdb')
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

#Reclassify (2 min/hourly image when out of a gdb, then 8 minutes when using gdb?)
arcpy.env.workspace= gdb
if not arcpy.ListRasters('reclas_*'):
    remap = RemapValue([[1, 0], [2, 0], [3, 0], [4, 0], [5, 2], [6, 3], [7, 0], [8, 1],['NODATA',0]])
    for f in clasfiles:
        t0 = time.time()
        print(f)
        outf = path.join(gdb, 'reclas_{}'.format(f[:-14]))
        if arcpy.Exists(outf):
            print('{} already exists'.format(outf))
        else:
            outReclass = Reclassify(path.join(res, f), 'Value', remap, "DATA")
            outReclass.save(outf)
        print(time.time() - t0)

#Create list of layers for odd and even hours to remove bing logo
reclasfiles = arcpy.ListRasters('reclas_*')
oddclasfiles = [path.join(gdb, f) for f in reclasfiles if (int(f[14:16]) % 2 > 0)]
evenclasfiles = [path.join(gdb, f) for f in reclasfiles if (int(f[14:16]) % 2 == 0)]

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
    traffic3am_raslist = arcpy.ListRasters('reclas_*_03_00')
    traffic3am_mean  = CellStatistics(traffic3am_raslist, statistics_type='MEAN', ignore_nodata='DATA')
    traffic3am_mean.save(traffic3am_ras)
if not arcpy.Exists(traffic2am_ras):
    traffic2am_raslist = arcpy.ListRasters('reclas_*_02_00')
    traffic2am_mean  = CellStatistics(traffic2am_raslist, statistics_type='MEAN', ignore_nodata='DATA')
    traffic2am_mean.save(traffic2am_ras)

#Clean out bing logo and road closure artefacts
#In arcmap: Int(Con(("bingmeanod"<1) & IsNull("boolean_logo.tif") & ("traffic3am"<0.2), 1000*"bingmeanod"))
arcpy.ResetEnvironments()
totalextent = arcpy.Describe(bingmean_odd).Extent
n_2d = 5 #Number of side tiles (total number of tiles = n_2d^2)
xpl = range(int(totalextent.XMin), int(math.ceil(totalextent.XMax)), int(math.ceil((totalextent.XMax - totalextent.XMin)/n_2d))) + [int(math.ceil(totalextent.XMax))]
ypl = range(int(totalextent.YMin), int(math.ceil(totalextent.YMax)), int(math.ceil((totalextent.YMax - totalextent.YMin)/n_2d))) + [int(math.ceil(totalextent.YMax))]
x=0
for xmin, xmax in zip(xpl[:-1], xpl[1:]):
    for ymin, ymax in zip(ypl[:-1], ypl[1:]):
        tilext = arcpy.Extent(xmin, ymin, xmax, ymax)
        print(tilext)
        arcpy.env.extent = tilext
        arcpy.env.snapRaster = bingmean_odd

        #Odd hour tile processing
        outodd = bingclean_odd + str(x)
        if not arcpy.Exists(outodd):
            print(outodd)
            bingclean = Int(Con((Raster(bingmean_odd)<1.0) & IsNull(Raster(logomask_odd)) & (Raster(traffic3am_ras)<0.2),
                        1000*Raster(bingmean_odd)))
            bingclean.save(outodd)
            arcpy.Delete_management(bingclean)

        #Even hour tile processing
        outeven = bingclean_even + str(x)
        if not arcpy.Exists(outeven):
            print(outeven)
            bingclean = Int(Con((Raster(bingmean_even)<1.0) & IsNull(Raster(logomask_even)) & (Raster(traffic2am_ras)<0.2),
                        1000*Raster(bingmean_even)))
            bingclean.save(outeven)
            arcpy.Delete_management(bingclean)
        x+=1

#Average over odd and even hours (mosaic to new raster with mean)
arcpy.ResetEnvironments()
arcpy.env.workspace = res
arcpy.MosaicToNewRaster_management(arcpy.ListRasters('bingclean*'), gdb, 'bingcleanmea', number_of_bands= 1,
                                   pixel_type = '16_BIT_UNSIGNED', mosaic_method='MEAN')
#Euclidean allocation
#Create raster of NoData value to set extent
bingeuc = path.join(gdb, 'bingeuc')
euc = arcpy.sa.EucAllocation(bingclean_mean, maximum_distance=50)
euc.save(bingeuc)

#Generate heat map for multiple kernel types (create for loops)
heat_bing = FocalStatistics(path.join(gdb, bingclean_mean), neighborhood=kernel_log100, statistics_type='SUM', ignore_nodata='DATA') #NEed to be in geodatabse
heat_bing.save('heat_bing')
heat_bing_int = Int(Raster('heat_bing')+0.5)
heat_bing_int.save('heat_bing_int')

#Integerize bing heat map, project, and compute heat index out of 100
arcpy.ProjectRaster_management(heat_bing, 'heat_bing_proj', UTM10, resampling_type='BILINEAR') #Project
heat_bing_int = Int(Raster('heat_bing_proj')+0.5)
heat_bing_int.save('heat_bing_int')

#Compute a heat index out of 100 (standardized, but non-transformed -- for communication)
bingmax = arcpy.GetRasterProperties_management('heat_bing_int', 'MAXIMUM')
bingheatindex = 100*Float(Raster('heat_bing_int'))/float(bingmax.getOutput(0))
bingheatindex.save('heat_bing_index')