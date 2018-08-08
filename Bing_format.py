import arcpy
import glob
import re
import time
from os import *
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

#Set up paths
rootdir = 'C:/Mathis/ICSL/stormwater/'
res= path.join(rootdir, 'results/bing/')
arcpy.env.workspace=res
mlclist = [filenames for (dirpath, dirnames, filenames) in walk(res)][0]
regex=re.compile(".+mlc\.tif$")
clasfiles = [m.group(0) for l in mlclist for m in [regex.search(l)] if m] #Look into this
kernel = NbrWeight(path.join(rootdir, 'results/logkernel100.txt')) #UPDATE

#Note that 180623_02_45 was deleted because too incomplete and 180625_22_45 is missing tiles (scheduled task crashed)

#Reclassify (2.5 min/hourly image - should take about 1h/day of data)
remap = RemapValue([[1, 0], [2, 0], [3, 0], [4, 0], [5, 2], [6, 3], [7, 0], [8, 1],['NODATA',0]])
for f in clasfiles:
    t0 = time.time()
    print(f)
    outf = '{}reclass.tif'.format(f[:-13])
    if arcpy.Exists(outf):
        print('{} already exists'.format(outf))
    else:
        outReclass = Reclassify(f, 'Value', remap, "DATA")
        outReclass.save(outf)
    print(time.time() - t0)

#Compute statistics across time - 5h for a week of data
mlclist = [filenames for (dirpath, dirnames, filenames) in walk(res)][0]
regex=re.compile(".+reclass\.tif$")
reclasfiles = [m.group(0) for l in mlclist for m in [regex.search(l)] if m]
t0 = time.time()
arcpy.env.extent = "MAXOF"
bingmeanout = 'bingmean1806.tif'
bingmean = CellStatistics(reclasfiles, statistics_type='MEAN', ignore_nodata='DATA')
bingmean.save(bingmeanout)
t1=time.time()
print(t1-t0)


#Majority across time
#Remove those that have a majority of 3

#Instead, use
bingcleanout = 'bingmeanclean1806.tif'
bingclean = Con(Raster(bingmeanout)<2.5, Raster(bingmeanout), 0)
bingclean.save(bingcleanout)

#Generate heat map
heat_bing = FocalStatistics(path.join(res,bingcleanout), neighborhood=kernel, statistics_type='SUM', ignore_nodata='DATA') #NEed to be in geodatabse
heat_bing.save('heat_bing')
heat_bing_int = Int(Raster('heat_bing')+0.5)
heat_bing_int.save('heat_bing_int')


#Integerize bing heat map
arcpy.ProjectRaster_management(heat_bing, 'heat_bing_proj', UTM10, resampling_type='BILINEAR') #Project
heat_bing_int = Int(Raster('heat_bing_proj')+0.5)
heat_bing_int.save('heat_bing_int')

#Compute a heat index out of 100 (standardized, but non-transformed -- for communication)
bingmax = arcpy.GetRasterProperties_management('heat_bing_int', 'MAXIMUM')
bingheatindex = 100*Float(Raster('heat_bing_int'))/float(bingmax.getOutput(0))
bingheatindex.save('heat_bing_index2')