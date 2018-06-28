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

#Note that 180623_02_45 was deleted because too incomplete and 180625_22_45 is missing tiles (scheduled task crashed)

#Reclassify (2.5 min/hourly image - should take about 1h/day of data)
remap = RemapValue([[1, 'NODATA'], [2, 0], [3, 'NODATA'], [4, 'NODATA'], [5, 2], [6, 3], [7, 'NODATA'], [8, 1]])
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

#Compute statistics across time (30 min for 19 rasters, expect ~4-5h for )
mlclist = [filenames for (dirpath, dirnames, filenames) in walk(res)][0]
regex=re.compile(".+reclass\.tif$")
reclasfiles = [m.group(0) for l in mlclist for m in [regex.search(l)] if m] #Look into this
t0 = time.time()
bingmean = CellStatistics(reclasfiles, statistics_type='MEAN', ignore_nodata='DATA')
bingmean.save('bingmean1806.tif')
print(time.time()-t0)



#bingmean = CellStatistics(reclasfiles, statistics_type='MEAN', ignore_nodata='DATA')
