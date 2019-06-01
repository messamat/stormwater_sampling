""""
Comment: in February 2019, changed Bing congestion post-processing to reclassify green:0, yellow:1, orange:2, red:3 to
green:1, yellow:2, orangeL3, red:4. This led to regenerating intermediate outputs with 2 as a suffix and bingeuc1902
or e.g. heat_bing1902log500"""

import arcpy
import math
import glob
import re
import time
from os import *
from arcpy.sa import *
from heatmap_custom import *
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater/'
res= path.join(rootdir, 'results/bing/')
arcpy.env.workspace=res
mlclist = [filenames for (dirpath, dirnames, filenames) in walk(res)][0]
regex=re.compile(".+mlc\.tif$")
clasfiles = [m.group(0) for l in mlclist for m in [regex.search(l)] if m] #Look into this
UTM10 = arcpy.SpatialReference(26910)

sites = path.join(rootdir, "data\\field_data\sampling_sites_edit.shp")

#Bing logo masks
logomask_odd = path.join(rootdir, 'results/boolean_logo.tif')
logomask_even = path.join(rootdir, 'results/boolean_logoalt.tif')

#Note that 180623_02_45 was deleted because too incomplete and 180625_22_45 is missing tiles (scheduled task crashed)
gdb = path.join(rootdir,'results/bing/postprocess.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir,'results/bing'), 'postprocess.gdb')

#Output variables (2 because changing from G0Y1O2R3 to G1Y2O3R4)
bingmean_odd = path.join(res,'bingmeanod2')
bingmean_even = path.join(res,'bingmeanev2')
traffic3am_ras = path.join(res, 'traffic3am2')
traffic2am_ras = path.join(res, 'traffic2am2')

bingclean_odd = path.join(res, 'bingcleanod2')
bingclean_even = path.join(res, 'bingcleanev2')
bingclean_mean = path.join(gdb, 'bingcleanmea2')
bingeuc = path.join(res, 'bingeucdir/bingeuc2')

# arcpy.Project_management(sites, path.join(rootdir, "data\\field_data\sampling_sites_edit_proj.shp"), out_coor_system=path.join(res, clasfiles[0]))
# sitesbuf =  path.join(rootdir, "data\\field_data\sampling_sites_edit1000mbuffer.shp")
# arcpy.Buffer_analysis(path.join(rootdir, "data\\field_data\sampling_sites_edit_proj.shp"), sitesbuf
#                      , '1000 meters', method='GEODESIC')

#Reclassify (2 min/hourly image when out of a gdb, then 8 minutes when using gdb?)
arcpy.env.workspace= gdb
# if not arcpy.ListRasters('reclas2_*'):
remap = RemapValue([[1, 0], [2, 1], [3, 0], [4, 0], [5, 3], [6, 4], [7, 0], [8, 2],['NODATA',0]])
for f in clasfiles:
    t0 = time.time()
    print(f)
    outf = path.join(gdb, 'reclas2_{}'.format(f[:-14]))
    # if arcpy.Exists(outf):
    #     print('{} already exists'.format(outf))
    # else:
    outReclass = Reclassify(path.join(res, f), 'Value', remap, "DATA")
    outReclass.save(outf)
    print(time.time() - t0)

#Create list of layers for odd and even hours to remove bing logo
reclasfiles = arcpy.ListRasters('reclas2_*')
oddclasfiles = [path.join(gdb, f) for f in reclasfiles if (int(f[15:17]) % 2 > 0)]
evenclasfiles = [path.join(gdb, f) for f in reclasfiles if (int(f[15:17]) % 2 == 0)]

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
xpl = range(int(totalextent.XMin), int(math.ceil(totalextent.XMax)),
            int(math.ceil((totalextent.XMax - totalextent.XMin)/n_2d))+1) + [int(math.ceil(totalextent.XMax))]
ypl = range(int(totalextent.YMin), int(math.ceil(totalextent.YMax)),
            int(math.ceil((totalextent.YMax - totalextent.YMin)/n_2d))+1) + [int(math.ceil(totalextent.YMax))]
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
                             Con((Raster(bingmean_odd)<2.0) & IsNull(Raster(logomask_odd)) & (Raster(traffic3am_ras)<0.4) &
                             (Raster(bingmean_even)<2.0) & IsNull(Raster(logomask_even)) & (Raster(traffic2am_ras)<0.4),
                             (Raster(bingmean_odd)+Raster(bingmean_even))/2,
                             Con((Raster(bingmean_odd)<2.0) & IsNull(Raster(logomask_odd)) & (Raster(traffic3am_ras)<0.4),
                                 Raster(bingmean_odd),
                                 Con((Raster(bingmean_even)<2.0) & IsNull(Raster(logomask_even)) & (Raster(traffic2am_ras)<0.4),
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
#arcpy.env.extent = sitesbuf
if not arcpy.Exists(path.join(gdb, 'bingeuc1902')):
    arcpy.env.workspace = os.path.join(res, 'bingeucdir')
    arcpy.MosaicToNewRaster_management(arcpy.ListRasters('bingeuc*'), gdb, 'bingeuc1902', coordinate_system_for_the_raster= UTM10,
                                       pixel_type = '16_BIT_UNSIGNED', number_of_bands=1, mosaic_method='LAST')

#Create heatmap - focal statistics
customheatmap(kernel_dir=path.join(rootdir, 'results/bing'), in_raster=path.join(gdb, 'bingeuc1902'),
              out_gdb = gdb, out_var='bing1902', divnum=100, keyw='log[3]00')


#Project, and compute heat index out of 100
arcpy.ProjectRaster_management(path.join(gdb, 'bingeuc2'), 'bingeuc2_proj', UTM10, resampling_type='BILINEAR') #Project
#Compute a heat index out of 100 (standardized, but non-transformed -- for communication)
bingmax = arcpy.GetRasterProperties_management('heat_bing_int', 'MAXIMUM')
bingheatindex = 100*Float(Raster('heat_bing_int'))/float(bingmax.getOutput(0))
bingheatindex.save('heat_bing_index')