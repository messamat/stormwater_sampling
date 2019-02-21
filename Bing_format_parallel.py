import multiprocessing
from functools import partial
import arcpy
import os
from collections import defaultdict
import time
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
rootdir = 'C:/Mathis/ICSL/stormwater'
res = os.path.join(rootdir, 'results/airdata/tiles')
'Useful resources:' \
'- https://pythongisandstuff.wordpress.com/2013/07/31/using-arcpy-with-multiprocessing-part-3/' \
'- https://gis.stackexchange.com/questions/140533/can-multiprocessing-with-arcpy-be-run-in-a-script-tool' \
'- https://gis.stackexchange.com/questions/17869/multiprocessing-errors-arcgis-implementation' \
'- https://www.esri.com/arcgis-blog/products/arcgis-desktop/analytics/multiprocessing-with-arcgis-approaches-and-considerations-part-1/?rmedium=blogs_esri_com&rsource=/esri/arcgis/2012/09/26/distributed-processing-with-arcgis-part-1/' \
'- https://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/parallel-processing-with-spatial-analyst.htm'
def bingmean(tile, N):
    outname = os.path.split(tile[0])[1][13:22]
    bingmean = arcpy.sa.Float(arcpy.sa.CellStatistics(tile, statistics_type='SUM', ignore_nodata='DATA')) / N
    bingmean.save(os.path.join(res, 'mean{}.tif'.format(outname))) #DO NOT OUTPUT TO GRID IN PARALLEL

if __name__ == '__main__':
    arcpy.env.workspace = res
    # Create list of layers for odd and even hours to remove bing logo
    mlclist = arcpy.ListRasters('*tif')
    print(res)
    oddclasfiles = [os.path.join(res, f) for f in mlclist if (int(f[7:9]) % 2 > 0)]
    evenclasfiles = [os.path.join(res, f) for f in mlclist if (int(f[7:9]) % 2 == 0)]

    tileset = defaultdict(list)
    for i in evenclasfiles + oddclasfiles:
        tileset[os.path.split(i)[1][13:22]].append(i)

    ntiles = max([len(v) for v in tileset.values()])

    print('Launch parallel processing')
    tic = time.time()
    p = multiprocessing.Pool(4)
    bingmean_partial = partial(bingmean, N=ntiles) #Set N parameter as constant
    p.map(bingmean_partial, tileset.values())
    p.close()
    print(time.time() - tic)

    # Compute mean for all days at 3 am to identify semi-permanent road closures
    traffic3am_ras = os.path.join(res, 'traffic3am')
    traffic2am_ras = os.path.join(res, 'traffic2am')
    if not arcpy.Exists(traffic3am_ras):
        arcpy.env.extent = "MAXOF"
        print('Produce 3 am traffic layer')
        traffic3am_raslist = arcpy.ListRasters('*_26_*.tif')
        ############### TO KEEP CODING #########################
        tileset_3am = defaultdict(list)
        for i in evenclasfiles + oddclasfiles:
            tileset_3am[os.path.split(i)[1][13:22]].append(i)

        ntiles = max([len(v) for v in tileset.values()])

        traffic3am_mean = arcpy.sa.CellStatistics(traffic3am_raslist, statistics_type='MEAN', ignore_nodata='DATA')
        traffic3am_mean.save(traffic3am_ras)

    if not arcpy.Exists(traffic2am_ras):
        arcpy.env.extent = "MAXOF"
        print('Produce 2 am traffic layer')
        traffic2am_raslist = arcpy.ListRasters('reclas_*_02_00')
        traffic2am_mean = CellStatistics(traffic2am_raslist, statistics_type='MEAN', ignore_nodata='DATA')
        traffic2am_mean.save(traffic2am_ras)