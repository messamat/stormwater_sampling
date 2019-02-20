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