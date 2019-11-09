#This

import multiprocessing
from functools import partial
import arcpy
from arcpy.sa import *
import traceback
import os
import re
import itertools
from collections import defaultdict
import time
from heatmap_custom import *

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
#rootdir = 'F:/Mathis/Levin_Lab/stormwater'
#res = #os.path.join(rootdir, 'results/airdata/tiles')
rootdir = 'F:/Mathis/Levin_Lab/stormwater'
res = os.path.join(rootdir, 'results/bing')

'Useful resources:' \
'- https://pythongisandstuff.wordpress.com/2013/07/31/using-arcpy-with-multiprocessing-part-3/' \
'- https://gis.stackexchange.com/questions/140533/can-multiprocessing-with-arcpy-be-run-in-a-script-tool' \
'- https://gis.stackexchange.com/questions/17869/multiprocessing-errors-arcgis-implementation' \
'- https://www.esri.com/arcgis-blog/products/arcgis-desktop/analytics/multiprocessing-with-arcgis-approaches-and-considerations-part-1/?rmedium=blogs_esri_com&rsource=/esri/arcgis/2012/09/26/distributed-processing-with-arcgis-part-1/' \
'- https://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/parallel-processing-with-spatial-analyst.htm'

def bingmean(tile, tilediceven, Neven, tiledicodd, Nodd, tiledic2am, N2am, tiledic3am, N3am,
              kernel_dir, keyw, outdir):
    tic =time.time()
    outras = 'bing{}'.format(tile)
    print(outras)

    if not arcpy.Exists(outras):
        tmpdir = os.path.join(os.path.dirname(outdir),'tmp_{}'.format(str(tile)))
        try:
            os.mkdir(tmpdir)
            arcpy.env.scratchWorkspace = tmpdir #Not having a separate scratch workspace can lead to bad locking issues

            arcpy.env.extent = "MAXOF"
            #Compute mean bing index
            bingodd = arcpy.sa.Float(arcpy.sa.CellStatistics(tiledicodd[tile], statistics_type='SUM', ignore_nodata='DATA')) / Nodd
            bingeven = arcpy.sa.Float(arcpy.sa.CellStatistics(tilediceven[tile], statistics_type='SUM', ignore_nodata='DATA')) / Neven

            #Remove road closure pixels (classified as maximum congestion - red - for a long duration often including at night) and make sure not to average over logo pixels as would half values
            #change threshold for bingmean3am and 2am to 3 when have more than one day of data
            #Think of simplifying logo approach with http://geoinformaticstutorial.blogspot.com/2012/09/reading-raster-data-with-python-and-gdal.html setgeotransform and getgeotransform with one signle logopix tile
            if len(tiledic2am[tile]) > 0 & len(tiledic3am[tile]) > 0:
                bingmean2am = arcpy.sa.Float(arcpy.sa.CellStatistics(tiledic2am[tile], statistics_type='SUM', ignore_nodata='DATA')) / N2am
                bingmean3am = arcpy.sa.Float(arcpy.sa.CellStatistics(tiledic3am[tile], statistics_type='SUM', ignore_nodata='DATA')) / N3am
                bingclean = Int(1000*
                                Con(IsNull(bingeven), #If even hours are null
                                    Con(IsNull(bingodd), bingodd,
                                        Con((bingodd > 0.1),
                                            Con((bingodd < 2),
                                                Con(IsNull(bingmean3am), bingodd,
                                                    Con((bingmean3am< 3), bingodd
                                                        )
                                                    )
                                                ),
                                            Con((bingodd < 2),
                                                Con(IsNull(bingmean3am), (bingodd*Nodd)/(Nodd+Neven),
                                                    Con((bingmean3am< 3), (bingodd*Nodd)/(Nodd+Neven)
                                                        )
                                                    )
                                                )
                                            )
                                        ),
                                    Con(IsNull(bingodd), #If odd hours are null
                                        Con((bingeven > 0.1),
                                            Con((bingeven < 2),
                                                Con(IsNull(bingmean2am), bingeven,
                                                    Con((bingmean2am< 3), bingeven
                                                        )
                                                    )
                                                ),
                                            Con((bingeven < 2),
                                                Con(IsNull(bingmean2am), (bingeven*Neven)/(Nodd+Neven),
                                                    Con((bingmean2am< 3), (bingeven*Neven)/(Nodd+Neven)
                                                        )
                                                    )
                                                )
                                            ),
                                        Con((bingeven < 2) & (bingodd < 2), #If neither even nor odd hours are null
                                            Con(IsNull(bingmean2am),
                                                Con(IsNull(bingmean3am), bingeven*Neven+bingodd*Nodd/(Nodd+Neven),
                                                    Con(bingmean3am < 3, bingeven*Neven+bingodd*Nodd/(Nodd+Neven)
                                                        )
                                                    ),
                                                Con(bingmean2am < 3,
                                                    Con(IsNull(bingmean3am), bingeven*Neven+bingodd*Nodd/(Nodd+Neven),
                                                        Con(bingmean3am < 3, bingeven*Neven+bingodd*Nodd/(Nodd+Neven)
                                                            )
                                                        )
                                                    )
                                                )
                                            )
                                        )
                                    )
                                )
            else:
                bingclean = Int(1000*
                                Con(IsNull(bingeven), #If even hours are null
                                    Con(IsNull(bingodd), bingodd,
                                        Con((bingodd > 0.1),
                                            Con((bingodd < 2), bingodd),
                                            Con((bingodd < 2), (bingodd*Nodd)/(Nodd+Neven))
                                            )
                                        ),
                                    Con(IsNull(bingodd), #If odd hours are null
                                        Con((bingeven > 0.1),
                                            Con((bingeven < 2), bingeven),
                                            Con((bingeven < 2),(bingeven*Neven)/(Nodd+Neven))
                                            ),
                                        Con((bingeven < 2) & (bingodd < 2), #If neither even nor odd hours are null
                                            bingeven*Neven+bingodd*Nodd/(Nodd+Neven)
                                            )
                                        )
                                    )
                                )
            #Compute heatmap
            customheatmap(kernel_dir=kernel_dir, in_raster=bingclean, scratch_dir = tmpdir,
                      out_gdb=outdir, out_var= outras, divnum=100, keyw=keyw, ext='.tif',
                      verbose = False)
            print('Done generating heatmap')
        except Exception:
            traceback.print_exc()
            arcpy.Delete_management(tmpdir)
        if arcpy.Exists(tmpdir):
            arcpy.Delete_management(tmpdir)
    else:
        print('{} already exists...'.format(outras))
    #bingclean.save(outras) #DO NOT OUTPUT TO GRID IN PARALLEL
    print(time.time()-tic)


if __name__ == '__main__':
    arcpy.env.workspace = res
    # Create list of layers
    regextile = re.compile("^[0-9]{6}_[0-9]{2}_[0-9]{2}.*mlc[.]tif$")
    mlclist = list(itertools.chain.from_iterable( #To unnest list
        [filter(regextile.search, filenames) for  (dirpath, dirnames, filenames) in os.walk(res)]))

    #Create dictionary associating each tile with all of the hourly data available for it and get a sorted list of all date times
    tileeven = defaultdict(list)
    tileodd = defaultdict(list)
    datetimes = []
    ntileeven = 0
    ntileodd = 0
    for i in mlclist:
        if (int(i[7:9]) % 2) > 0:
            tileodd[re.sub('[a-zA-Z]', '', os.path.split(i)[1][13:22])].append(os.path.join(res, i))
        elif (int(i[7:9]) % 2) == 0:
            tileeven[re.sub('[a-zA-Z]', '', os.path.split(i)[1][13:22])].append(os.path.join(res, i))
        datetimes.append(os.path.split(i)[1][0:12])
    daterange = list(sorted(set(datetimes)))
    ntileeven = len([h for h in daterange if (int(h[7:9]) % 2) == 0])
    ntileodd = len([h for h in daterange if (int(h[7:9]) % 2) > 0])

    #Get list of tiles at 2 am and 3 am to detect road closures
    traffic2am_ras = os.path.join(res, 'traffic2am')
    traffic3am_ras = os.path.join(res, 'traffic3am')
    tileset2am = {}
    tileset3am = {}
    re2am = re.compile(".*[0-9]{6}_02_[0-9]{2}.*[.]tif$")
    re3am = re.compile(".*[0-9]{6}_03_[0-9]{2}.*[.]tif$")
    for key in tileeven:
        tileset2am[key] = filter(re2am.search, tileeven[key])
        tileset3am[key] = filter(re3am.search, tileodd[key])

    ntiles2am = len([h for h in daterange if h[7:9]=='02'])
    ntiles3am = len([h for h in daterange if h[7:9]=='03'])

    #Get all keys
    tilelist = list(sorted(set(tileeven.keys()+tileodd.keys())))

    print('Launch parallel processing')
    tic = time.time()
    p = multiprocessing.Pool(int(multiprocessing.cpu_count()/2))
    bingmean_partial = partial(bingmean,
                               tilediceven=tileeven, Neven=ntileeven,
                               tiledicodd=tileodd, Nodd=ntileodd,
                               tiledic2am=tileset2am, N2am=ntiles2am,
                               tiledic3am=tileset3am, N3am=ntiles3am,
                               kernel_dir=os.path.join(rootdir, 'results/bing'),
                               keyw='log300',
                               outdir=res) #Set N parameter as constant
    p.map(bingmean_partial, tilelist)
    p.close()
    print(time.time() - tic)