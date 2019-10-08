__author__ = 'Mathis Messager'
import multiprocessing
from functools import partial
import arcpy
import os
from collections import defaultdict
import time
from heatmap_custom import *
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

def roadtoheat(site, inshp, res, kernel_dir, keyw, inFID, heatfield, sitedic, outdir):
    expr =  """{0} = {1}""".format(arcpy.AddFieldDelimiters(inshp, inFID), str(site))
    print(expr)
    arcpy.MakeFeatureLayer_management(in_features=inshp, out_layer='lyr')
    arcpy.SelectLayerByAttribute_management('lyr', selection_type='NEW_SELECTION', where_clause=expr)
    #print(site)
    nshp = len([row[0] for row in arcpy.da.SearchCursor('lyr', [inFID])])  #int(arcpy.GetCount_management('lyr').getOutput(0)) would be faster but often has a bug
    #print(nshp)
    if nshp > 0:
        #print('{} features'.format(nshp))
        arcpy.ResetEnvironments()
        arcpy.env.extent = sitedic[site]
        outras =  os.path.join(outdir, 'hpmstiger_{0}{1}.tif'.format(heatfield, site))
        if not arcpy.Exists(outras):
            print('{} does not exist, generate heatmap'.format(outras))
            tmpdir = os.path.join(os.path.dirname(outdir),'tmp_{}'.format(str(site)))
            os.mkdir(tmpdir)
            arcpy.env.scratchWorkspace = tmpdir
            arcpy.PolylineToRaster_conversion('lyr', value_field=heatfield, out_rasterdataset= outras,
                                              priority_field=heatfield,cellsize=res)
            customheatmap(kernel_dir=kernel_dir, in_raster=outras, scratch_dir = tmpdir,
                          out_gdb=outdir, out_var= heatfield+str(site), divnum=100, keyw=keyw, ext='.tif',
                          verbose = True)
            arcpy.Delete_management(tmpdir)
    arcpy.SelectLayerByAttribute_management('lyr', selection_type='CLEAR_SELECTION') #might not be necessary
    arcpy.Delete_management('lyr') #might not be necessary

if __name__ == '__main__':
    rootdir = 'D:/Mathis/ICSL//stormwater'
    template_ras = os.path.join(rootdir,'results/bing/181204_02_00_class_mlc.tif')
    res = arcpy.GetRasterProperties_management(template_ras, 'CELLSIZEX').getOutput(0)

    AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')
    AQIdir = os.path.join(rootdir, 'results/airdata')
    AQIsites_bufunion = os.path.join(rootdir, 'results/airdata/airsites_550bufunion.shp')
    roadAQI = os.path.join(AQIgdb, 'AQI_hpmstigerinters')

    #[f.name for f in arcpy.ListFields(OSMAQIproj)]
    AQIbuflist = defaultdict(list)
    with arcpy.da.SearchCursor(AQIsites_bufunion,['FID','SHAPE@']) as cursor:
        for row in cursor:
            ex = row[1].extent
            AQIbuflist[row[0]] = " ".join([str(c) for c in [ex.XMin, ex.YMin, ex.XMax, ex.YMax]])
    keylist = AQIbuflist.keys()
    print('{} sites to process'.format(str(len(keylist))))

    #For troubleshooting:
    # OSMtoheat_partial = partial(OSMtoheat,
    #                             inshp = OSMAQIproj,
    #                             res=res,
    #                             kernel_dir = os.path.join(rootdir, 'results/bing'),
    #                             keyw = 'log300',
    #                             inFID = 'FID_airsites_600bufunion',
    #                             heatfield = 'fclassSPD',
    #                             sitedic = AQIbuflist,
    #                             outdir = AQIdir)
    #
    # for key in keylist[1:20]:
    #     OSMtoheat_partial(key)


    tic = time.time()
    p = multiprocessing.Pool(int(multiprocessing.cpu_count()/2))
    roadtoheat_partial = partial(roadtoheat,
                                inshp = roadAQI,
                                res=res,
                                kernel_dir = os.path.join(rootdir, 'results/bing'),
                                keyw = 'log100',
                                inFID = 'FID_airsites_550bufunion',
                                heatfield = 'aadt_filled',
                                sitedic = AQIbuflist,
                                outdir = AQIdir)
    print('Launch parallel processing')
    p.map(roadtoheat_partial, keylist)
    p.close()
    print(time.time() - tic)