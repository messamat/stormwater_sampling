#Author: Mathis Messager
#Date created: 2018/09/20

#Download effective, preliminary, and pending FEMA flood hazard layers
#See https://hazards.fema.gov/femaportal/wps/portal/NFHLWMS for available FEMA GIS data
#See https://support.esri.com/en/technical-article/000013950 for info on spatial reference

import arcpy
from os import *
import itertools
import glob
arcpy.env.overwriteOutput = True

#Folder structure
rootdir = 'C:/Mathis/ICSL/flood_vulnerability'
datadir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

FEMAdraftdir = os.path.join(datadir, 'FEMAdraft_20180920')
FEMAdraftgdb = os.path.join(FEMAdraftdir, 'FEMAdraft.gdb')
if not (os.path.exists(FEMAdraftdir)):
    print('Create new directory: {}'.format(FEMAdraftdir))
    os.mkdir(FEMAdraftdir)
if not arcpy.Exists(FEMAdraftgdb):
    print('Create FEMA geodatabase')
    arcpy.CreateFileGDB_management(FEMAdraftdir, 'FEMAdraft.gdb')

counties = os.path.join(datadir,"WA_County_Boundaries\WA_County_Boundaries.shp")

arcpy.env.workspace = FEMAdraftgdb



#Input geometry


def APIdownload(URL, name, itersize, ingeom, proj):
    if ingeom is str and


    IDrange = range(0,10000000, itersize)
    for i,j in itertools.izip(IDrange, IDrange[1:]):
        where = "OBJECTID>{0} AND OBJECTID<{1}".format(i,j)
        #&geometryType=esriGeometryPoint
        query = "?where={}&&text=&objectIds=&time=&geometry=&geometryType=esriGeometryPolygon&inSR=&\
        spatialRel=esriSpatialRelIntersects&relationParam=&outFields=*&returnGeometry=true&returnTrueCurves=false&\
        maxAllowableOffset=&geometryPrecision=&outSR=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&\
        groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&returnDistinctValues=false&\
        resultOffset=&resultRecordCount=&queryByDistance=&returnExtentsOnly=false&datumTransformation=&parameterValues=&\
        rangeValues=&f=pjson".format(where)
        fsURL = URL + query
        fs = arcpy.FeatureSet()
        fs.load(fsURL)
        if long(arcpy.GetCount_management(fs)[0])>0 :
            outname = '{0}_{1}_{2}'.format(name,i,j)
            arcpy.CopyFeatures_management(fs, outname)
            print(outname)
        else:
            break
            print('Downloaded up to {} features'.format(i))
    arcpy.Merge_management(arcpy.ListFeatureClasses('{}_*'.format(name)),output='{}'.format(name))
    for fc in arcpy.ListFeatureClasses('{}_*'.format(name)):
        arcpy.Delete_management(fc)

APIdownload("https://hazards.fema.gov/gis/nfhl/rest/services/PrelimPending/Prelim_NFHL/MapServer/24/query",
            'S_Fld_Haz_Ar', 100)