#Author: Mathis Messager
#Date created: 2018/07/17

#http://docs.python-guide.org/en/latest/dev/virtualenvs/
#https://github.com/openaddresses/pyesridump

import arcpy
from os import *
import itertools
import glob
arcpy.env.overwriteOutput = True

outdir = "C:/Mathis/ICSL/stormwater/data/CitySeattle_20180717"
arcpy.env.workspace = outdir

def APIdownload(URL, name):
    IDrange = range(0,100000, 1000)
    for i,j in itertools.izip(IDrange, IDrange[1:]):
        where = "OBJECTID>{0} AND OBJECTID<{1}".format(i,j)
        #&geometryType=esriGeometryPoint
        query = "?where={}&returnGeometry=true&f=json".format(where)
        fsURL = URL + query
        fs = arcpy.FeatureSet()
        fs.load(fsURL)
        if long(arcpy.GetCount_management(fs)[0])>0 :
            name = '{0}_{1}_{2}'.format(name,i,j)
            arcpy.CopyFeatures_management(fs, name)
            print(name)
        else:
            break
            print('Downloaded up to {} features'.format(i))
    arcpy.Merge_management(arcpy.ListFeatureClasses('{}_*.shp'.format(name)),output='{}.shp'.format(name))
    for fc in arcpy.ListFeatureClasses('{}_*.shp'.format(name)):
        arcpy.Delete_management(fc)

APIdownload("https://gisrevprxy.seattle.gov/arcgis/rest/services/ext/WM_CityGISLayers/MapServer/66/query", 'creekwtshd')
APIdownload("https://gisrevprxy.seattle.gov/arcgis/rest/services/ext/WM_CityGISLayers/MapServer/56/query", 'urbanwtshd')



