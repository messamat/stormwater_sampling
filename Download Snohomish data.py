#Author: Mathis Messager
#Date created: 2018/06/02

#http://docs.python-guide.org/en/latest/dev/virtualenvs/
#https://github.com/openaddresses/pyesridump

import arcpy
from os import *
import itertools
import glob
arcpy.env.overwriteOutput = True

#Download road data for snohomish county
baseURL = "http://gismaps.snoco.org/snocogis/rest/services/transportation/transportation/MapServer/11/query"
outdir = "C:/Mathis/ICSL/stormwater/data/Snohomish_20180602"
arcpy.env.workspace = outdir
name='snocoroads'
IDrange = range(0,100000, 1000)
for i,j in itertools.izip(IDrange, IDrange[1:]):
    where = "OBJECTID>{0} AND OBJECTID<{1}".format(i,j)
    #&geometryType=esriGeometryPoint
    query = "?where={}&returnGeometry=true&f=json".format(where)
    fsURL = baseURL + query
    fs = arcpy.FeatureSet()
    fs.load(fsURL)
    if long(arcpy.GetCount_management(fs)[0])>0 :
        name = '{0}_{1}_{2}'.format(name,i,j)
        arcpy.CopyFeatures_management(fs, name)
        print(name)
    else:
        break
        print('Downloaded up to {} features'.format(i))

arcpy.Merge_management(arcpy.ListFeatureClasses('snocoroads_*.shp'),output='snocoroads.shp')
for fc in arcpy.ListFeatureClasses('snocoroads_*.shp'):
    arcpy.Delete_management(fc)
#

#Download
baseURL = "http://gismaps.snoco.org/snocogis/rest/services/transportation/transportation_infrastructure/MapServer/19/query"
outdir = "C:/Mathis/ICSL/stormwater/data/Snohomish_20180602"
arcpy.env.workspace = outdir
name='snocorow'
IDrange = range(0,100000, 1000)
for i,j in itertools.izip(IDrange, IDrange[1:]):
    where = "OBJECTID>{0} AND OBJECTID<{1}".format(i,j)
    #&geometryType=esriGeometryPoint
    query = "?where={}&returnGeometry=true&f=json".format(where)
    fsURL = baseURL + query
    fs = arcpy.FeatureSet()
    fs.load(fsURL)
    if long(arcpy.GetCount_management(fs)[0])>0 :
        name = '{0}_{1}_{2}'.format(name,i,j)
        arcpy.CopyFeatures_management(fs, name)
        print(name)
    else:
        break
        print('Downloaded up to {} features'.format(i))





################### DUMP ##########################
""""#In terminal
pip install virtualenv
virtualenv esridump
cd .\esridump
.\Scripts\activate
pip install esridump

esri2geojson http://gismaps.snoco.org/snocogis/rest/services/transportation/transportation/MapServer/11

import json
from esridump.dumper import EsriDumper
d = EsriDumper('http://example.com/arcgis/rest/services/Layer/MapServer/1')
# Iterate over each feature
for feature in d:
    print(json.dumps(feature))
d = EsriDumper('http://example.com/arcgis/rest/services/Layer/MapServer/2')
# Or get all features in one list
all_features = list(d)


#Then
from osgeo import ogr
points='C:/Users/install/esridump/traffic.geojson'
pointogr = ogr.CreateGeometryFromJson(points)
points

"""
#####################