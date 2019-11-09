import arcpy
import imp
from arcpy.sa import *
import numpy as np
from os import *
from SpatialJoinLines_LargestOverlap import *

rootdir = 'C:/Mathis/ICSL/stormwater/'
resdir  = path.join(rootdir,'results')
arcpy.env.workspace = resdir
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

wa_wtshd = path.join(rootdir, 'data/NOAA_ERMA_20180601/PSwatersheds/shapefile_4047.shp')
USland = path.join(rootdir, 'data/USA_adm_gdb/USA_adm.gdb/USA_adm0')
OSMroads = path.join(rootdir, 'data/OSM_WA_20180601/gis_osm_roads_free_1.shp')
citylims = path.join(rootdir, 'data/WSDOT_GIS_20180508/CityLimits/CityLimits.gdb/CityLimits')
counties = path.join(rootdir, 'data/TIGER2017/tl_2018_us_county/tl_2018_us_county.shp')
NLCD_imp = os.path.join(rootdir, 'data/NLCD_2016_Impervious_L48_20190405.img') #Based on 2016 dara
ref_cs =


#Create gdb
gdb=path.join(rootdir,'results/PSOSM.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir, 'results'), 'PSOSM.gdb')

#Output variables
PSOSM='PSwtshd_OSMroads.shp'
PSOSM_all='PSwtshd_OSMroads_all.shp'

#######################################################################################################################
#Subset puget sound watersheds
arcpy.MakeFeatureLayer_management(wa_wtshd, out_layer='wtshd_lyr')
arcpy.SelectLayerByAttribute_management('wtshd_lyr', 'NEW_SELECTION', "PUGETSOUND = 'Y'")
arcpy.CopyFeatures_management('wtshd_lyr', 'PSwtshd.shp')
arcpy.Delete_management('wtshd_lyr')

#Dissolve Puget Sound watersheds into one polygon
arcpy.Dissolve_management('PSwtshd.shp', 'PSwtshd_dissolve.shp')

#Erase non-terrestrial areas from Puget Sound polygon
#proj = arcpy.Describe(wa_wtshd).spatialReference
arcpy.Project_management(USland, 'USA_adm_proj.shp', out_coor_system=UTM10)
arcpy.env.extent = arcpy.Describe(wa_wtshd).extent
arcpy.Buffer_analysis('USA_adm_proj.shp', out_feature_class='USbuf.shp', buffer_distance_or_field='100 meters')
arcpy.Intersect_analysis(in_features=['PSwtshd_dissolve.shp',USland], out_feature_class='PSwtshd_extrude')

#Intersect Puget Sound with road network from OSM
arcpy.MakeFeatureLayer_management(OSMroads, 'OSMroads_lyr')
#Do not include service streets as generally not enough through-traffic to be on Bing map, also excludes tracks as mainly used for forestry and agricultural purpose
sel = "{} IN ('motorway','motorway_link','living_street','primary','primary_link','residential','secondary','secondary_link'," \
      "'tertiary','tertiary_link','trunk','trunk_link','unclassified','unknown')".format('"fclass"')
arcpy.SelectLayerByAttribute_management('OSMroads_lyr', 'NEW_SELECTION', sel)
arcpy.Intersect_analysis(['OSMroads_lyr', 'PSwtshd_dissolve.shp'],out_feature_class=PSOSM)
arcpy.Dissolve_management(PSOSM,'PSwtshd_OSMroads_dissolve.shp')
arcpy.Delete_management('OSMroads_lyr')

#Select cities that intersect PS watershed
arcpy.MakeFeatureLayer_management(citylims, 'citylims_lyr')
arcpy.SelectLayerByLocation_management('citylims_lyr', 'INTERSECT', 'PSwtshd_dissolve.shp', selection_type='NEW_SELECTION')
arcpy.CopyFeatures_management('citylims_lyr', 'PScitylimits.shp')
arcpy.Delete_management('citylims_lyr')

#Select counties that intersect PS watershed
arcpy.MakeFeatureLayer_management(counties, 'counties_lyr')
arcpy.SelectLayerByLocation_management('counties_lyr', 'INTERSECT', 'PScitylimits.shp', selection_type='NEW_SELECTION')
arcpy.SelectLayerByAttribute_management('counties_lyr', 'ADD_TO_SELECTION', "COUNTYNS='01531927'") #Add Lewis county
arcpy.CopyFeatures_management('counties_lyr', 'PScounties.shp')
arcpy.Delete_management('counties_lyr')