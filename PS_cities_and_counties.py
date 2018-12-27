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
Kingroads = path.join(rootdir, 'data/King_201806\Metro_Transportation_Network_TNET_in_King_County_for_Car_Mode__trans_network_car_line\Metro_Transportation_Network_TNET_in_King_County_for_Car_Mode__trans_network_car_line.shp')
Pierceroads = path.join(rootdir, 'data/Pierce_20180611/Mobility_Data/Mobility_Data.shp')
counties = path.join(rootdir, 'data/TIGER2017/tl_2018_us_county/tl_2018_us_county.shp')
citylims = path.join(rootdir, 'data/WSDOT_GIS_20180508/CityLimits/CityLimits.gdb/CityLimits')
traffic_wsdot = path.join(rootdir, 'data/WSDOT_TPTTraffic_20180508/2016_TrafficCounts/2016TrafficCounts.gdb/TrafficCounts2016')

kernel = NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernell00.txt') #UPDATE
UTM10 = arcpy.SpatialReference(26910)

#Create gdb
gdb=path.join(rootdir,'results/PSOSM.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir, 'results'), 'PSOSM.gdb')

#Output variables
PSOSM='PSwtshd_OSMroads.shp'
PSOSM_all='PSwtshd_OSMroads_all.shp'
OSMPierce = 'PSwtshd_OSMroads_Pierce.shp'
intersect= path.join(gdb, "buf_inters")
intersect_stats = path.join(gdb, "buf_inters_stats")
OSMPierce_datajoin = path.join(gdb, 'OSMPierce_datajoin')
OSMWSDOT_datajoin = path.join(gdb, 'OSM_WSDOT_joinstats')

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
arcpy.CopyFeatures_management('counties_lyr', 'PScounties')
arcpy.Delete_management('counties_lyr')

##########################################################################################################################################
#Create heatmap based on roads functional class for all Puget Sound OSM roads
##########################################################################################################################################
arcpy.env.workspace = gdb

#Select OSM roads, but this map with service roads, as enough through traffic to potentially have some impact
arcpy.MakeFeatureLayer_management(OSMroads, 'OSMroads_lyr')
np.unique([row[0] for row in arcpy.da.SearchCursor('OSMroads_lyr', ['fclass'])])
sel = "{} IN ('motorway','motorway_link','living_street','primary','primary_link','residential','secondary','secondary_link'," \
      "'tertiary','tertiary_link','trunk','trunk_link','unclassified','unknown', 'service')".format('"fclass"')
arcpy.SelectLayerByAttribute_management('OSMroads_lyr', 'NEW_SELECTION', sel)
arcpy.Intersect_analysis(['OSMroads_lyr', 'PSwtshd_dissolve.shp'],out_feature_class=PSOSM_all)
arcpy.Delete_management('OSMroads_lyr')

#Convert OSM functional categories to numbers
numdic = {'unknown':1,'service':1, 'residential':2, 'unclassified':3,'tertiary_link':3.5,'tertiary':4, 'secondary_link':4.5,
 'secondary':5, 'primary_link':5.5, 'primary':6, 'trunk_link':6.5, 'trunk':7,'motorway_link':7.5,'motorway':8}
arcpy.AddField_management(PSOSM_all, 'fclassnum', 'FLOAT')
with arcpy.da.UpdateCursor(PSOSM_all, ['fclass','fclassnum']) as cursor:
    for row in cursor:
        if row[0] in numdic.keys():
            row[1] = numdic[row[0]]
        else:
            row[1]=0
        cursor.updateRow(row)


# Join OSM and Pierce County + WSDOT traffic counts data to improve interpolation of speed limits
# and traffic volume within road fclasses

#Subselect OSM roads for Pierce County
arcpy.MakeFeatureLayer_management(counties, 'counties_lyr')
arcpy.SelectLayerByAttribute_management('counties_lyr', 'NEW_SELECTION', "COUNTYNS='01529159'")
arcpy.Clip_analysis(PSOSM_all, 'counties_lyr', OSMPierce)

SpatialJoinLines_LargestOverlap(target_features= OSMPierce, join_features=Pierceroads, out_fc = OSMPierce_datajoin,
                                outgdb=gdb, bufsize='10 meters', keep_all=True,
                                fields_select=['RoadNumber', 'RoadName', 'FFC', 'FFCDesc', 'ADTSource', 'ADT',
                                               'ADTYear', 'SpeedLimit'])
#Join OSM with WSDOT traffic counts
arcpy.SpatialJoin_analysis(traffic_wsdot, PSOSM_all, path.join(gdb, 'OSM_WSDOT_join'), 'JOIN_ONE_TO_ONE', 'KEEP_COMMON',
                           match_option='CLOSEST_GEODESIC', search_radius='20 meters', distance_field_name='joindist')
arcpy.Statistics_analysis(path.join(gdb, 'OSM_WSDOT_join'),OSMWSDOT_datajoin ,
                          statistics_fields= [['AADT', 'MEAN'],['DirectionOfTravel', 'FIRST'],['fclass', 'FIRST']],
                          case_field = 'osm_id')