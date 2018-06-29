import arcpy
from arcpy.sa import *
import numpy as np
from os import *

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
counties = path.join(rootdir, 'data/TIGER2017/cb_2017_us_county_500k/cb_2017_us_county_500k.shp')

kernel = NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernell00.txt') #UPDATE
UTM10 = arcpy.SpatialReference(26910)

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
PSOSM='PSwtshd_OSMroads.shp'
arcpy.Intersect_analysis(['OSMroads_lyr', 'PSwtshd_dissolve.shp'],out_feature_class=PSOSM)
arcpy.Dissolve_management(PSOSM,'PSwtshd_OSMroads_dissolve.shp')
arcpy.Delete_management('OSMroads_lyr')

#Select cities that intersect PS watershed
citylims = path.join(rootdir, 'data/WSDOT_GIS_20180508/CityLimits/CityLimits.gdb/CityLimits')
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

#Cut out OSM roads within King County boundaries and KC metro roads outside King County boundaries
arcpy.MakeFeatureLayer_management(counties, 'counties_lyr')
arcpy.SelectLayerByAttribute_management('counties_lyr', 'NEW_SELECTION', "COUNTYNS='01531933'")
arcpy.Erase_analysis(PSOSM, 'counties_lyr','PSwtshd_OSMroads_noKing.shp')
arcpy.Intersect_analysis([Kingroads,'counties_lyr'], 'Kingroads_cut.shp')

#-----------------------------------------------------------------------------------------------------------------------
# Join OSM and Pierce County data to improve interpolation of speed limits within road fclasses
#-----------------------------------------------------------------------------------------------------------------------
#Import Spatial Join Largest Overlap function from ESRI: http://www.arcgis.com/home/item.html?id=e9cccd343bf84916bda1910c31e5eab2
#spjoin = imp.load_source('module.name',  path.join(rootdir,'bin/SpatialJoinLargestOverlap/Scripts/SpatialJoinLargestOverlap.py'))
# [f.name for f in arcpy.ListFields(Pierceroads)]
target_features=path.join(resdir,'PSwtshd_OSMroads_noKing.shp')
join_features=Pierceroads
out_fc=path.join(resdir,'PSwtshd_OSMroads_noKing_Piercejoin.shp')
keep_all=True
fields_select = ['RoadNumber','RoadName','FFC','FFCDesc','ADTSource','ADT','ADTYear','SpeedLimit']
minoverlap=0.75

#Split OSM and county road lines at intersections
arcpy.FeatureToLine_management(target_features, 'target_split.shp')
arcpy.FeatureToLine_management(join_features, 'joinfeat_split.shp')
#Bufferize both road datasets
arcpy.Buffer_analysis('target_split.shp', 'target_buf.shp','10 meters', method='GEODESIC')
arcpy.Buffer_analysis('joinfeat_split.shp', 'joinfeat_buf.shp','10 meters', method='GEODESIC')
#Get buffer area for target feature
arcpy.AddGeometryAttributes_management('target_buf.shp', 'AREA_GEODESIC', Area_Unit='SQUARE_METERS')

#Spatial join with largest overlap
# Calculate intersection between Target Feature and Join Features
intersect="buf_inters.shp"
arcpy.Intersect_analysis(['target_buf.shp', 'joinfeat_buf.shp'], path.join(resdir,intersect), join_attributes='ONLY_FID')
arcpy.AddGeometryAttributes_management(intersect, 'AREA_GEODESIC', Area_Unit='SQUARE_METERS')
arcpy.JoinField_management(intersect, 'FID_target', join_table='target_buf.shp', join_field='FID', fields='AREA_GEO')

[f.name for f in arcpy.ListFields(intersect)]
# Find which Join Feature has the largest overlap with each Target Feature
fields = ["FID_{0}".format(path.splitext(path.basename(target_features))[0]),
          "FID_{0}".format(path.splitext(path.basename(join_features))[0]),
          "AREA_GEO"]
overlap_dict = {}
with arcpy.da.SearchCursor(intersect, fields) as scur:
    for row in scur:
        try:
            if row[2] > overlap_dict[row[0]][1]:
                overlap_dict[row[0]] = [row[1], row[2]]
        except:
            overlap_dict[row[0]] = [row[1], row[2]]

# Copy the target features and write the largest overlap join feature ID to each record
# Set up all fields from the target features + ORIG_FID
fieldmappings = arcpy.FieldMappings()
fieldmappings.addTable(target_features)
fieldmap = arcpy.FieldMap()
fieldmap.addInputField(target_features, arcpy.Describe(target_features).OIDFieldName)
fld = fieldmap.outputField
fld.type, fld.name, fld.aliasName = "LONG", "ORIG_FID", "ORIG_FID"
fieldmap.outputField = fld
fieldmappings.addFieldMap(fieldmap)
# Perform the copy
arcpy.conversion.FeatureClassToFeatureClass(target_features, path.dirname(out_fc), path.basename(out_fc), "", fieldmappings)
# Add a new field JOIN_FID to contain the fid of the join feature with the largest overlap
arcpy.AddField_management(out_fc, "JOIN_FID", "LONG")
# Calculate the JOIN_FID field
with arcpy.da.UpdateCursor(out_fc, ["ORIG_FID", "JOIN_FID"]) as ucur:
    for row in ucur:
        try:
            row[1] = overlap_dict[row[0]][0]
            ucur.updateRow(row)
        except:
            if not keep_all:
                ucur.deleteRow()
# Join all attributes from the join features to the output
joinfields = [x.name for x in arcpy.ListFields(join_features) if not x.required and x.name in fields_select]
arcpy.JoinField_management(out_fc, "JOIN_FID", join_features, arcpy.Describe(join_features).OIDFieldName, joinfields)

# Need to know the Target Features shape type, to know to read the SHAPE_AREA oR SHAPE_LENGTH property
#geom = "AREA" if arcpy.Describe(target_features).shapeType.lower() == "polygon" and arcpy.Describe(join_features).shapeType.lower() == "polygon" else "LENGTH"

##########################################################################################################################################
#Create heatmap based on roads functional class for all Puget Sound OSM roads
##########################################################################################################################################
#Create gdb
gdb=path.join(rootdir,'results/PSOSM')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir, 'results'), 'PSOSM')
arcpy.env.workspace = gdb

#Select OSM roads, but this map with service roads, as enough through traffic to potentially have some impact
arcpy.MakeFeatureLayer_management(OSMroads, 'OSMroads_lyr')
np.unique([row[0] for row in arcpy.da.SearchCursor('OSMroads_lyr', ['fclass'])])
sel = "{} IN ('motorway','motorway_link','living_street','primary','primary_link','residential','secondary','secondary_link'," \
      "'tertiary','tertiary_link','trunk','trunk_link','unclassified','unknown', 'service')".format('"fclass"')
arcpy.SelectLayerByAttribute_management('OSMroads_lyr', 'NEW_SELECTION', sel)

PSOSM_all='PSwtshd_OSMroads_all.shp'
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

#Create heatmap
res = arcpy.GetRasterProperties_management(path.join(rootdir,'results/bing/180620_09_30_class_mlc.tif'), 'CELLSIZEX')
PSOSM_allproj = 'PSwtshd_OSMroads_all_proj.shp'
arcpy.Project_management(PSOSM_all, PSOSM_allproj, out_coor_system=UTM10)
arcpy.PolylineToRaster_conversion(PSOSM_allproj, value_field='fclassnum', out_rasterdataset='osmfclass', priority_field='fclassnum',cellsize=res)
heat_osmfclas = FocalStatistics(path.join(gdb,'osmfclass'), neighborhood=NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernel100.txt'),
                                statistics_type='SUM', ignore_nodata='DATA') #It seems that full paths and using a raster within a file geodatabase is required to run function
heat_osmfclas.save('heat_osmfclas')
heat_osmfclas_int = Int(Raster('heat_osmfclas'+0.5))
heat_osmfclas_int.save('heat_osmfclas_int')