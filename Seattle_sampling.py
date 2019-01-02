#Author: Mathis Messager
#Purpose: create a heuristic model of traffic-related heavy metal pollution to select maple trees on which to sample moss
# in the city of Seattle, WA

import arcpy
from arcpy.sa import *
from os import *
import re
import numpy as np
from collections import defaultdict
import pandas as pd
import numpy as np
import time

#Custom modules
from SpatialJoinLines_LargestOverlap import *
from heatmap_custom import *
from GTFStoSHP import *
from explode_overlapping import *

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

#Set up paths
rootdir = 'F:/Levin_Lab/stormwater/'

roads = path.join(rootdir, 'data/CitySeattle_20180601/Seattle_Streets/Seattle_Streets.shp')
traffic_seattle = path.join(rootdir, 'data/CitySeattle_20180601/2016_Traffic_Flow_Counts/2016_Traffic_Flow_Counts.shp')
traffic_wsdot = path.join(rootdir, 'data/WSDOT_TPTTraffic_20180508/2016_TrafficCounts/2016TrafficCounts.gdb/TrafficCounts2016')
trees = path.join(rootdir, 'data/CitySeattle_20180601/Trees/Trees.shp')
zoning = path.join(rootdir, 'data/CitySeattle_20180626/City_of_Seattle_Zoning/WGS84/City_of_Seattle_Zoning.shp')
censustract = path.join(rootdir, 'data/TIGER2017/Profile-County_Tract/Profile-County_Tract.gdb/Tract_2010Census_DP1')
heat_bing = path.join(rootdir, 'results/bing/bingmean1806_Seattle_heat.tif')
NLCD_reclass = path.join(rootdir, 'results/LU.gdb/NLCD_reclass_final')
NLCD_imp = path.join(rootdir, 'data/nlcd_2011_impervious_2011_edition_2014_10_10/nlcd_2011_impervious_2011_edition_2014_10_10.img')
PSwatershed = path.join(rootdir, 'results/PSwtshd_dissolve.shp')
cities = path.join(rootdir, 'results/PScitylimits.shp')
counties = path.join(rootdir, 'data/TIGER2017/tl_2018_us_county/tl_2018_us_county.shp')

template_ras = path.join(rootdir,'results/bing/181204_02_00_class_mlc.tif')
res = arcpy.GetRasterProperties_management(template_ras, 'CELLSIZEX')

PSgdb=path.join(rootdir,'results/PSOSM.gdb')
PSOSM_all= path.join(rootdir, 'results/PSwtshd_OSMroads_all.shp')
OSMPierce_datajoin = path.join(PSgdb, 'OSMPierce_datajoin')
OSMWSDOT_datajoin = path.join(PSgdb, 'OSM_WSDOT_joinstats')

PSgtfs = os.path.join(rootdir, 'data\SoundTransit_201812\gtfs_puget_sound_consolidated.zip')
Binggdb = os.path.join(rootdir, 'results/bing/postprocess.gdb')
#NED data paths to be updated when transfer to laptop
hardrive = 'D:/Processing'
NED19proj = os.path.join(hardrive, 'ned19_psproj')
NED13proj = os.path.join(rootdir, 'results/ned13_psproj')
NED19smooth = os.path.join(hardrive, 'ned19_smooth')
NED13smooth = os.path.join(hardrive, 'ned13_smooth')
rangetab19 = path.join(PSgdb, 'OSM_elv19range')
rangetab13 = path.join(PSgdb, 'OSM_elv13range')
rangetab19_smooth = path.join(PSgdb, 'OSM_elv19range_smooth')
rangetab13_smooth = path.join(PSgdb, 'OSM_elv13range_smooth')
#########

gdb = path.join(rootdir,'results/Seattle_sampling.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir,'results'), 'Seattle_sampling.gdb')
arcpy.env.workspace = gdb

#New variables
roadstraffic = 'Seattle_roadstraffic'
roadstraffic_avg =roadstraffic+'_AADT'
NLCD_reclass_PS = path.join(rootdir, 'results/NLCD_reclass_final_PS.tif')
NLCD_imp_PS = path.join(rootdir, 'results/nlcd_imp_ps')
UTM10 = arcpy.SpatialReference(26910)
OSMSeattle = path.join(PSgdb, 'PSwtshd_OSMroads_Seattle')
OSMSeattle_datajoin = path.join(PSgdb, 'OSMSeattle_datajoin')
PSOSM_allproj = path.join(PSgdb, 'PSwtshd_OSMroads_all_proj')
OSM_AADT = path.join(PSgdb,'OSM_AADT')
PSOSMras = path.join(PSgdb, 'PSwtshd_OSMroads_all_projras')
PSOSM_elv = path.join(PSgdb, 'PSwtshd_OSMroads_all_proj_elv')
OSM_gradient = path.join(PSgdb, 'OSM_gradient')

PStransit = os.path.join(rootdir, 'results/transit.gdb/PStransit')
PStransitbus = PStransit + '_busroutes'
PStransitbus_proj = PStransit + '_busroutes_proj'
PStransitbus_splitdiss = PStransitbus_proj + '_splitv_diss'
PStransitras = os.path.join(rootdir, 'results/transit.gdb/PStransit_ras')

########################################################################################################################
# GET TRAFFIC COUNT FOR EVERY ROAD SEGMENT IN SEATTLE
#   Average Annual Daily Traffic (AADT) is available for most arterial streets in the city of Seattle (here for 2016).
#   However, there is no common street identifier between the AADT dataset and the broader street dataset
#   of the city of Seattle. Therefore, spatially join the AADT dataset and the City of Seattle general street dataset.
#   Note that:
#   1. AADT segments are larger than City of Seattle streets, hence no need to split the street dataset.
#   2. Here average week day traffic is used as it has a fuller coverage of streets.
########################################################################################################################

#Check for common fields between Seattle streets dataset and traffic_count dataset
roads_fields = [f.name for f in arcpy.ListFields(roads)]
traffic_fields = [f2.name for f2 in arcpy.ListFields(traffic_seattle)]
list(set(roads_fields).intersection(traffic_fields))

#Project datasets
print('Project datasets')
arcpy.Project_management(roads, 'roads_proj', out_coor_system=UTM10)
arcpy.Project_management(traffic_seattle, 'traffic_proj', out_coor_system=UTM10)

#Join Seattle traffic data to street data
print('Join Seattle traffic data to street data')
arcpy.SpatialJoin_analysis('roads_proj', 'traffic_proj', roadstraffic, join_operation='JOIN_ONE_TO_ONE',
                           join_type='KEEP_ALL', match_option='HAVE_THEIR_CENTER_IN', search_radius='2 meters')

#Get rid of erroneous matches at intersections between highways and arterials
print('Get rid of erroneous matches')
traffic_fields = [f.name for f in arcpy.ListFields('Seattle_roadstraffic')][39:]
traffic_fields.insert(0,'ARTDESCRIP')
traffic_fields.insert(1,'FIRST_ARTD')
with arcpy.da.UpdateCursor(roadstraffic,traffic_fields) as cursor:
    for row in cursor:
        if row[1] is not None and row[0] != row[1]:
            row[2:len(traffic_fields)] = np.repeat(None, len(traffic_fields)-2)
        cursor.updateRow(row)

#Create unique ID field to help in joining
print('Create CUSTOM_ID')
arcpy.AddField_management(roadstraffic, 'CUSTOM_ID', 'LONG')
with arcpy.da.UpdateCursor(roadstraffic, ['CUSTOM_ID']) as cursor:
    x=0
    for row in cursor:
        row[0]=x
        cursor.updateRow(row)
        x=x+1

#-----------------------------------------------------------------------------------------------------------------------
# JOIN DATA SOURCES
#   The Seattle AADT dataset does not include data for state highways. The Washington State Department of Transportation,
#   on the other hand, has point AADT data for most highways and several arterials.
#   Here, we spatially join WSDOT point data to the closest Seattle street segment dataset
#-----------------------------------------------------------------------------------------------------------------------

#Join WSDOT traffic data to Seattle streets
print('Join WSDOT traffic data to Seattle streets')
arcpy.SpatialJoin_analysis(traffic_wsdot, roadstraffic, 'WSDOT_streets_join', 'JOIN_ONE_TO_ONE', 'KEEP_COMMON',
                           match_option='CLOSEST_GEODESIC', search_radius='50 meters', distance_field_name='joindist')
arcpy.Dissolve_management('WSDOT_streets_join', 'WSDOT_streets_join_diss', dissolve_field='CUSTOM_ID',
                          statistics_fields=[['AADT', 'MEAN'],['joindist', 'MEAN'],['Join_Count','SUM']])
arcpy.JoinField_management(roadstraffic, 'CUSTOM_ID', 'WSDOT_streets_join_diss', 'CUSTOM_ID', ['MEAN_joindist','MEAN_AADT', 'SUM_Join_Count'])
#Not a big deal if the match isn't perfect as the final output will be rasterized

#Average between WSDOT and Seattle estimates (as some roads have AADT estimates from both sources)
print('Compute average between WSDOT and Seattle estimates')
arcpy.AddField_management(roadstraffic, 'AADT_avg', 'DOUBLE')
[f.name for f in arcpy.ListFields(roadstraffic)]
arcpy.CalculateField_management(roadstraffic, 'AADT_avg', '!MEAN_AADT!', expression_type='PYTHON')
arcpy.MakeFeatureLayer_management(roadstraffic, 'roadstraffic_lyr')
arcpy.SelectLayerByAttribute_management('roadstraffic_lyr', 'NEW_SELECTION', "NOT COUNTAAWDT IS NULL")
arcpy.CalculateField_management('roadstraffic_lyr', 'AADT_avg', '(float(!COUNTAAWDT!)+!MEAN_AADT!)/2', expression_type='PYTHON')
arcpy.SelectLayerByAttribute_management('roadstraffic_lyr', 'NEW_SELECTION', "MEAN_AADT IS NULL AND NOT COUNTAAWDT IS NULL")
arcpy.CalculateField_management('roadstraffic_lyr', 'AADT_avg', expression='float(!COUNTAAWDT!)', expression_type='PYTHON')
arcpy.SelectLayerByAttribute_management('roadstraffic_lyr', 'CLEAR_SELECTION')
arcpy.CopyFeatures_management('roadstraffic_lyr', roadstraffic_avg)
arcpy.Delete_management('roadstraffic_lyr')

#-----------------------------------------------------------------------------------------------------------------------
# INTERPOLATE
#   After joining the two sources, AADT values are interpolated for missing segments based on linear inverse distance
#   weighting from road segments with the same name.
#-----------------------------------------------------------------------------------------------------------------------

#Interpolate for designated lines
print('Interpolate AADT values along roads based on inverse distance weighting along lines')
arcpy.MakeFeatureLayer_management(roadstraffic_avg, 'roadstraffic_lyr')
arcpy.SelectLayerByAttribute_management('roadstraffic_lyr', 'NEW_SELECTION',
                                        "NOT (ARTDESCRIP = ' ' OR ARTDESCRIP = 'Not Designated')")
routes = 'roadstraffic_route'
arcpy.CreateRoutes_lr('roadstraffic_lyr', route_id_field='STNAME_ORD', out_feature_class=routes, measure_source='LENGTH',
                      ignore_gaps='NO_IGNORE')
arcpy.LocateFeaturesAlongRoutes_lr('roadstraffic_lyr', routes, route_id_field = 'STNAME_ORD', out_table= 'routes_loc')
#For each line with a designated type and no AADT, get the AADT of the closest segment on both sides on the same street name
arcpy.AddField_management('routes_loc', field_name='AADT_interp', field_type='double')
interpdicplus = defaultdict(list)
interpdicminus = defaultdict(list)
[f.name for f in arcpy.ListFields('routes_loc')]
with arcpy.da.UpdateCursor('routes_loc', ['CUSTOM_ID','STNAME_ORD','FMEAS', 'ARTDESCRIP', 'AADT_avg','AADT_interp']) as cursorOuter:
    for rowOuter in cursorOuter:
        if rowOuter[4] is None: #If AADT is null
            print(rowOuter[0])
            #For each segment of artierial street
            with arcpy.da.SearchCursor('routes_loc',['CUSTOM_ID', 'STNAME_ORD', 'FMEAS', 'ARTDESCRIP', 'AADT_avg']) as cursorInner:
                for rowInner in cursorInner:
                    #Loop over segments of the same street that have a non-zero AADT value
                    if rowOuter[1] == rowInner[1] and rowInner[4] > 0:
                        #Compute the difference in position along the street of each segment to the reference segment
                        fdif = rowOuter[2] - rowInner[2]
                        if fdif>0:
                            try:
                                if fdif < interpdicplus[rowOuter[0]][0]: #If segment is closer to reference segment than previous closest segment
                                    interpdicplus[rowOuter[0]] = [fdif, rowInner[4]] #Then write [distance, AADT]
                            except:
                                interpdicplus[rowOuter[0]] = [fdif, rowInner[4]] #If no record was in dic for that segment, then write [distance, AADT]
                        if fdif<0:
                            try:
                                if fdif > interpdicminus[rowOuter[0]][0]:  # If segment is closer to reference segment than previous closest segment
                                    interpdicminus[rowOuter[0]] = [-fdif,rowInner[4]]  # Then write [distance, AADT]
                            except:
                                interpdicminus[rowOuter[0]] = [-fdif, rowInner[4]] #If no record was in dic for that segment, then write [distance, AADT]
            #Interpolate AADT
            #if type(interpdicplus) is list and type(interpdicminus) is list:
            if interpdicplus[rowOuter[0]] != [] and interpdicminus[rowOuter[0]] != []: #If segments with AADT both up and down
                # Compute inverse distance weighted average of AADT
                rowOuter[5] = (interpdicplus[rowOuter[0]][1]/interpdicplus[rowOuter[0]][0])+(interpdicminus[rowOuter[0]][1]/interpdicminus[rowOuter[0]][0])/\
                                (1/interpdicplus[rowOuter[0]][0]+1/interpdicminus[rowOuter[0]][0])
            else:
                try: #If not, assign the value of the nearest segment up or down the road
                    rowOuter[5] = interpdicplus[rowOuter[0]][1]
                except:
                    try:
                        rowOuter[5] = interpdicminus[rowOuter[0]][1]
                    except:
                        print('There is no AADT data for this road')
            cursorOuter.updateRow(rowOuter)

#-----------------------------------------------------------------------------------------------------------------------
#FILL IN AADT FOR REMAINING ROADS AND CORRECT ERRORS
#   For freeway ramps, get the city-wide median AADT for freeway ramps
#   For arterial roads with no AADT estimates, get the city-wide median AADT for their given road type (minor vs major arterial, etc.
#   For all non-arterial roads with no estimates, assign AADT of 1000
#-----------------------------------------------------------------------------------------------------------------------

#Compute average AADT on freeway ramps
print('Compute averafe AADT on freeway ramps')
arcpy.MakeFeatureLayer_management(roadstraffic_avg, 'roadstraffic_lyr')
arcpy.SelectLayerByAttribute_management('roadstraffic_lyr', 'NEW_SELECTION', where_clause="STNAME_ORD LIKE '% RP'")
RP_ADDT=np.median([int(round(row[0])) for row in arcpy.da.SearchCursor('roadstraffic_lyr', 'AADT_avg') if row[0] is not None])
arcpy.Delete_management('roadstraffic_lyr')

#Compute median AADT by road type
print('Compute median AADT by road type')
arcpy.MakeFeatureLayer_management(roadstraffic_avg, 'roadstraffic_lyr')
arcpy.SelectLayerByAttribute_management('roadstraffic_lyr', 'NEW_SELECTION','NOT AADT_avg IS NULL')
statdic = defaultdict(list)
for row in arcpy.da.SearchCursor('roadstraffic_lyr', ['ARTDESCRIP','AADT_avg']): #Get dictionary of key: road type, value: AADT
    statdic[row[0]].append(int(round(row[1])))
for k, v in statdic.iteritems(): #For each road type, compute median AADT
    statdic[k] = np.median(v)

#Fill AADT interpolation for remaining arterial roads
print('Fill AADT interpolation for remaining arterial roads')
with arcpy.da.UpdateCursor('routes_loc', ['STNAME_ORD','ARTDESCRIP','AADT_avg','AADT_interp','CUSTOM_ID']) as cursor:
    for row in cursor:
        if row[2] is None and row[3] is None:
            if row[0] == 'I90 WB' or row[0] == 'I90 EXPRESS': #Check value of adjacent portions of I90 in that area
                print(row[0])
                row[3] = 15000
            if row[0][-3:] == ' RP': #Ramps get median AADT for all ramps that have data
                row[3]= RP_ADDT
            if row[4]==16130: #Portion of the Burke-Gilman trail was erroneously assigned 200,000 AADT
                row[3] = 0
            if row[3] is None: #Assign median AADT by road type
                if row[1] not in ['','Not Designated']:
                    row[3]=statdic[row[1]]
                else:
                    row[3]=1000
        if row[3] is None: #Transfer AADT_avg to AADT_interp
            row[3] = row[2]
        cursor.updateRow(row)

#For all non-arterial roads that do not already have an AADT value, assign 1000
arcpy.JoinField_management(roadstraffic_avg, 'CUSTOM_ID', 'routes_loc', 'CUSTOM_ID', 'AADT_interp')
with arcpy.da.UpdateCursor(roadstraffic_avg, ['ARTDESCRIP', 'AADT_avg', 'AADT_interp']) as cursor:
    for row in cursor:
        if row[1] is not None:
            row[2] = row[1]
        if row[1] is None and row[2] is None and row[0] in ['', 'Not Designated']:
            row[2] = 1000
        cursor.updateRow(row)

#Write out individual roads' attributes
arcpy.CopyRows_management(roadstraffic_avg, path.join(rootdir, 'results/Seattle_roads.dbf'))

########################################################################################################################
# PREPARE OTHER VARIABLES TO CREATE HEATMAPS: FUNCTIONAL-CLASS BASED AADT AND SPEED LIMIT, SLOPE, AND TRANSIT ROUTES
########################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------
#ASSIGN FUNCTIONAL CLASS-BASED MEDIANS OF AADT TO OSM FOR PUGET SOUND
#-----------------------------------------------------------------------------------------------------------------------
#Subselect OSM roads for Seattle
arcpy.MakeFeatureLayer_management(cities, 'cities_lyr')
arcpy.SelectLayerByAttribute_management('cities_lyr', selection_type='NEW_SELECTION',
                                        where_clause="CityFIPSLo ='5363000'")
arcpy.Clip_analysis(PSOSM_all, 'cities_lyr', OSMSeattle)

#Join Seattle AADT data to OSM Seattle
SpatialJoinLines_LargestOverlap(target_features=OSMSeattle, join_features=path.join(gdb, roadstraffic_avg),
                                out_fc=OSMSeattle_datajoin, outgdb=PSgdb, bufsize='10 meters', keep_all=True,
                                fields_select=['AADT_avg', 'SPEEDLIMIT', 'ARTDESCRIP'])

#Create a compound table for all data on traffic counts
#(only keep OSM street segments whose buffer overlapped at least 80% with original data & whose length > 5 m
# & post-2000 Pierce County data)
fclass_ADT = pd.DataFrame(
    [list(row[:3]) + [None, 'Seattle']
     for row in arcpy.da.SearchCursor(OSMSeattle_datajoin, ['osm_id', 'fclass', 'AADT_avg', 'intersper', 'LENGTH_GEO']) \
     if row[2] is not None and row[3] > 0.8 and row[4] > 10] +
    [list(row[:3])+ [int(row[3]), 'Pierce']
     for row in arcpy.da.SearchCursor(OSMPierce_datajoin, ['osm_id', 'fclass', 'ADT', 'ADTYear', 'intersper', 'LENGTH_GEO']) \
     if row[2] is not None and row[3] != u' ' and row[4] > 0.8 and row[5] > 10] +
    [list(row) + [None, 'WSDOT']
     for row in arcpy.da.SearchCursor(OSMWSDOT_datajoin, ['osm_id', 'FIRST_fclass', 'MEAN_AADT']) \
     if row[2] is not None],
    columns=['osm_id', 'fclass', 'ADT', 'year', 'agency'])
fclass_ADT = fclass_ADT[(fclass_ADT['year'] > 2000) | (fclass_ADT['year'].isnull())]
#Check segments that have multiple ADT counts but different values
duplis = fclass_ADT[(fclass_ADT.duplicated('osm_id', keep=False)) &
                    (~(fclass_ADT.duplicated(['osm_id','ADT'], keep=False)))].sort_values(by=['osm_id'])
#Keep most recent measurement for Pierce
fclass_ADT[fclass_ADT['agency'] == 'Pierce'] = fclass_ADT[fclass_ADT['agency'] == 'Pierce'].\
    sort_values('year', ascending=False).\
    drop_duplicates('osm_id').\
    sort_index()
#Average for each segment
fclass_ADTnodupli = fclass_ADT.groupby('osm_id').agg({'ADT' : ['mean'], 'fclass' : ['first', 'last']})
check = fclass_ADTnodupli[fclass_ADTnodupli['fclass','first'] != fclass_ADTnodupli['fclass','last']]
#Compute median ADT for each OSM functional class
fclass_ADTnodupli.columns = fclass_ADTnodupli.columns.droplevel()
fclass_ADTmedian = fclass_ADTnodupli.groupby('first').median()

#Convert OSM functional categories to numbers
#Service road ADT doesn't really make sense, as includes mostly alleys.
#Unclassified seems to be misused in OSM and often appear to stand for unknwon
#Change both to residential level
arcpy.AddField_management(PSOSM_all, 'fclassADT', 'LONG')
with arcpy.da.UpdateCursor(PSOSM_all, ['fclass','fclassADT']) as cursor:
    for row in cursor:
        if row[0] in fclass_ADTmedian.index:
            if row[0] in ['service','unclassified']:
                row[1] = int(fclass_ADTmedian.loc['residential'])
            else:
                row[1] = int(fclass_ADTmedian.loc[row[0]])
        else:
            row[1]=0
        cursor.updateRow(row)

#Project OSM data
arcpy.Project_management(PSOSM_all, PSOSM_allproj, out_coor_system=UTM10)

#-----------------------------------------------------------------------------------------------------------------------
# PREPARE TRANSIT DATA TO CREATE HEATMAP BASED ON BUS ROUTES
#-----------------------------------------------------------------------------------------------------------------------
#Convert GTFS to routes with number of trips per week on each line
GTFStoSHPweeklynumber(gtfs_dir= PSgtfs, out_gdb=os.path.dirname(PStransit), out_fc = os.path.basename(PStransit),
                      keep = True)
#Only keep buses with trips
arcpy.MakeFeatureLayer_management(PStransit + '_routes', 'PStransit_lyr',
                                  where_clause= '(route_type = 3) AND (MIN_service_len > 0) AND (SUM_adjustnum > 0)')
arcpy.CopyFeatures_management('PStransit_lyr', PStransitbus)
arcpy.Project_management(PStransitbus, PStransitbus_proj, UTM10)

#Create raster of weekly number of buses at the same resolution as bing data
# Convert weekly number of buses to integer
arcpy.AddField_management(PStransitbus_proj, 'adjustnum_int', 'SHORT')
arcpy.CalculateField_management(PStransitbus_proj, 'adjustnum_int',
                                expression='int(10*!SUM_adjustnum!+0.5)', expression_type='PYTHON')

arcpy.SplitLine_management(PStransitbus_proj, PStransitbus_proj + '_splitv')
arcpy.FindIdentical_management(PStransitbus_proj + '_splitv', "explFindID", "Shape")
arcpy.MakeFeatureLayer_management(PStransitbus_proj + '_splitv', "intlyr")
arcpy.AddJoin_management("intlyr", arcpy.Describe("intlyr").OIDfieldName, "explFindID", "IN_FID", "KEEP_ALL")
arcpy.Dissolve_management("intlyr", PStransitbus_splitdiss, dissolve_field='explFindID.FEAT_SEQ',
                          statistics_fields=[[os.path.split(PStransitbus_proj)[1] + '_splitv.adjustnum_int', 'SUM']])
arcpy.RepairGeometry_management(PStransitbus_splitdiss, delete_null = 'DELETE_NULL') #sometimes creates empty geom
tolerance = (2.0**0.5)*float(res.getOutput(0))/2 #length of half pixel diagonal
ExplodeOverlappingLines(PStransitbus_splitdiss, tolerance)

#For each set of non-overlapping lines, create its own raster
tilef = 'expl'
tilelist = list(set([row[0] for row in arcpy.da.SearchCursor(PStransitbus_splitdiss, [tilef])]))
outras_base = path.join(rootdir, 'results/transit.gdb/busnum_')
arcpy.env.snapRaster = template_ras
for tile in tilelist:
    outras = outras_base + str(tile)
    if not arcpy.Exists(outras):
        selexpr = '{0} = {1}'.format(tilef, tile)
        print(selexpr)
        arcpy.MakeFeatureLayer_management(PStransitbus_splitdiss, 'bus_lyr', where_clause= selexpr)
        arcpy.PolylineToRaster_conversion('bus_lyr', value_field='adjustnum_int', out_rasterdataset=outras, cellsize=res)

#Mosaic to new raster
arcpy.env.workspace = os.path.split(outras_base)[0]
transitras_tiles = arcpy.ListRasters('busnum_*')
arcpy.MosaicToNewRaster_management(transitras_tiles, arcpy.env.workspace, os.path.split(PStransitras)[1],
                                   pixel_type='32_BIT_UNSIGNED', number_of_bands= 1, mosaic_method = 'SUM')
for tile in transitras_tiles:
    print('Deleting {}...'.format(tile))
    arcpy.Delete_management(tile)
arcpy.ClearEnvironment('Workspace')

#-----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA ON ROAD GRADIENTS
#-----------------------------------------------------------------------------------------------------------------------
#Run on NED19 max-min for line itself
arcpy.PolylineToRaster_conversion(PSOSM_allproj, 'osm_id', PSOSMras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'fclassADT', cellsize = NED19proj)
ZonalStatisticsAsTable(PSOSMras, 'osm_id', NED19proj, out_table = rangetab19, statistics_type= 'RANGE', ignore_nodata='NODATA')


arcpy.AddField_management(rangetab19, 'osm_id', 'TEXT')
arcpy.CalculateField_management(rangetab19, 'osm_id', '!osm_id!', 'PYTHON')

#Run on NED13 max-min for line itself
ZonalStatisticsAsTable(PSOSMras, 'Value', NED13proj, out_table = rangetab13, statistics_type= 'RANGE', ignore_nodata='NODATA')
arcpy.AddField_management(rangetab13, 'osm_id', 'TEXT')
arcpy.CalculateField_management(rangetab13, 'osm_id', '!Value!', 'PYTHON')

#Run on smoothed NED 19 max-min for line itself
ned19_smooth = FocalStatistics(NED19proj, NbrRectangle(3, 3, 'cells'), statistics_type='MEDIAN', ignore_nodata='DATA')
ned19_smooth.save(NED19smooth)
ZonalStatisticsAsTable(PSOSMras, 'Value', NED19smooth, out_table = rangetab19_smooth, statistics_type= 'RANGE', ignore_nodata='NODATA')
arcpy.AddField_management(rangetab19_smooth, 'osm_id', 'TEXT')
arcpy.CalculateField_management(rangetab19_smooth, 'osm_id', '!Value!', 'PYTHON')

#Run on smoothed NED 13 max-min for line itself
ned13_smooth = FocalStatistics(NED13proj, NbrRectangle(3, 3, 'cells'), statistics_type='MEDIAN', ignore_nodata='DATA')
ned13_smooth.save(NED13smooth)
ZonalStatisticsAsTable(PSOSMras, 'Value', NED13smooth, out_table = rangetab13_smooth, statistics_type= 'RANGE', ignore_nodata='NODATA')
arcpy.AddField_management(rangetab13_smooth, 'osm_id', 'TEXT')
arcpy.CalculateField_management(rangetab13_smooth, 'osm_id', '!Value!', 'PYTHON')

#Join all data to road vector
arcpy.MakeFeatureLayer_management(PSOSM_allproj, 'osmroads')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab19, 'osm_id')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab13, 'osm_id')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab19_smooth, 'osm_id')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab13_smooth, 'osm_id')
arcpy.CopyFeatures_management('osmroads', PSOSM_elv) #Often stays stuck in Python

arcpy.AlterField_management(PSOSM_elv, 'RANGE', 'RANGE19', 'RANGE19')
arcpy.AlterField_management(PSOSM_elv, 'RANGE_1', 'RANGE13', 'RANGE13')
arcpy.AlterField_management(PSOSM_elv, 'RANGE_12', 'RANGE19smooth', 'RANGE19smooth')
arcpy.AlterField_management(PSOSM_elv, 'RANGE_12_13', 'RANGE13smooth', 'RANGE13smooth')

########
arcpy.AddField_management(PSOSM_elv, 'gradient19', 'FLOAT')
arcpy.AddField_management(PSOSM_elv, 'gradient13', 'FLOAT')
arcpy.AddField_management(PSOSM_elv, 'gradient19_smooth', 'FLOAT')
arcpy.AddField_management(PSOSM_elv, 'gradient13_smooth', 'FLOAT')
arcpy.CalculateField_management(PSOSM_elv, 'gradient19', '!RANGE19!/!Shape_Length!', 'PYTHON')
arcpy.CalculateField_management(PSOSM_elv, 'gradient13', '!RANGE13!/!Shape_Length!', 'PYTHON')
arcpy.CalculateField_management(PSOSM_elv, 'gradient19_smooth', '!RANGE19smooth!/!Shape_Length!', 'PYTHON')
arcpy.CalculateField_management(PSOSM_elv, 'gradient13_smooth', '!RANGE13smooth!/!Shape_Length!', 'PYTHON')

#Compare to Seattle roads slope values. See whether 1/3 arc-sec is just less precise or also bias compared to 19
SpatialJoinLines_LargestOverlap(target_features=PSOSM_elv, join_features=roads,
                                out_fc=PSOSM_elv + 'Seattlejoin', outgdb=PSgdb, bufsize='10 meters', keep_all=True,
                                fields_select=['SLOPE_PCT', 'ARTDESCRIP'])

############# TO DO:
# Fill in and adjust values for those roads outside of NED 1/9 extent
#############

########################################################################################################################
# CREATE HEATMAPS
# of speedlimit(for Seattle), raw AADT (for Seattle), functional class-based  AADT (for Puget Sound) and
# BING (see Bing_format.py)
# Use a decay function to 'simulate' the pollution spread of various levels of traffic volume, speed, and congestion
########################################################################################################################
#SPEED LIMIT
arcpy.PolylineToRaster_conversion(roadstraffic_avg, value_field='SPEEDLIMIT', out_rasterdataset='Seattle_spdlm', priority_field='SPEEDLIMIT',cellsize=res)
heat_spdlm = FocalStatistics(path.join(gdb,'Seattle_spdlm'), neighborhood=NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernel100.txt'),
                             statistics_type='SUM', ignore_nodata='DATA') #It seems that full paths are needed to make this work
heat_spdlm.save('heat_spdlm')
heat_spdlm_int = Int(Raster('heat_spdlm')+0.5) #Constantly result in overall python crash?
heat_spdlm_int.save('heat_spdlm_int')
arcpy.CopyRaster_management('heat_spdlm_int', path.join(rootdir, 'results/heatspdlm_int'))

#Raw AADT (Seattle)
arcpy.PolylineToRaster_conversion(roadstraffic_avg, value_field='AADT_interp', out_rasterdataset='Seattle_AADT', priority_field='AADT_interp',cellsize=res)
customheatmap(kernel_dir=path.join(rootdir, 'results/bing'), in_raster=path.join(gdb, 'Seattle_AADT'),
              out_gdb = gdb, out_var='AADT', divnum=100, keyw='')

#OSM functional class-based AADT
arcpy.PolylineToRaster_conversion(PSOSM_allproj, value_field='fclassADT', out_rasterdataset=OSM_AADT,
                                  priority_field='fclassADT',cellsize=res)
customheatmap(kernel_dir=path.join(rootdir, 'results/bing'), in_raster=OSM_AADT,
              out_gdb = PSgdb, out_var='OSMAADT', divnum=100, keyw='')

#Bus transit
customheatmap(kernel_dir=path.join(rootdir, 'results/bing'), in_raster=PStransitras,
              out_gdb = path.join(rootdir, 'results/transit.gdb'), out_var='bustransit', divnum=100, keyw='log200')

############################################################################################################################
################################  TO DO:
#Road gradient
arcpy.PolylineToRaster_conversion(PSOSM_elv, value_field='gradient##', out_rasterdataset=OSM_gradient,
                                  priority_field='fclassADT', cellsize=res)
customheatmap(kernel_dir=path.join(rootdir, 'results/bing'), in_raster=OSM_gradient,
              out_gdb = PSgdb, out_var='OSMgradient', divnum=100, keyw='log[12]00')

########################################################################################################################
# GET TREES HEATMAP VALUES
# Select candidate species of trees from the City of Seattle's street-tree dataset and extract heatmap values at their
# location
########################################################################################################################
#Project
arcpy.Project_management(trees, 'trees_proj', UTM10)
#Get heat values for all trees
def Iter_ListRaster(workspaces, wildcard):
    outlist = []
    for ws in workspaces:
        arcpy.env.workspace = ws
        rlist = arcpy.ListRasters(wildcard)
        if rlist is not None:
            outlist.extend([path.join(ws, r) for r in rlist if path.join(ws, r) not in outlist])
    return outlist

heatlist = Iter_ListRaster([PSgdb, gdb, Binggdb], 'heat*')

ExtractMultiValuesToPoints('trees_proj', heatlist, bilinear_interpolate_values='BILINEAR')

#Get zoning
arcpy.Project_management(zoning, 'zoning_proj', UTM10)
arcpy.SpatialJoin_analysis('trees_proj', 'zoning_proj', 'trees_zoning', join_operation='JOIN_ONE_TO_ONE', match_option='WITHIN')
#Get census data
arcpy.Project_management(censustract, 'Tract_2010Census_proj', UTM10)
arcpy.SpatialJoin_analysis('trees_zoning', 'Tract_2010Census_proj', 'trees_zoning_census', join_operation='JOIN_ONE_TO_ONE', match_option='WITHIN')

########################################################################################################################
# EXPORT DATA
########################################################################################################################
#Export NLCD data to Puget Sound scale
arcpy.CopyRaster_management(NLCD_reclass, NLCD_reclass_PS)
#Export NLCD impervious data
imp = arcpy.sa.ExtractByMask(NLCD_imp, PSwatershed)
imp.save(NLCD_imp_PS)
#Compute focal stats
imp_mean = arcpy.sa.FocalStatistics(NLCD_imp_PS, neighborhood = NbrCircle(3, "CELL"), statistics_type= 'MEAN')
imp_mean.save(NLCD_imp_PS + '_mean.tif')