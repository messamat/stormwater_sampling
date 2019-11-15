#Author: Mathis Messager
#Purpose: create a heuristic model of traffic-related heavy metal pollution to select maple trees on which to sample moss
# in the city of Seattle, WA

import arcpy
from arcpy.sa import *
import os
import zipfile
import re
import numpy as np
import pandas as pd
import numpy as np
import time
import logging
import itertools

#Custom modules
from explode_overlapping import *
from SpatialJoinLines_LargestOverlap import *
from heatmap_custom import *
from GTFStoSHP import *

# create logger with GTFStoSHP
logger = logging.getLogger('GTFStoSHP')
logger.setLevel(logging.INFO)

#Set options
pd.set_option('display.max_columns', 10)
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'

seattleroads = os.path.join(rootdir, 'data/CitySeattle_20180601/Seattle_Streets/Seattle_Streets.shp')
kingroads = os.path.join(rootdir, 'data/King_201806/Metro_Transportation_Network_TNET_in_King_County_for_Car_Mode__trans_network_car_line/'
                               'Metro_Transportation_Network_TNET_in_King_County_for_Car_Mode__trans_network_car_line.shp')
traffic_seattle = os.path.join(rootdir, 'data/CitySeattle_20180601/2016_Traffic_Flow_Counts/2016_Traffic_Flow_Counts.shp')
traffic_wsdot = os.path.join(rootdir, 'data/WSDOT_TPTTraffic_20180508/2016_TrafficCounts/2016TrafficCounts.gdb/TrafficCounts2016')
trees = os.path.join(rootdir, 'data/CitySeattle_20180601/Trees/Trees.shp')
zoning = os.path.join(rootdir, 'data/CitySeattle_20180626/City_of_Seattle_Zoning/WGS84/City_of_Seattle_Zoning.shp')
censustract = os.path.join(rootdir, 'data/TIGER2017/Profile-County_Tract/Profile-County_Tract.gdb/Tract_2010Census_DP1')
heat_bing = os.path.join(rootdir, 'results/bing/bingmean1806_Seattle_heat.tif')
NLCD_reclass = os.path.join(rootdir, 'results/LU.gdb/NLCD_reclass_final') #Based on 2011 data
NLCD_imp = os.path.join(rootdir, 'data/NLCD_2016_Impervious_L48_20190405.img') #Based on 2016 dara
PSwatershed = os.path.join(rootdir, 'results/PSwtshd_dissolve.shp')
pscities = os.path.join(rootdir, 'results/PScitylimits.shp')
pscounties = os.path.join(rootdir, 'results/PScounties.shp')

template_ras = os.path.join(rootdir,'results/bing/bingeuc1902')
restemplate = arcpy.GetRasterProperties_management(template_ras, 'CELLSIZEX')
cs_ref = arcpy.Describe(NLCD_imp).SpatialReference

PSgdb=os.path.join(rootdir,'results/PSOSM.gdb')
OSMroads = os.path.join(rootdir, 'data/OSM_WA_20180601/gis_osm_roads_free_1.shp')
Kingroads = os.path.join(rootdir, 'data/King_201806\Metro_Transportation_Network_TNET_in_King_County_for_Car_Mode__trans_network_car_line\Metro_Transportation_Network_TNET_in_King_County_for_Car_Mode__trans_network_car_line.shp')
Pierceroads = os.path.join(rootdir, 'data/Pierce_20180611/Mobility_Data/Mobility_Data.shp')
traffic_wsdot = os.path.join(rootdir, 'data/WSDOT_TPTTraffic_20180508/2016_TrafficCounts/2016TrafficCounts.gdb/TrafficCounts2016')
kernel = NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernell00.txt') #UPDATE

STgtfs = os.path.join(rootdir, 'data\SoundTransit_201812\gtfs_puget_sound_consolidated.zip')
transitwiki_dir = os.path.join(rootdir, 'data/TransitWiki201812')

Binggdb = os.path.join(rootdir, 'results/bing/postprocess.gdb')
XRFsites = os.path.join(rootdir, 'data/field_data/sampling_sites_edit.shp')

#NED data paths to be updated when transfer to laptop
NED19proj = os.path.join(rootdir, 'results/ned19_psproj')
NED13proj = os.path.join(rootdir, 'results/ned13_psproj')
NED19int = os.path.join(rootdir, 'results/ned19_psint')
NED13int = os.path.join(rootdir, 'results/ned13_psint')
NED19smooth = os.path.join(rootdir, 'results/ned19_smooth')
NED13smooth = os.path.join(rootdir, 'results/ned13_smooth')
rangetab19 = os.path.join(PSgdb, 'OSM_elv19range')
rangetab13 = os.path.join(PSgdb, 'OSM_elv13range')
rangetab19_smooth = os.path.join(PSgdb, 'OSM_elv19range_smooth')
rangetab13_smooth = os.path.join(PSgdb, 'OSM_elv13range_smooth')
sroadsras = os.path.join(PSgdb, 'Seattle_roadras')
srangetab19 = os.path.join(PSgdb, 'Seattle_elv19range')
srangetab13 = os.path.join(PSgdb, 'Seattle_elv13range')

gdb = os.path.join(rootdir,'results/Seattle_sampling.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'Seattle_sampling.gdb')
arcpy.env.workspace = gdb

#New variables
OSMPierce = 'PSwtshd_OSMroads_Pierce.shp'
intersect= os.path.join(gdb, "buf_inters")
intersect_stats = os.path.join(gdb, "buf_inters_stats")
OSMPierce_datajoin = os.path.join(gdb, 'OSMPierce_datajoin')
OSMWSDOT_datajoin = os.path.join(gdb, 'OSM_WSDOT_joinstats')

roadstraffic = 'Seattle_roadstraffic'
roadstraffic_avg =roadstraffic+'_AADT'
NLCD_reclass_PS = os.path.join(rootdir, 'results/NLCD_reclass_final_PS.tif')
NLCD_imp_PS = os.path.join(rootdir, 'results/nlcd_imp_ps')
OSMSeattle = os.path.join(PSgdb, 'PSwtshd_OSMroads_Seattle')
OSMSeattle_datajoin = os.path.join(PSgdb, 'OSMSeattle_datajoin')
OSMKing_datajoin = os.path.join(PSgdb, 'OSMKing_datajoin')
OSMroads_proj = 'gis_osm_roads_free_1_aea'
PSOSM_all= os.path.join(rootdir, 'results/PSwtshd_OSMroads_all.shp')
OSM_AADT = os.path.join(PSgdb,'OSM_AADT')
OSM_SPD = os.path.join(PSgdb,'OSM_SPD')
PSOSMras = os.path.join(PSgdb, 'PSwtshd_OSMroads_all_projras')
PSOSM_elv = os.path.join(PSgdb, 'PSwtshd_OSMroads_all_proj_elv2')
OSM_gradient = os.path.join(PSgdb, 'OSM_gradient')
soundtransit = os.path.join(rootdir, 'results/transit.gdb/SoundTransit')

PStransit = os.path.join(rootdir, 'results/transit.gdb/PStransit')
PStransitbus = PStransit + '_busroutes'
PStransitbus_proj = PStransit + '_busroutes_proj'
PStransitbus_splitdiss = PStransitbus_proj + '_splitv_diss'
PStransitduplitab = os.path.join(rootdir, 'results/transit.gdb/explFindID')
PStransitras = os.path.join(rootdir, 'results/transit.gdb/PStransit_ras')

trees_aea = os.path.join(gdb, 'trees_aea')

########################################################################################################################
# PREPARE LAND USE DATA FOR PUGET SOUND
########################################################################################################################
#Export NLCD data to Puget Sound scale
arcpy.env.snapRaster = NLCD_imp
arcpy.Clip_management(in_raster=NLCD_reclass, rectangle=PSwatershed, out_raster=NLCD_reclass_PS)
#Export NLCD impervious data
arcpy.Clip_management(in_raster=NLCD_imp, rectangle=PSwatershed, out_raster=NLCD_imp_PS)
#Compute focal stats
imp_mean = arcpy.sa.FocalStatistics(NLCD_imp_PS, neighborhood = NbrCircle(3, "CELL"), statistics_type= 'MEAN')
imp_mean.save(NLCD_imp_PS + '_mean.tif')

########################################################################################################################
# PREPARE VARIABLES TO CREATE HEATMAPS (Puget Sound): FUNCTIONAL-CLASS BASED AADT AND SPEED LIMIT, SLOPE, AND TRANSIT ROUTES
# These outputs were used for initial moss sampling
########################################################################################################################
#-----------------------------------------------------------------------------------------------------------------------
# Prepare OSM data to create heatmap based on roads functional class for all Puget Sound OSM roads
#-----------------------------------------------------------------------------------------------------------------------
arcpy.env.workspace = gdb

#Select OSM roads, but this map with service roads, as enough through traffic to potentially have some impact
arcpy.Project_management(OSMroads, OSMroads_proj, out_coor_system= cs_ref)
arcpy.MakeFeatureLayer_management(OSMroads_proj, 'OSMroads_lyr')
np.unique([row[0] for row in arcpy.da.SearchCursor('OSMroads_lyr', ['fclass'])])
sel = "{} IN ('motorway','motorway_link','living_street','primary','primary_link','residential','secondary','secondary_link'," \
      "'tertiary','tertiary_link','trunk','trunk_link','unclassified','unknown', 'service')".format('"fclass"')
arcpy.SelectLayerByAttribute_management('OSMroads_lyr', 'NEW_SELECTION', sel)
arcpy.Intersect_analysis(['OSMroads_lyr', 'PSwtshd_dissolve.shp'],out_feature_class=PSOSM_all)
arcpy.Delete_management('OSMroads_lyr')

# Join OSM and Pierce County + WSDOT traffic counts data to improve interpolation of speed limits
# and traffic volume within road fclasses
#Subselect OSM roads for Pierce County
arcpy.MakeFeatureLayer_management(pscounties, 'counties_lyr')
arcpy.SelectLayerByAttribute_management('counties_lyr', 'NEW_SELECTION', "COUNTYNS='01529159'")
arcpy.Clip_analysis(PSOSM_all, 'counties_lyr', OSMPierce)

SpatialJoinLines_LargestOverlap(target_features= OSMPierce, join_features=Pierceroads, out_fc = OSMPierce_datajoin,
                                outgdb=gdb, bufsize='10 meters', keep_all=True,
                                fields_select=['RoadNumber', 'RoadName', 'FFC', 'FFCDesc', 'ADTSource', 'ADT',
                                               'ADTYear', 'SpeedLimit'])
#Join OSM with WSDOT traffic counts
arcpy.SpatialJoin_analysis(traffic_wsdot, PSOSM_all, os.path.join(gdb, 'OSM_WSDOT_join'), 'JOIN_ONE_TO_ONE', 'KEEP_COMMON',
                           match_option='CLOSEST_GEODESIC', search_radius='20 meters', distance_field_name='joindist')
arcpy.Statistics_analysis(os.path.join(gdb, 'OSM_WSDOT_join'),OSMWSDOT_datajoin,
                          statistics_fields= [['AADT', 'MEAN'],['DirectionOfTravel', 'FIRST'],['fclass', 'FIRST']],
                          case_field = 'osm_id')


#-----------------------------------------------------------------------------------------------------------------------
#ASSIGN FUNCTIONAL CLASS-BASED MEDIANS OF AADT AND SPEED LIMIT TO OSM FOR PUGET SOUND
#-----------------------------------------------------------------------------------------------------------------------
#Subselect OSM roads for Seattle
arcpy.MakeFeatureLayer_management(pscities, 'cities_lyr')
arcpy.SelectLayerByAttribute_management('cities_lyr', selection_type='NEW_SELECTION',
                                        where_clause="CityFIPSLo ='5363000'")
arcpy.Clip_analysis(PSOSM_all, 'cities_lyr', OSMSeattle)

#Join Seattle AADT data to OSM Seattle
SpatialJoinLines_LargestOverlap(target_features=OSMSeattle, join_features=os.path.join(gdb, roadstraffic_avg),
                                out_fc=OSMSeattle_datajoin, outgdb=PSgdb, bufsize='10 meters', keep_all=True,
                                fields_select=['AADT_avg', 'SPEEDLIMIT', 'ARTDESCRIP'])

#Join King County speed limit data to OSM Seattle
SpatialJoinLines_LargestOverlap(target_features=PSOSM_all, join_features=kingroads,
                                out_fc=OSMKing_datajoin, outgdb=PSgdb, bufsize='10 meters', keep_all=True,
                                fields_select=['SPEED_LIM'])
Klyr = arcpy.MakeFeatureLayer_management(OSMKing_datajoin, 'NOT SPEEDLIMIT IS NULL')
arcpy.CopyFeatures_management(Klyr, OSMKing_datajoin + '_sel')

#Create a compound table for all data on traffic counts
#(only keep OSM street segments whose buffer overlapped at least 80% with original data & whose length > 10 m
# & post-2000 Pierce County data)
fclass_ADTSPD = pd.DataFrame(
    [list(row[:4]) + [None,'Seattle'] + list(row[4:])
     for row in arcpy.da.SearchCursor(OSMSeattle_datajoin, ['osm_id', 'fclass', 'AADT_avg', 'SPEEDLIMIT', 'intersper', 'LENGTH_GEO']) \
     if row[2] is not None and row[4] > 0.8 and row[5] > 10] +
    \
    [list(row[:3])+ [row[4], int(row[3]), 'Pierce'] + list(row[5:])
     for row in arcpy.da.SearchCursor(OSMPierce_datajoin, ['osm_id', 'fclass', 'ADT', 'ADTYear', 'SpeedLimit', 'intersper', 'LENGTH_GEO']) \
     if row[2] is not None and row[3] != u' ' and row[5] > 0.8 and row[6] > 10] + \
 \
    [list(row[:2]) + [None, int(row[2]), None, 'King'] + list(row[3:])
     for row in
     arcpy.da.SearchCursor(OSMKing_datajoin, ['osm_id', 'fclass', 'SPEED_LIM', 'intersper', 'LENGTH_GEO']) \
     if row[3] > 0.8 and row[4] > 10] +
        \
    [list(row) + [None, None, 'WSDOT', None, None]
     for row in arcpy.da.SearchCursor(OSMWSDOT_datajoin, ['osm_id', 'FIRST_fclass', 'MEAN_AADT']) \
     if row[2] is not None],
    columns=['osm_id', 'fclass', 'ADT', 'SpeedLimit', 'year', 'agency', 'intersper', 'LENGTH_GEO'])

#------------- Process ADT data ---------------------
fclass_ADT = fclass_ADTSPD[((fclass_ADTSPD['year'] > 2000) | (fclass_ADTSPD['year'].isnull())) & \
                           (fclass_ADTSPD['ADT'].notnull())]
#Check segments that have multiple ADT counts but different values
duplis = fclass_ADT[(fclass_ADT.duplicated('osm_id', keep=False)) &
                    (~(fclass_ADT.duplicated(['osm_id','ADT'], keep=False)))].sort_values(by=['osm_id'])
#Keep most recent measurement for Pierce
fclass_ADT.loc[fclass_ADT['agency'] == 'Pierce'] = fclass_ADT.loc[fclass_ADT['agency'] == 'Pierce'].\
    sort_values('year', ascending=False).\
    drop_duplicates('osm_id').\
    sort_index()
#Average for each segment
fclass_ADTnodupli = fclass_ADT.groupby('osm_id').agg({'ADT' : ['mean'], 'fclass' : ['first', 'last']})
check = fclass_ADTnodupli[fclass_ADTnodupli['fclass','first'] != fclass_ADTnodupli['fclass','last']] #Check whether there are duplicate
#Compute median ADT for each OSM functional class
fclass_ADTnodupli.columns = fclass_ADTnodupli.columns.droplevel()
fclass_ADTmedian = fclass_ADTnodupli.groupby('first').median()

#------------- Process speed limit data ---------------------
fclass_ADTSPD['SpeedLimit'] = pd.to_numeric(fclass_ADTSPD['SpeedLimit'], errors='coerce')
fclass_SPD = fclass_ADTSPD[(fclass_ADTSPD['SpeedLimit'].notnull()) & (fclass_ADTSPD['SpeedLimit']>0)].\
    drop(['ADT', 'year'], axis=1)

# Check segments that have multiple speed limit counts but different values
duplis = fclass_SPD[(fclass_SPD.duplicated('osm_id', keep=False)) &
                    (~(fclass_SPD.duplicated(['osm_id','SpeedLimit'], keep=False)))].sort_values(by=['osm_id'])
# Average for each segment
fclass_SPDnodupli = fclass_SPD.groupby('osm_id').agg({'SpeedLimit' : ['mean'], 'fclass' : ['first', 'last']})
check = fclass_SPDnodupli[fclass_SPDnodupli['fclass','first'] != fclass_SPDnodupli['fclass','last']]
#Compute mean and median SPD for each OSM functional class
fclass_SPDnodupli.columns = fclass_SPDnodupli.columns.droplevel()
fclass_SPDmedian = fclass_SPDnodupli.groupby('first').median()
fclass_SPDmean = fclass_SPDnodupli.groupby('first').mean()

#Join SPD and AADT
fclass_SPDADTmedian = fclass_SPDmedian.join(fclass_ADTmedian, lsuffix='_SPD', rsuffix='_ADT')
fclass_SPDADTmedian.to_csv(os.path.join(rootdir, 'results/OSM_SPDADTmedian.csv'))


#Convert OSM functional categories to numbers
#Service road ADT doesn't really make sense, as includes mostly alleys.
#Unclassified seems to be misused in OSM and often appear to stand for unknwon
#Change both to residential level
arcpy.AddField_management(PSOSM_all, 'fclassADT', 'LONG')
arcpy.AddField_management(PSOSM_all, 'fclassSPD', 'LONG')
with arcpy.da.UpdateCursor(PSOSM_all, ['fclass','fclassADT', 'fclassSPD']) as cursor:
    for row in cursor:
        if row[0] in fclass_ADTmedian.index:
            if row[0] in ['service','unclassified']: #Value for service and unclassified seem overestimated
                row[1] = int(fclass_ADTmedian.loc['residential'])
            else:
                row[1] = int(fclass_ADTmedian.loc[row[0]])

            row[2] = int(fclass_SPDmedian.loc[row[0]])
        else:
            row[1]=0
            row[2]=0
        cursor.updateRow(row)

#-----------------------------------------------------------------------------------------------------------------------
# PREPARE TRANSIT DATA FOR PUGET SOUND TO CREATE HEATMAP BASED ON BUS ROUTES
#-----------------------------------------------------------------------------------------------------------------------
#Convert GTFS to routes with number of trips per week on each line
GTFStoSHPweeklynumber(gtfs_dir= STgtfs, out_gdb=os.path.dirname(soundtransit), out_fc = os.path.basename(soundtransit),
                      keep = False)

for gtfsdir in os.listdir(transitwiki_dir):
    indir = os.path.join(transitwiki_dir, gtfsdir)
    outname = re.sub('\W','_', os.path.splitext(gtfsdir)[0])
    if (os.path.isdir(indir) or zipfile.is_zipfile(indir)):
        if not arcpy.Exists(os.path.join(os.path.dirname(PStransit), outname + '_routes')):
            print(outname)
            # Create log to write out errors (https://docs.python.org/3/howto/logging.html#logging-basic-tutorial)
            errorlog = os.path.join(transitwiki_dir,
                                    datetime.now().strftime('errorlog_{}_%Y%m%d%H%M%S.log'.format(outname)))
            fh = logging.FileHandler(errorlog) # create file handler which logs even debug messages
            fh.setLevel(logging.WARNING) #Set handler level
            GTFStoSHPweeklynumber(gtfs_dir= indir, out_gdb=os.path.dirname(PStransit),
                                  out_fc = outname, keep=False)
            fh.close() #close handler
            if os.stat(errorlog).st_size == 0L: #Delete log if empty
                os.remove(errorlog)

#Merge all transit datasets
arcpy.Merge_management(arcpy.ListFeatureClasses('*_routes'), output = PStransit)

#Only keep buses with trips and whose schedule lasts more than 1 day
arcpy.MakeFeatureLayer_management(PStransit, 'PStransit_lyr',
                                  where_clause= '(route_type = 3) AND (MIN_service_len > 1) AND (SUM_adjustnum > 0)')
arcpy.CopyFeatures_management('PStransit_lyr', PStransitbus)
arcpy.Project_management(PStransitbus, PStransitbus_proj, cs_ref)

#Create raster of weekly number of buses at the same resolution as bing data
# Convert weekly number of buses to integer
arcpy.AddField_management(PStransitbus_proj, 'adjustnum_int', 'SHORT')
arcpy.CalculateField_management(PStransitbus_proj, 'adjustnum_int',
                                expression='int(10*!SUM_adjustnum!+0.5)', expression_type='PYTHON')

#Split lines at all intersections so that small identical overlapping segments can be dissolved
arcpy.SplitLine_management(PStransitbus_proj, PStransitbus_proj + '_split') #Split at intersection
arcpy.FindIdentical_management(in_dataset=PStransitbus_proj + '_split', out_dataset=PStransitduplitab, fields="Shape") #Find overlapping segments and make them part of a group (FEAT_SEQ)
arcpy.MakeFeatureLayer_management(PStransitbus_proj + '_split', "intlyr")
arcpy.AddJoin_management("intlyr", arcpy.Describe("intlyr").OIDfieldName, PStransitduplitab, "IN_FID", "KEEP_ALL")
arcpy.Dissolve_management("intlyr", PStransitbus_splitdiss, dissolve_field='explFindID.FEAT_SEQ',
                          statistics_fields=[[os.path.split(PStransitbus_proj)[1] + '_split.adjustnum_int', 'SUM']]) #Dissolve overlapping segments
arcpy.RepairGeometry_management(PStransitbus_splitdiss, delete_null = 'DELETE_NULL') #sometimes creates empty geom

#Get the length of a half pixel diagonal to create buffers for
#guaranteeing that segments potentially falling within the same pixel are rasterized separately
tolerance = (2.0**0.5)*float(restemplate.getOutput(0))/2
arcpy.env.workspace = os.path.dirname(soundtransit)
ExplodeOverlappingLines(PStransitbus_splitdiss, tolerance)

#For each set of non-overlapping lines, create its own raster
tilef = 'expl'
tilelist = list(set([row[0] for row in arcpy.da.SearchCursor(PStransitbus_splitdiss, [tilef])]))
outras_base = os.path.join(rootdir, 'results/transit.gdb/busnum_')
arcpy.env.snapRaster = template_ras
for tile in tilelist:
    outras = outras_base + str(tile)
    if not arcpy.Exists(outras):
        selexpr = '{0} = {1}'.format(tilef, tile)
        print(selexpr)
        arcpy.MakeFeatureLayer_management(PStransitbus_splitdiss, 'bus_lyr', where_clause= selexpr)
        arcpy.PolylineToRaster_conversion('bus_lyr', value_field='SUM_PStransit_busroutes_proj_split_adjustnum_int',
                                          out_rasterdataset=outras, cellsize=restemplate)

#Mosaic to new raster
arcpy.env.workspace = os.path.split(outras_base)[0]
transitras_tiles = arcpy.ListRasters('busnum_*')
arcpy.MosaicToNewRaster_management(transitras_tiles, arcpy.env.workspace, os.path.split(PStransitras)[1],
                                   pixel_type='32_BIT_UNSIGNED', number_of_bands= 1, mosaic_method = 'SUM')
for tile in transitras_tiles:
    print('Deleting {}...'.format(tile))
    arcpy.Delete_management(tile)
arcpy.ClearEnvironment('Workspace')

#----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA ON ROAD GRADIENTS (layers generated for initial sampling)
#-----------------------------------------------------------------------------------------------------------------------
#Run on NED19 max-min for line itself
arcpy.PolylineToRaster_conversion(PSOSM_all, 'osm_id', PSOSMras + '19', cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'fclassADT', cellsize = NED19proj)
ZonalStatisticsAsTable(PSOSMras + '19', 'osm_id', NED19proj, out_table = rangetab19, statistics_type= 'RANGE', ignore_nodata='NODATA')

#Run on NED13 max-min for line itself
arcpy.PolylineToRaster_conversion(PSOSM_all, 'osm_id', PSOSMras + '13', cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'fclassADT', cellsize = NED13proj)
ZonalStatisticsAsTable(PSOSMras + '13', 'osm_id', NED13proj, out_table = rangetab13, statistics_type= 'RANGE', ignore_nodata='NODATA')

#Run on smoothed NED 19 max-min for line itself
int19 = Int(Raster(NED19proj) + 0.5)
int19.save(NED19int)
ned19_smooth = FocalStatistics(NED19int, NbrRectangle(3, 3, 'cells'), statistics_type='MEDIAN', ignore_nodata='DATA')
ned19_smooth.save(NED19smooth)
ZonalStatisticsAsTable(PSOSMras + '19', 'osm_id', NED19smooth, out_table = rangetab19_smooth, statistics_type= 'RANGE', ignore_nodata='NODATA')

#Run on smoothed NED 13 max-min for line itself
int13 = Int(Raster(NED13proj) + 0.5)
int13.save(NED13int)
ned13_smooth = FocalStatistics(NED13int, NbrRectangle(3, 3, 'cells'), statistics_type='MEDIAN', ignore_nodata='DATA')
ned13_smooth.save(NED13smooth)
ZonalStatisticsAsTable(PSOSMras + '13', 'osm_id', NED13smooth, out_table = rangetab13_smooth, statistics_type= 'RANGE', ignore_nodata='NODATA')

#Join all data to road vector
arcpy.MakeFeatureLayer_management(PSOSM_all, 'osmroads')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab19, 'osm_id')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab13, 'osm_id')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab19_smooth, 'osm_id')
arcpy.AddJoin_management('osmroads', 'osm_id', rangetab13_smooth, 'osm_id')
arcpy.CopyFeatures_management('osmroads', PSOSM_elv) #Often stays stuck in Python

arcpy.AlterField_management(PSOSM_elv, 'RANGE', 'RANGE19', 'RANGE19')
arcpy.AlterField_management(PSOSM_elv, 'RANGE_1', 'RANGE13', 'RANGE13')
arcpy.AlterField_management(PSOSM_elv, 'RANGE_12', 'RANGE19smooth', 'RANGE19smooth')
arcpy.AlterField_management(PSOSM_elv, 'RANGE_12_13', 'RANGE13smooth', 'RANGE13smooth')

arcpy.AddField_management(PSOSM_elv, 'gradient19', 'FLOAT')
arcpy.AddField_management(PSOSM_elv, 'gradient13', 'FLOAT')
arcpy.AddField_management(PSOSM_elv, 'gradient19_smooth', 'FLOAT')
arcpy.AddField_management(PSOSM_elv, 'gradient13_smooth', 'FLOAT')
arcpy.CalculateField_management(PSOSM_elv, 'gradient19', '!RANGE19!/!Shape_Length!', 'PYTHON')
arcpy.CalculateField_management(PSOSM_elv, 'gradient13', '!RANGE13!/!Shape_Length!', 'PYTHON')
arcpy.CalculateField_management(PSOSM_elv, 'gradient19_smooth', '!RANGE19smooth!/!Shape_Length!', 'PYTHON')
arcpy.CalculateField_management(PSOSM_elv, 'gradient13_smooth', '!RANGE13smooth!/!Shape_Length!', 'PYTHON')

#Compare to Seattle roads slope values. Apply same method to Seattle road dataset
roadproj = arcpy.Project_management(seattleroads, os.path.join(PSgdb, 'Seattle_roadproj'), cs_ref)
arcpy.env.snapRaster = NED19proj
arcpy.PolylineToRaster_conversion(roadproj, 'OBJECTID', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED19proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED19proj, out_table = srangetab19,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED19smooth, out_table = srangetab19 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.env.snapRaster = NED13proj
arcpy.PolylineToRaster_conversion(roadproj, 'OBJECTID', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED13proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED13proj, out_table = srangetab13,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED13smooth, out_table = srangetab13 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
arcpy.ClearEnvironment('snapRaster')

# Fill in and adjust values for those roads outside of NED 1/9 extent
arcpy.AddField_management(PSOSM_elv, 'gradient_composite', 'FLOAT')
with arcpy.da.UpdateCursor(PSOSM_elv, ['gradient_composite','gradient19_smooth', 'gradient13_smooth']) as cursor:
    for row in cursor:
        if row[1] is not None:
            row[0] = min(row[1], 0.5)
        elif row[2] is not None:
            row[0] = min(row[2], 0.5)
        else:
            pass
        cursor.updateRow(row)

########################################################################################################################
# CREATE HEATMAPS
# of speedlimit(for Seattle), raw AADT (for Seattle), functional class-based  AADT (for Puget Sound) and
# BING (see Bing_format.py)
# Use a decay function to 'simulate' the pollution spread of various levels of traffic volume, speed, and congestion
########################################################################################################################
#Seattle SPEED LIMIT
arcpy.env.snapRaster = template_ras
arcpy.PolylineToRaster_conversion(roadstraffic_avg, value_field='SPEEDLIMIT', out_rasterdataset='Seattle_spdlm',
                                  priority_field='SPEEDLIMIT',cellsize=restemplate)
heat_spdlm = FocalStatistics(os.path.join(gdb,'Seattle_spdlm'), neighborhood=NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernel100.txt'),
                             statistics_type='SUM', ignore_nodata='DATA') #It seems that full paths are needed to make this work
heat_spdlm.save('heat_spdlm')
heat_spdlm_int = Int(Raster('heat_spdlm')+0.5) #Constantly result in overall python crash?
heat_spdlm_int.save('heat_spdlm_int')
arcpy.CopyRaster_management('heat_spdlm_int', os.path.join(rootdir, 'results/heatspdlm_int'))

#Seattle AADT
arcpy.PolylineToRaster_conversion(roadstraffic_avg, value_field='AADT_interp', out_rasterdataset='Seattle_AADT',
                                  priority_field='AADT_interp',cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(gdb, 'Seattle_AADT'),
              out_gdb = gdb, out_var='AADT', divnum=100, keyw='')

#OSM functional class-based AADT
arcpy.PolylineToRaster_conversion(PSOSM_all, value_field='fclassADT', out_rasterdataset=OSM_AADT,
                                  priority_field='fclassADT',cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=OSM_AADT,
              out_gdb = PSgdb, out_var='OSMAADT', divnum=100, keyw='')

#OSM functional class-based SPEED LIMIT
arcpy.PolylineToRaster_conversion(PSOSM_all, value_field='fclassSPD', out_rasterdataset=OSM_SPD,
                                  priority_field='fclassSPD',cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=OSM_SPD,
              out_gdb = PSgdb, out_var='OSMSPD', divnum=100, keyw='log[2]00')

#Bus transit
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=PStransitras,
              out_gdb = os.path.join(rootdir, 'results/transit.gdb'), out_var='bustransit', divnum=100, keyw='log500')

#Road gradient
arcpy.PolylineToRaster_conversion(PSOSM_elv, value_field='gradient_composite', out_rasterdataset=OSM_gradient,
                                  priority_field='fclassADT', cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=OSM_gradient,
              out_gdb = PSgdb, out_var='OSMgradient', divnum=0.01, keyw='log[123]00')

########################################################################################################################
# GET TREES HEATMAP VALUES
########################################################################################################################
#Project
arcpy.Project_management(trees, trees_aea, cs_ref)

#Get heat values for all trees
def Iter_ListRaster(workspaces, wildcard):
    outlist = []
    for ws in workspaces:
        arcpy.env.workspace = ws
        rlist = arcpy.ListRasters(wildcard)
        if rlist is not None:
            outlist.extend([os.path.join(ws, r) for r in rlist if os.path.join(ws, r) not in outlist])
    arcpy.ClearEnvironment('Workspace')
    return outlist

heatlist = Iter_ListRaster([PSgdb, gdb, os.path.join(rootdir, 'results/transit.gdb'), Binggdb], 'heat*') + \
           [NLCD_imp_PS, NLCD_imp_PS + '_mean.tif', NLCD_reclass_PS]

ExtractMultiValuesToPoints(trees_aea, heatlist, bilinear_interpolate_values='NONE')
arcpy.AddGeometryAttributes_management(trees_aea, 'POINT_X_Y_Z_M', Coordinate_System= arcpy.SpatialReference(4326))

#Get zoning
arcpy.Project_management(zoning, 'zoning_proj', cs_ref)
arcpy.SpatialJoin_analysis('trees_aea', 'zoning_proj', 'trees_zoning', join_operation='JOIN_ONE_TO_ONE', match_option='WITHIN')
#Get census data
arcpy.Project_management(censustract, 'Tract_2010Census_proj', cs_ref)
arcpy.SpatialJoin_analysis('trees_zoning', 'Tract_2010Census_proj', 'trees_zoning_census', join_operation='JOIN_ONE_TO_ONE', match_option='WITHIN')