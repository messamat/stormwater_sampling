#Author: Mathis Messager
#Purpose: create a heuristic model of traffic-related heavy metal pollution to select maple trees on which to sample moss
# in the city of Seattle, WA

import arcpy
from arcpy.sa import *
from os import *
import numpy as np
from collections import defaultdict

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

#Set up paths
rootdir = 'C:/Mathis/ICSL/stormwater/'

roads = path.join(rootdir, 'data/CitySeattle_20180601/Seattle_Streets/Seattle_Streets.shp')
traffic_seattle = path.join(rootdir, 'data/CitySeattle_20180601/2016_Traffic_Flow_Counts/2016_Traffic_Flow_Counts.shp')
traffic_wsdot = path.join(rootdir, 'data/WSDOT_TPTTraffic_20180508/2016_TrafficCounts/2016TrafficCounts.gdb/TrafficCounts2016')
trees = path.join(rootdir, 'data/CitySeattle_20180601/Trees/Trees.shp')
zoning = path.join(rootdir, 'data/CitySeattle_20180626/City_of_Seattle_Zoning/WGS84/City_of_Seattle_Zoning.shp')
censustract = path.join(rootdir, 'data/TIGER2017/Profile-County_Tract/Profile-County_Tract.gdb/Tract_2010Census_DP1')

kernel = NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernell00.txt') #UPDATE

gdb = path.join(rootdir,'results/Seattle_sampling.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir,'results'), 'Seattle_sampling.gdb')
arcpy.env.workspace = gdb

#New variables
roadstraffic = 'Seattle_roadstraffic'
roadstraffic_avg =roadstraffic+'_AADT'
UTM10 = arcpy.SpatialReference(26910)

########################################################################################################################
# GET TRAFFIC COUNT FOR EVERY ROAD SEGMENT IN SEATTLE
#   Average Annual Daily Traffic (AADT) is available for most arterial streets in the city of Seattle (here for 2016).
#   However, there is no common street identifier between the AADT dataset and the broader street dataset
#   of the city of Seattle. Therefore, spatially join the AADT dataset and the City of Seattle general street dataset.
#   Note that:
#   1. AADT segments are larger than City of Seattle streets, hence no need to split the street dataset.
#   2. Here average week day traffic is used as it has a fuller coverage of streets.

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
with arcpy.da.UpdateCursor(roadstraffic_avg, ['ARTDESCRIP','AADT_avg','AADT_interp']) as cursor:
    for row in cursor:
        if row[1] is not None:
            row[2] = row[1]
        if row[1] is None and row[2] is None and row[0] in ['','Not Designated']:
            row[2]=1000
        cursor.updateRow(row)

#Write out individual roads' attributes
arcpy.CopyRows_management(roadstraffic_avg, path.join(rootdir, 'results/Seattle_roads.dbf'))

########################################################################################################################

########################################################################################################################
# CREATE HEATMAPS OF SPEEDLIMIT AND AADT
# Use a logarithmic decay function to 'simulate' the pollution spread of various levels of traffic volume, speed,
# and congestion
########################################################################################################################
res = arcpy.GetRasterProperties_management(path.join(rootdir,'results/bing/180620_09_30_class_mlc.tif'), 'CELLSIZEX')
#SPEED LIMIT
arcpy.PolylineToRaster_conversion(roadstraffic_avg, value_field='SPEEDLIMIT', out_rasterdataset='Seattle_spdlm', priority_field='SPEEDLIMIT',cellsize=res)
heat_spdlm = FocalStatistics(path.join(gdb,'Seattle_spdlm'), neighborhood=NbrWeight('C:/Mathis/ICSL/stormwater/results/logkernel100.txt'),
                             statistics_type='SUM', ignore_nodata='DATA') #It seems that full paths are needed to make this work
heat_spdlm.save('heat_spdlm')
heat_spdlm_int = Int(Raster('heat_spdlm')+0.5) #Constantly result in overall python crash?
heat_spdlm_int.save('heat_spdlm_int')

#AADT
arcpy.PolylineToRaster_conversion(roadstraffic_avg, value_field='AADT_interp', out_rasterdataset='Seattle_AADT', priority_field='AADT_interp',cellsize=res)
heat_aadt = FocalStatistics(path.join(gdb,'Seattle_AADT'), neighborhood=kernel, statistics_type='SUM', ignore_nodata='DATA')
heat_aadt.save('heat_AADT')
heat_aadt_int = Int(Raster('heat_AADT'+0.5))
heat_aadt_int.save('heat_aadt_int')

#Get overall distribution of values in rasters
arcpy.BuildRasterAttributeTable_management('heat_AADT_int')
arcpy.BuildRasterAttributeTable_management('heat_spdlm_int')
arcpy.CopyRows_management('heat_AADT_int',  path.join(rootdir, 'results/heat_AADT.dbf'))
arcpy.CopyRows_management('heat_spdlm_int',  path.join(rootdir, 'results/heat_spdlm.dbf'))

########################################################################################################################
# GET TREES HEATMAP VALUES
# Select candidate species of trees from the City of Seattle's street-tree dataset and extract heatmap values at their
# location
########################################################################################################################
#Project
arcpy.Project_management(trees, 'trees_proj', UTM10)
#Get values
ExtractMultiValuesToPoints('trees_proj', ['heat_AADT','heat_spdlm'], bilinear_interpolate_values='BILINEAR')
#Get zoning
arcpy.Project_management(zoning, 'zoning_proj', UTM10)
arcpy.SpatialJoin_analysis('trees_proj', 'zoning_proj', 'trees_zoning', join_operation='JOIN_ONE_TO_ONE', match_option='WITHIN')
#Get census data
arcpy.Project_management(censustract, 'Tract_2010Census_proj', UTM10)
arcpy.SpatialJoin_analysis('trees_zoning', 'Tract_2010Census_proj', 'trees_zoning_census', join_operation='JOIN_ONE_TO_ONE', match_option='WITHIN')
#Export table
arcpy.CopyRows_management('trees_zoning_census', out_table=path.join(rootdir, 'results/trees_tab.dbf'))