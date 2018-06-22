import arcpy
from os import *
import numpy as np
from collections import defaultdict

rootdir = 'C:/Mathis/ICSL/stormwater/'
arcpy.env.overwriteOutput=True

roads = path.join(rootdir, 'data/CitySeattle_20180601/Seattle_Streets/Seattle_Streets.shp')
traffic_seattle = path.join(rootdir, 'data/CitySeattle_20180601/2016_Traffic_Flow_Counts/2016_Traffic_Flow_Counts.shp')
traffic_wsdot = path.join(rootdir, 'data/WSDOT_TPTTraffic_20180508/2016_TrafficCounts/2016TrafficCounts.gdb/TrafficCounts2016')

gdb = path.join(rootdir,'results/Seattle_sampling.gdb')
if arcpy.Exists(gdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(path.join(rootdir,'results'), 'Seattle_sampling.gdb')
arcpy.env.workspace = gdb

########################################################################################################################
# GET TRAFFIC COUNT FOR EVERY ROAD SEGMENT IN SEATTLE
########################################################################################################################
#Check for common fields between Seattle streets dataset and traffic_count dataset
roads_fields = [f.name for f in arcpy.ListFields(roads)]
traffic_fields = [f2.name for f2 in arcpy.ListFields(traffic_seattle)]
list(set(roads_fields).intersection(traffic_fields))

#Project datasets
print('Project datasets')
UTM10 = arcpy.SpatialReference(26910)
arcpy.Project_management(roads, 'roads_proj', out_coor_system=UTM10)
arcpy.Project_management(traffic_seattle, 'traffic_proj', out_coor_system=UTM10)

#Join Seattle traffic data to street data
print('Join Seattle traffic data to street data')
roadstraffic = 'Seattle_roadstraffic'
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

###################################################################################
#Join WSDOT traffic data to Seattle streets
print('Join WSDOT traffic data to Seattle streets')
arcpy.SpatialJoin_analysis(traffic_wsdot, roadstraffic, 'WSDOT_streets_join', 'JOIN_ONE_TO_ONE', 'KEEP_COMMON',
                           match_option='CLOSEST_GEODESIC', search_radius='50 meters', distance_field_name='joindist')
arcpy.Dissolve_management('WSDOT_streets_join', 'WSDOT_streets_join_diss', dissolve_field='CUSTOM_ID',
                          statistics_fields=[['AADT', 'MEAN'],['joindist', 'MEAN'],['Join_Count','SUM']])
arcpy.JoinField_management(roadstraffic, 'CUSTOM_ID', 'WSDOT_streets_join_diss', 'CUSTOM_ID', ['MEAN_joindist','MEAN_AADT', 'SUM_Join_Count'])
#Not a big deal if the match isn't perfect as the final output will be rasterized

#Average between WSDOT and Seattle estimates
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
roadstraffic_avg =roadstraffic+'_AADT'
arcpy.CopyFeatures_management('roadstraffic_lyr', roadstraffic_avg)
arcpy.Delete_management('roadstraffic_lyr')

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
#Fill in AADT interpolation for remaining arterial roads
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





########################################################################################################################
# CREATE HEATMAPS OF SPEEDLIMIT AND AADT
########################################################################################################################
res = arcpy.Describe(file.path(rootdir,'results/bing/180620_09_30_class_mlc.tif')).cellSize
arcpy.FeatureToRaster_conversion(roadstraffic_avg, field='AADT_interp', out_raster='Seattle_AADT', )
