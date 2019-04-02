'''
Author: Mathis Messager
Contact info: messamat@uw.edu
Creation date: March 2019

Purpose: test various methods to compute road segment gradient using Seattle reference values for comparison.
The basic approach is to extract the minimum and maximum elevation values that each segment overlaps in a Digital
Elevation Model (DE) and then compute slope as the ratio of this range to the segment's length.

Three different sets of approaches were tested resulting in 12 different possible values:
- 1/9 vs 1/3 arcsecond DEM
- Splitting road segments into 10 m, 25 m subsegments and averaging across subsegments or taking the min-max for whole
    segments.
- Using raw DEM or first smoothing DEM with a 3x3 moving average window
'''

import arcpy
from arcpy.sa import *
import os
from collections import defaultdict

arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = False
arcpy.CheckOutExtension("Spatial")

rootdir = "D:\Mathis\ICSL\stormwater"
PSgdb=os.path.join(rootdir,'results/PSOSM.gdb')
NED19proj = os.path.join(rootdir, 'results/ned19_psproj')
NED13proj = os.path.join(rootdir, 'results/ned13_psproj')
NED19smooth = os.path.join(rootdir, 'results/ned19_smooth')
NED13smooth = os.path.join(rootdir, 'results/ned13_smooth')

Seattleroads = os.path.join(rootdir, 'data/CitySeattle_20180601/Seattle_Streets/Seattle_Streets.shp')
UTM10 = arcpy.SpatialReference(26910)
sroadsras = os.path.join(PSgdb, 'Seattle_roadras')
srangetab19 = os.path.join(PSgdb, 'Seattle_elv19range')
srangetab13 = os.path.join(PSgdb, 'Seattle_elv13range')
Seattleroadsproj = os.path.join(PSgdb, 'Seattle_roadproj')

#########################################################################################################################
#ANALYSIS

#Compare to Seattle roads slope values. Apply same method to Seattle road dataset
roadproj = arcpy.Project_management(Seattleroads, Seattleroadsproj, UTM10)

#Densify roads
arcpy.CopyFeatures_management(Seattleroadsproj, Seattleroadsproj + 'dens10m')
arcpy.Densify_edit(Seattleroadsproj + 'dens10m', densification_method='DISTANCE', distance='10', max_deviation='1.5')
#Split at vertices
arcpy.SplitLine_management(Seattleroadsproj + 'dens10m',Seattleroadsproj + 'split10')

#Compute statistics
arcpy.PolylineToRaster_conversion(Seattleroadsproj + 'split10', 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED19proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED19proj, out_table = srangetab19,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED19smooth, out_table = srangetab19 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.PolylineToRaster_conversion(Seattleroadsproj + 'split10', 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED13proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED13proj, out_table = srangetab13,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED13smooth, out_table = srangetab13 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

#Get all range values for zonal statistics tables
tablist =  [srangetab19, srangetab13, srangetab19 + '_smooth', srangetab13 + '_smooth']
elvdic = defaultdict(lambda: [0] * len(tablist))
for i in range(0, len(tablist)):
    print(tablist[i])
    for row in arcpy.da.SearchCursor(tablist[i], ['Value', 'RANGE']):
        elvdic[row[0]][i] = row[1]

#Compute slope
arcpy.AddField_management(Seattleroadsproj + 'split10', 'gradient19', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split10', 'gradient13', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split10', 'gradient19_smooth', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split10', 'gradient13_smooth', 'FLOAT')

with arcpy.da.UpdateCursor(Seattleroadsproj + 'split10', ['OBJECTID_1', 'Shape_Length',
                                                        'gradient19', 'gradient13',
                                                        'gradient19_smooth' ,'gradient13_smooth']) as cursor:
    for row in cursor:
        row[2] = elvdic[row[0]][0]/row[1] # gradient19 = RANGE_19/Shape_Length
        row[3] = elvdic[row[0]][1]/row[1] # gradient13 = RANGE_13/Shape_Length
        row[4] = elvdic[row[0]][2]/row[1] # gradient19_smooth = RANGE_19smooth/Shape_Length
        row[5] = elvdic[row[0]][3]/row[1] # gradient13_smooth = RANGE_13smooth/Shape_Length
        cursor.updateRow(row)

#Compute average slope for full roads
slopedic = defaultdict(lambda: [0]*5)
with arcpy.da.SearchCursor(Seattleroadsproj + 'split10', ['OBJECTID', 'Shape_Length',
                                              'gradient19', 'gradient13',
                                              'gradient19_smooth' ,'gradient13_smooth']) as cursor:
    for row in cursor:
        slopedic[row[0]][0] += row[1] #sum(Shape_Length)
        slopedic[row[0]][1] += row[1]*row[2] #sum(gradient19*Shape_Length)
        slopedic[row[0]][2] += row[1]*row[3]
        slopedic[row[0]][3] += row[1]*row[4]
        slopedic[row[0]][4] += row[1]*row[5]

arcpy.AddField_management(Seattleroadsproj, 'gradient19m_split10', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj, 'gradient13m_split10', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj, 'gradient19_smoothm_split10', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj, 'gradient13_smoothm_split10', 'FLOAT')
with arcpy.da.UpdateCursor(Seattleroadsproj, ['OBJECTID', 'gradient19m_split10', 'gradient13m_split10',
                                              'gradient19_smoothm_split10' ,'gradient13_smoothm_split10']) as cursor:
    for row in cursor:
        row[1] = slopedic[row[0]][1]/ slopedic[row[0]][0]  #gradient19m = sum(gradient19*Shape_Length)/sum(Shape_Length)
        row[2] = slopedic[row[0]][2]/ slopedic[row[0]][0]  #gradient13m = sum(gradient13*Shape_Length)/sum(Shape_Length)
        row[3] = slopedic[row[0]][3]/ slopedic[row[0]][0]
        row[4] = slopedic[row[0]][4]/ slopedic[row[0]][0]
        cursor.updateRow(row)

#########################################################################################################################
#TRY WITH 25 M SPLIT
#Densify roads
arcpy.CopyFeatures_management(Seattleroadsproj, Seattleroadsproj + 'dens25m')
arcpy.Densify_edit(Seattleroadsproj + 'dens25m', densification_method='DISTANCE', distance='25', max_deviation='1.5')
#Split at vertices
arcpy.SplitLine_management(Seattleroadsproj + 'dens25m',Seattleroadsproj + 'split25')

#Compute statistics
arcpy.PolylineToRaster_conversion(Seattleroadsproj + 'split25', 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED19proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED19proj, out_table = srangetab19,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED19smooth, out_table = srangetab19 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.PolylineToRaster_conversion(Seattleroadsproj + 'split25', 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED13proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED13proj, out_table = srangetab13,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED13smooth, out_table = srangetab13 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

#Get all range values for zonal statistics tables
tablist =  [srangetab19, srangetab13, srangetab19 + '_smooth', srangetab13 + '_smooth']
elvdic = defaultdict(lambda: [0] * len(tablist))
for i in range(0, len(tablist)):
    print(tablist[i])
    for row in arcpy.da.SearchCursor(tablist[i], ['Value', 'RANGE']):
        elvdic[row[0]][i] = row[1]

#Compute slope
arcpy.AddField_management(Seattleroadsproj + 'split25', 'gradient19', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split25', 'gradient13', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split25', 'gradient19_smooth', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split25', 'gradient13_smooth', 'FLOAT')

with arcpy.da.UpdateCursor(Seattleroadsproj + 'split25', ['OBJECTID_1', 'Shape_Length',
                                                        'gradient19', 'gradient13',
                                                        'gradient19_smooth' ,'gradient13_smooth']) as cursor:
    for row in cursor:
        row[2] = elvdic[row[0]][0]/row[1] # gradient19 = RANGE_19/Shape_Length
        row[3] = elvdic[row[0]][1]/row[1] # gradient13 = RANGE_13/Shape_Length
        row[4] = elvdic[row[0]][2]/row[1] # gradient19_smooth = RANGE_19smooth/Shape_Length
        row[5] = elvdic[row[0]][3]/row[1] # gradient13_smooth = RANGE_13smooth/Shape_Length
        cursor.updateRow(row)

#Compute average slope for full roads
slopedic = defaultdict(lambda: [0]*5)
with arcpy.da.SearchCursor(Seattleroadsproj + 'split25', ['OBJECTID', 'Shape_Length',
                                              'gradient19', 'gradient13',
                                              'gradient19_smooth' ,'gradient13_smooth']) as cursor:
    for row in cursor:
        slopedic[row[0]][0] += row[1] #sum(Shape_Length)
        slopedic[row[0]][1] += row[1]*row[2] #sum(gradient19*Shape_Length)
        slopedic[row[0]][2] += row[1]*row[3]
        slopedic[row[0]][3] += row[1]*row[4]
        slopedic[row[0]][4] += row[1]*row[5]

arcpy.AddField_management(Seattleroadsproj, 'gradient19m_split25', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj, 'gradient13m_split25', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj, 'gradient19_smoothm_split25', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj, 'gradient13_smoothm_split25', 'FLOAT')
with arcpy.da.UpdateCursor(Seattleroadsproj, ['OBJECTID', 'gradient19m_split25', 'gradient13m_split25',
                                              'gradient19_smoothm_split25' ,'gradient13_smoothm_split25']) as cursor:
    for row in cursor:
        row[1] = slopedic[row[0]][1]/ slopedic[row[0]][0]  #gradient19m = sum(gradient19*Shape_Length)/sum(Shape_Length)
        row[2] = slopedic[row[0]][2]/ slopedic[row[0]][0]  #gradient13m = sum(gradient13*Shape_Length)/sum(Shape_Length)
        row[3] = slopedic[row[0]][3]/ slopedic[row[0]][0]
        row[4] = slopedic[row[0]][4]/ slopedic[row[0]][0]
        cursor.updateRow(row)

#########################################################################################################################
#TRY WITH NO SPLIT
[f.name for f in arcpy.ListFields(Seattleroadsproj)]
#Compute statistics
arcpy.PolylineToRaster_conversion(Seattleroadsproj, 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED19proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED19proj, out_table = srangetab19,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED19smooth, out_table = srangetab19 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.PolylineToRaster_conversion(Seattleroadsproj , 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED13proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED13proj, out_table = srangetab13,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED13smooth, out_table = srangetab13 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

#Get all range values for zonal statistics tables
tablist =  [srangetab19, srangetab13, srangetab19 + '_smooth', srangetab13 + '_smooth']
elvdic = defaultdict(lambda: [0] * len(tablist))
for i in range(0, len(tablist)):
    print(tablist[i])
    for row in arcpy.da.SearchCursor(tablist[i], ['Value', 'RANGE']):
        elvdic[row[0]][i] = row[1]

#Compute slope
arcpy.AddField_management(Seattleroadsproj , 'gradient19', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj , 'gradient13', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj , 'gradient19_smooth', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj , 'gradient13_smooth', 'FLOAT')

with arcpy.da.UpdateCursor(Seattleroadsproj, ['OBJECTID_1', 'Shape_Length',
                                                        'gradient19', 'gradient13',
                                                        'gradient19_smooth' ,'gradient13_smooth']) as cursor:
    for row in cursor:
        row[2] = elvdic[row[0]][0]/row[1] # gradient19 = RANGE_19/Shape_Length
        row[3] = elvdic[row[0]][1]/row[1] # gradient13 = RANGE_13/Shape_Length
        row[4] = elvdic[row[0]][2]/row[1] # gradient19_smooth = RANGE_19smooth/Shape_Length
        row[5] = elvdic[row[0]][3]/row[1] # gradient13_smooth = RANGE_13smooth/Shape_Length
        cursor.updateRow(row)