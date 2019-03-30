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

#Compare to Seattle roads slope values. Apply same method to Seattle road dataset
roadproj = arcpy.Project_management(Seattleroads, Seattleroadsproj, UTM10)
#Densify roads
arcpy.Densify_edit(Seattleroadsproj, densification_method='DISTANCE', distance='10', max_deviation='1.5')
#Split at vertices
arcpy.SplitLine_management(Seattleroadsproj,Seattleroadsproj + 'split10')

#Compute statistics
arcpy.PolylineToRaster_conversion(Seattleroadsproj + 'split', 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'SURFACEWID', cellsize = NED19proj)
ZonalStatisticsAsTable(sroadsras, 'Value', NED19proj, out_table = srangetab19,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(sroadsras, 'Value', NED19smooth, out_table = srangetab19 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.PolylineToRaster_conversion(Seattleroadsproj + 'split', 'OBJECTID_1', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
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
arcpy.AddField_management(Seattleroadsproj + 'split', 'gradient19', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split', 'gradient13', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split', 'gradient19_smooth', 'FLOAT')
arcpy.AddField_management(Seattleroadsproj + 'split', 'gradient13_smooth', 'FLOAT')

with arcpy.da.UpdateCursor(Seattleroadsproj + 'split', ['OBJECTID_1', 'Shape_Length',
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
with arcpy.da.SearchCursor(Seattleroadsproj + 'split', ['OBJECTID', 'Shape_Length',
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
arcpy.Densify_edit(Seattleroadsproj, densification_method='DISTANCE', distance='25', max_deviation='1.5')
#Split at vertices
arcpy.SplitLine_management(Seattleroadsproj,Seattleroadsproj + 'split25')

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