from arcpy.sa import *
import geopandas as gpd #use pd because arcgis is a pain with geojson
import arcpy
import os
import numpy as np

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False
arcpy.env.cellSize = "MINOF"

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'
resdir = os.path.join(rootdir, 'results')
datdir = os.path.join(rootdir, 'data')

#Create new gdb for predictions
predgdb = os.path.join(rootdir,'results/PSpredictions.gdb')
if arcpy.Exists(predgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'PSpredictions.gdb')

#Define variables
PSdissolve = os.path.join(resdir, 'PSwtshd_dissolve.shp')
PSdissolveras = os.path.join(predgdb, 'PSwtshd_dissolve')
nlcd_ps = os.path.join(rootdir, 'results/nlcd_ps.tif')

predzn36index = os.path.join(predgdb, 'predzn36index')
predcu19index = os.path.join(predgdb, 'predcu19index')
predpi30 = os.path.join(predgdb, 'predpi30')
predzn36land = os.path.join(predgdb, 'predzn36indexland')
predcu19land = os.path.join(predgdb, 'predcu19indexland')
predpi30land = os.path.join(predgdb, 'predpi30land')

buildings_json = os.path.join(datdir, 'MicrosoftFootprints_20181012\Washington', 'Washington.geojson')
buildings = os.path.join(datdir, 'MicrosoftFootprints_20181012', 'Washington_buildings.shp')
buildingsproj = os.path.join(predgdb, 'buildingsproj')
buildingsbuf100 = os.path.join(predgdb, 'buildingsbuf100')
buildingsbuf100ras = os.path.join(predgdb, 'buildingsbuf100')


EJall = os.path.join(datdir,'EJSCREEN','EJSCREEN_2019_USPR.gdb','EJSCREEN_Full')
EJWA = os.path.join(predgdb, 'EJWA')
EJWAproj = os.path.join(predgdb, 'EJWAproj')
EJWAunion = os.path.join(predgdb, 'EJWA_buildingsunion')
EJWAidras = os.path.join(predgdb, 'EJWAidras')
EJWAidras_sub = os.path.join(predgdb, 'EJWAidras_sub')
EJWAidras_subtab = os.path.join(predgdb, 'EJWAidras_subtab')
EJWApixelpopras = os.path.join(predgdb, 'EJWApixelpop100')
EJWApixeldipopras  = os.path.join(predgdb, 'EJWApixeldipop100')

EJWApixelpopras_sub = os.path.join(predgdb, 'EJWApixelpop100sub')
EJWApixeldipopras_sub  = os.path.join(predgdb, 'EJWApixeldipop100sub')

####--------------- Export data to compute area vs total pollution ----------------
#Only keep land pixels and export attribute tables
#Export attribute table
predzn_land = arcpy.sa.Con((Raster(nlcd_ps) != 11) & (Raster(PSdissolveras) == 1), Raster(predzn36index))
predzn_land.save(predzn36land)
predcu_land = arcpy.sa.Con((Raster(nlcd_ps) != 11) & (Raster(PSdissolveras) == 1), Raster(predcu19index))
predcu_land.save(predcu19land)
predpi_land = arcpy.sa.Con((Raster(nlcd_ps) != 11) & (Raster(PSdissolveras) == 1), Raster(predpi30))
predpi_land.save(predpi30land)

arcpy.CopyRows_management(predzn36land, out_table = os.path.join(predgdb, 'predzn36_tab'))
arcpy.CopyRows_management(predcu19land, out_table = os.path.join(predgdb, 'predcu_tab'))
arcpy.CopyRows_management(predpi30land, out_table = os.path.join(predgdb, 'predpi_tab'))

#Produce boolean vector showing 1 - land pixels with Cu or Zn beyond background levels and 0 - pixels at background level
predallbool_land = arcpy.sa.Con((Raster(predcu19index)>0) | (Raster(predzn36index)>0 & (Raster(nlcd_ps) != 11)), 1, 0)
predallbool_land.save(os.path.join(predgdb, 'predcuznbool_land'))

####--------------- Intersect with EJScreen -----------------------------
#-- Compute WA excess demographic index (Demographic Index for Block Group –Demographic Index for Washington)
#Subset EJSCREEN for WA
arcpy.MakeFeatureLayer_management(EJall, 'EJlyr', where_clause="SUBSTRING(ID, 1,2) = '53'")
arcpy.CopyFeatures_management('EJlyr', EJWA)

#Compute demographic index for Washington
DIWA = np.divide(sum([row[0] for row in arcpy.da.SearchCursor(EJWA, ['VULEO'])]),
                 float(sum([row[0] for row in arcpy.da.SearchCursor(EJWA, ['ACSTOTPOP'])])))

#For each block group, compute excess demographic index (as population)
if not 'DISPEOWA' in [f.name for f in arcpy.ListFields(EJWA)]:
    arcpy.AddField_management(EJWA, 'DISPEOWA', 'FLOAT')
with arcpy.da.UpdateCursor(EJWA, ['DISPEOWA', 'VULEO', 'ACSTOTPOP']) as cursor:
    for row in cursor:
        row[0] = row[1] - (row[2]*DIWA)
        cursor.updateRow(row)

arcpy.Project_management(EJWA, EJWAproj, predcu19land)

#-- Spread demographic excess index to populated areas
#Create 100 m buffer around buildings
buildingsdf = gpd.read_file(buildings_json)
buildingsdf.to_file(buildings)
arcpy.Project_management(buildings, buildingsproj, predcu19land)

arcpy.MakeFeatureLayer_management(buildings, 'buildingslyr')
arcpy.SelectLayerByLocation_management('buildingslyr', 'WITHIN', PSdissolve, selection_type='NEW_SELECTION')
arcpy.Buffer_analysis('buildingslyr', buildingsbuf100, '100 meters', dissolve_option='NONE')

arcpy.env.mask = predcu19land
arcpy.env.snapRaster = predcu19land
arcpy.PolygonToRaster_conversion(EJWAproj, value_field='ID', out_rasterdataset=EJWAidras,
                                 cell_assignment='MAXIMUM_AREA', cellsize = predcu19land)
arcpy.PolygonToRaster_conversion(buildingsbuf100, value_field='BUFF_DIST', out_rasterdataset=buildingsbuf100ras,
                                 cell_assignment='CELL_CENTER', cellsize = predcu19land)

################### TO REWRITE ################################## GIVEN THAT BUILDING AND EJWAIDRAS ARE SEPARATE
EJWAidras_sub = os.path.join(predgdb, 'EJWAidras_sub')
#Subset EJWAidras to only keep areas that are (within 100 m of building OR in developed area) AND not water:
Con(((Mod(Raster(EJWAidras), 10)==1) | (InList(Raster(nlcd_ps),[21,22,23,24])>0)) &
    (Raster(nlcd_ps) != 11),
    EJWAidras).save(EJWAidras_sub)
EJWAid_dict = {row[0]:row[1] for row in arcpy.da.SearchCursor(EJWAidras_sub, ['IDCUSTOM', 'Count'])}
#EJWAidras_subtab = os.path.join(predgdb, 'EJWAidras_subtab')
#arcpy.CopyRows_management(EJWAidras_sub, out_table=EJWAidras_subtab)

#Compute per-pixel population in each blockgroup (x100 to get rid of decimals)
if not 'PIXELPOP100' in [f.name for f in arcpy.ListFields(EJWAproj)]:
    arcpy.AddField_management(EJWAproj, 'PIXELPOP100', 'SHORT')
if not 'PIXELDISPEOWA100' in [f.name for f in arcpy.ListFields(EJWAproj)]:
    arcpy.AddField_management(EJWAproj, 'PIXELDISPEOWA100', 'SHORT')

with arcpy.da.UpdateCursor(EJWAproj, ['CUSTOMID', 'ACSTOTPOP', 'PIXELPOP100', 'DISPEOWA', 'PIXELDISPEOWA100']) as cursor:
    for row in cursor:
        row[2] = row[1]/EJWAid_dict[row[0]]
        row[4]= row[3]/EJWAid_dict[row[0]]
        cursor.updateRow(row)

arcpy.PolygonToRaster_conversion(EJWAproj, value_field='PIXELPOP100', out_rasterdataset=EJWApixelpopras,
                                 cell_assignment='MAXIMUM_AREA', cellsize = predcu19land)
arcpy.PolygonToRaster_conversion(EJWAproj, value_field='PIXELDISPEOWA100', out_rasterdataset=EJWApixeldipopras,
                                 cell_assignment='MAXIMUM_AREA', cellsize = predcu19land)

EJWApixelpopras_sub = os.path.join(predgdb, 'EJWApixelpop100sub')
EJWApixeldipopras_sub  = os.path.join(predgdb, 'EJWApixeldipop100sub')
Con(Raster(EJWAidras_sub)>0, Raster(EJWApixelpopras)).save(EJWApixelpopras_sub)
Con(Raster(EJWAidras_sub)>0, Raster(EJWApixeldipopras)).save(EJWApixeldipopras_sub)

#Create raster of population by uniformally spreading through EJScreen block groups (same resolution and snapped to pollution rasters)

#Create raster of WA excess demographic index

######## Think about looking only at urban areas.
######## Outside of urban areas, could use all urban pixels + 100 m buffer around buildings to spread population

#Multiply pollution index by population and export attribute table

#Multiply pollution by EJ screen demographic index and export attribute table

#Multiple pollution index by population by EJ screen demographic index and export attribute table