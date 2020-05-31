from arcpy.sa import *
#import geopandas as gpd #use pd because arcgis is a pain with geojson 0.4.0
import arcpy
import os
import numpy as np

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False
arcpy.env.cellSize = "MINOF"

#Set up paths
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
resdir = os.path.join(rootdir, 'results')
datdir = os.path.join(rootdir, 'data')

#Create new gdb for predictions
predgdb = os.path.join(rootdir,'results/PSpredictions.gdb')
if arcpy.Exists(predgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'PSpredictions.gdb')

#Define variables
template_ras = os.path.join(rootdir, 'results/bing/heatbing1902log300proj.tif')
PSdissolve = os.path.join(resdir, 'PSwtshd_dissolve.shp')
PSdissolveras = os.path.join(predgdb, 'PSwtshd_dissolve')
NLCD = os.path.join(rootdir, 'data', 'NLCD_2016_Land_Cover_L48_20190424.img')

usdotgdb = os.path.join(resdir, 'usdot/usdot.gdb')
hpms_ps = os.path.join(usdotgdb, 'hpmstiger_PS')
hpms_ps_uhigh = os.path.join(usdotgdb, 'hpmstiger_PS_urbanhighways')
hpms_ps_uhighbuf = os.path.join(usdotgdb, 'hpmstiger_PS_urbanhighways_buf500')
hpms_ps_uhighbufras = os.path.join(usdotgdb, 'hpmstiger_PS_urbanhighways_buf500ras')

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
buildingsbuf100ras = os.path.join(predgdb, 'buildingsbuf100ras')

EJall = os.path.join(datdir,'EJSCREEN','EJSCREEN_2019_USPR.gdb','EJSCREEN_Full')
EJWA = os.path.join(predgdb, 'EJWA')
EJWAproj = os.path.join(predgdb, 'EJWAproj')
EJWA_ps = os.path.join(predgdb, 'EJWA_PS')
EJWAidras = os.path.join(predgdb, 'EJWAidras')
EJWAidras_sub = os.path.join(predgdb, 'EJWAidras_sub')
EJWAidras_subtab = os.path.join(predgdb, 'EJWAidras_subtab')
EJWApixelpopras = os.path.join(predgdb, 'EJWApixelpop100')
EJWApixeldipopras  = os.path.join(predgdb, 'EJWApixeldipop100')

EJWApixelpopras_sub = os.path.join(predgdb, 'EJWApixelpop10000sub')
EJWApixeldipopras_sub = os.path.join(predgdb, 'EJWApixeldipop10000sub')
predznpop = os.path.join(predgdb, 'predznpop')
predcupop = os.path.join(predgdb, 'predcupop')
predznexcessdipop = os.path.join(predgdb, 'predznexcessdipop')
predcuexcessdipop = os.path.join(predgdb, 'predcuexcessdipop')

samplepoints = os.path.join(predgdb, 'randompoints_ps')

def exportrastab(rasdat, out_table):
    if isinstance(rasdat, str):
        rasdat = Raster(rasdat)
    print('Create raster attribute table for {}...'.format(rasdat))
    arcpy.BuildRasterAttributeTable_management(rasdat)
    arcpy.CopyRows_management(rasdat, out_table)

arcpy.env.extent = PSdissolveras
arcpy.env.snapRaster = template_ras

####--------------- Create buffer around urban highways to exclude from analysis ----------------
#Subselect urban highways, buffer, rasterize them
if not arcpy.Exists(hpms_ps_uhighbufras):
    arcpy.MakeFeatureLayer_management(hpms_ps, 'hpmshighways',
                                      where_clause='(f_system = 1 OR f_system = 2) AND '
                                                   '(NOT (urban_code = 99999 OR urban_code = 0))')
    arcpy.CopyFeatures_management('hpmshighways', hpms_ps_uhigh)
    arcpy.Buffer_analysis(hpms_ps_uhigh, hpms_ps_uhighbuf, '500 meters', method='GEODESIC')


    arcpy.PolygonToRaster_conversion(hpms_ps_uhighbuf, 'state_code', hpms_ps_uhighbufras,
                                     cellsize=template_ras)

####--------------- Export data to compute area vs total pollution ----------------
#Only keep land pixels and export attribute tables
#Export attribute table
arcpy.env.extent = PSdissolveras
arcpy.env.snapRaster = template_ras

if not arcpy.Exists(predzn36land):
    predzn_land = Int(100*arcpy.sa.Con((Raster(NLCD) != 11) & (Raster(PSdissolveras) == 1), Raster(predzn36index)))
    predzn_land.save(predzn36land)
if not arcpy.Exists(predcu19land):
    predcu_land = Int(100 *arcpy.sa.Con((Raster(NLCD) != 11) & (Raster(PSdissolveras) == 1), Raster(predcu19index)))
    predcu_land.save(predcu19land)
if not arcpy.Exists(predpi30land):
    predpi_land = Int(100 *arcpy.sa.Con((Raster(NLCD) != 11) & (Raster(PSdissolveras) == 1), Raster(predpi30)))
    predpi_land.save(predpi30land)

if not all(arcpy.Exists(os.path.join(predgdb, t)) for t in ['predzn_tab', 'predcu_tab','predpi_tab']):
    exportrastab(predzn36land, out_table = os.path.join(predgdb, 'predzn_tab'))
    exportrastab(predcu19land, out_table = os.path.join(predgdb, 'predcu_tab'))
    exportrastab(predpi30land, out_table = os.path.join(predgdb, 'predpi_tab'))

#Produce boolean vector showing 1 - land pixels with Cu or Zn beyond background levels and 0 - pixels at background level
if not arcpy.Exists(os.path.join(predgdb, 'predcuznbool_land')):
    predallbool_land = arcpy.sa.Con((Raster(predcu19index)>0) | (Raster(predzn36index)>0 & (Raster(NLCD) != 11)), 1, 0)
    predallbool_land.save(os.path.join(predgdb, 'predcuznbool_land'))

arcpy.ClearEnvironment('extent')

####--------------- Compute area vs total pollution around urban highways ----------------
if not all(arcpy.Exists(os.path.join(predgdb, t)) for t in
           ['predznhighway_tab', 'predcuhighway_tab','predpihighway_tab']):
    exportrastab(ExtractByMask(Raster(predzn36land), Raster(hpms_ps_uhighbufras)),
                 out_table=os.path.join(predgdb,  'predznhighway_tab'))
    exportrastab(ExtractByMask(predcu19land, hpms_ps_uhighbufras),
                 out_table=os.path.join(predgdb,  'predcuhighway_tab'))
    exportrastab(ExtractByMask(predpi30land, hpms_ps_uhighbufras),
                 out_table=os.path.join(predgdb,  'predpihighway_tab'))

####--------------- Intersect with EJScreen -----------------------------
#-- Compute WA excess demographic index (Demographic Index for Block Group –Demographic Index for Washington)
#Subset EJSCREEN for WA
if not arcpy.Exists(EJWAproj):
    arcpy.MakeFeatureLayer_management(EJall, 'EJlyr', where_clause="SUBSTRING(ID, 1,2) = '53'")
    arcpy.CopyFeatures_management('EJlyr', EJWA)

    #Compute demographic index for Washington
    DIWA = np.divide(sum([row[0] for row in arcpy.da.SearchCursor(EJWA, ['VULEO'])]),
                     float(sum([row[0] for row in arcpy.da.SearchCursor(EJWA, ['ACSTOTPOP'])])))

    #For each block group, compute excess demographic index (as population)
    if not 'DIWEIGHTWA' in [f.name for f in arcpy.ListFields(EJWA)]:
        arcpy.AddField_management(EJWA, 'DIWEIGHTWA', 'FLOAT')

    #Compute block group population*(vulnerable population in block group (count)/vulnerable population if block group had average washington demographic index)
    with arcpy.da.UpdateCursor(EJWA, ['DIWEIGHTWA', 'ACSTOTPOP', 'VULEOPCT']) as cursor:
        for row in cursor:
            row[0] = row[1]*row[2]/DIWA
            cursor.updateRow(row)

    arcpy.Project_management(EJWA, EJWAproj, predcu19land)

#Intersect blockgroups with Puget Sound watershed polygon
if not arcpy.Exists(EJWA_ps):
    arcpy.MakeFeatureLayer_management(EJWAproj, 'EJWAproj_lyr')
    arcpy.SelectLayerByLocation_management('EJWAproj_lyr', 'INTERSECT', PSdissolve, selection_type='NEW_SELECTION')
    arcpy.CopyFeatures_management('EJWAproj_lyr', EJWA_ps)

# Convert blockgroup IDs to raster
if not arcpy.Exists(EJWAidras):
    arcpy.PolygonToRaster_conversion(EJWA_ps, value_field='ID', out_rasterdataset=EJWAidras,
                                     cell_assignment='MAXIMUM_AREA', cellsize=predcu19land)

#-- Spread demographic excess index to populated areas
#Create 100 m buffer around buildings and rasterize
if not arcpy.Exists(buildingsproj):
    buildingsdf = gpd.read_file(buildings_json)
    buildingsdf.to_file(buildings)
    arcpy.Project_management(buildings, buildingsproj, predcu19land)

if not arcpy.Exists(buildingsbuf100):
    arcpy.MakeFeatureLayer_management(buildings, 'buildingslyr')
    arcpy.SelectLayerByLocation_management('buildingslyr', 'WITHIN', EJWA_ps, selection_type='NEW_SELECTION')
    arcpy.Buffer_analysis('buildingslyr', buildingsbuf100, '100 meters', dissolve_option='NONE')

if not arcpy.Exists(buildingsbuf100ras)
    arcpy.env.extent = EJWA_ps
    arcpy.PolygonToRaster_conversion(buildingsbuf100, value_field='BUFF_DIST', out_rasterdataset=buildingsbuf100ras,
                                     cell_assignment='CELL_CENTER', cellsize=predcu19land)

#Subset EJWAidras to only keep areas that are (within 100 m of building OR in developed area) AND not water:
if not arcpy.Exists(EJWAidras_sub):
    Con((((20 < Raster(NLCD)) & (Raster(NLCD) < 25)) | (~IsNull(Raster(buildingsbuf100ras)))) & (
                Raster(NLCD) != 11),
        Raster(EJWAidras)).save(EJWAidras_sub)

if not (arcpy.Exists(EJWApixelpopras) and arcpy.Exists(EJWApixeldipopras)):
    #Count the number of pixels within selected areas for each block group
    EJWAvalueID_dict = {str(row[0]): row[1] for row in arcpy.da.SearchCursor(EJWAidras, ['Value', 'ID'])}
    EJWAvaluecount_dict = {str(row[0]):row[1] for row in arcpy.da.SearchCursor(EJWAidras_sub, ['Value', 'Count'])}
    EJWAIDcount_dict = {v: EJWAvaluecount_dict[k] for k, v in EJWAvalueID_dict.iteritems() if k in EJWAvaluecount_dict}

    #Compute per-pixel population in each blockgroup (x100 to get rid of decimals)
    if not 'PIXELPOP10000' in [f.name for f in arcpy.ListFields(EJWA_ps)]:
        arcpy.AddField_management(EJWA_ps, 'PIXELPOP10000', 'SHORT')
    if not 'PIXELDIWEIGHTWA10000' in [f.name for f in arcpy.ListFields(EJWA_ps)]:
        arcpy.AddField_management(EJWA_ps, 'PIXELDIWEIGHTWA10000', 'SHORT')

    with arcpy.da.UpdateCursor(EJWA_ps, ['ID', 'ACSTOTPOP', 'PIXELPOP10000', 'DIWEIGHTWA', 'PIXELDIWEIGHTWA10000']) as cursor:
        for row in cursor: #for each block group
            if 'ACSTOTPOP'==0: #if total population == 0
                row[2]=0 #per-pixel pop == 0
                row[4]=0 #per-pixel demographic index excess == 0
            elif row[0] in EJWAIDcount_dict.keys():
                # population per pixel = 10000*population/number of populated pixels in BG
                row[2] = int(10000*float(row[1])/EJWAIDcount_dict[row[0]]+0.5)
                row[4]= int(10000*float(row[3])/EJWAIDcount_dict[row[0]]+0.5)
            cursor.updateRow(row)

    #Rasterize per-pixel total population and per-pixel excess demographic index (based on WA average)
    arcpy.env.extent = PSdissolveras
    arcpy.PolygonToRaster_conversion(EJWA_ps, value_field='PIXELPOP10000', out_rasterdataset=EJWApixelpopras,
                                     cell_assignment='MAXIMUM_AREA', cellsize = predcu19land)
    arcpy.PolygonToRaster_conversion(EJWA_ps, value_field='PIXELDIWEIGHTWA10000', out_rasterdataset=EJWApixeldipopras,
                                     cell_assignment='MAXIMUM_AREA', cellsize = predcu19land)

if not (arcpy.Exists(EJWApixelpopras_sub) and arcpy.Exists(EJWApixeldipopras_sub)):
    Con(Raster(EJWAidras_sub)>0, Raster(EJWApixelpopras)).save(EJWApixelpopras_sub)
    Con(Raster(EJWAidras_sub)>0, Raster(EJWApixeldipopras)).save(EJWApixeldipopras_sub)

#Multiply pollution index by population and export attribute table
if not arcpy.Exists(predznpop):
    (Int(Raster(predzn36land)* Raster(EJWApixelpopras_sub)/100.0)).save(predznpop)
if not arcpy.Exists(predcupop):
    (Int(Raster(predcu19land) * Raster(EJWApixelpopras_sub)/100.0)).save(predcupop)

#Multiply pollution by EJ screen by excess demographic index and export attribute table
if not arcpy.Exists(predznexcessdipop):
    (Int(Raster(predzn36land) * Raster(EJWApixeldipopras_sub)/100.0)).save(predznexcessdipop)
if not arcpy.Exists(predcuexcessdipop):
    (Int(Raster(predcu19land) * Raster(EJWApixeldipopras_sub)/100.0)).save(predcuexcessdipop)

############### check output and run from here on #################################
if not all(arcpy.Exists(os.path.join(predgdb, t)) for t in
           ['EJWApixelpop10000sub_tab', 'EJWApixeldipop10000sub_tab',
            'predznpop10000_tab', 'predcupop10000_tab',
            'predcudiratio10000_tab', 'predzndiratio10000_tab']):
    exportrastab(EJWApixelpopras_sub, out_table=os.path.join(predgdb, 'EJWApixelpop10000sub_tab'))
    exportrastab(EJWApixeldipopras_sub, out_table=os.path.join(predgdb, 'EJWApixeldipop10000sub_tab'))
    exportrastab(predznpop, out_table=os.path.join(predgdb, 'predznpop10000_tab'))
    exportrastab(predcupop, out_table=os.path.join(predgdb, 'predcupop10000_tab'))
    exportrastab(predznexcessdipop, out_table=os.path.join(predgdb, 'predcudiratio10000_tab'))
    exportrastab(predcuexcessdipop, out_table=os.path.join(predgdb, 'predzndiratio10000_tab'))


if not all(arcpy.Exists(os.path.join(predgdb, t)) for t in
           ['EJWApixelpop10000subhighway_tab', 'EJWApixeldipop10000subhighway_tab',
            'predznpop10000highway_tab', 'predcupop10000highway_tab',
            'predcudiratio10000highway_tab', 'predzndiratio10000highway_tab']):
    exportrastab(ExtractByMask(Raster(EJWApixelpopras_sub), Raster(hpms_ps_uhighbufras)),
                               out_table=os.path.join(predgdb, 'EJWApixelpop10000subhighway_tab'))
    exportrastab(ExtractByMask(Raster(EJWApixeldipopras_sub), Raster(hpms_ps_uhighbufras)),
                               out_table=os.path.join(predgdb, 'EJWApixeldipop10000subhighway_tab'))
    exportrastab(ExtractByMask(Raster(predznpop), Raster(hpms_ps_uhighbufras)),
                               out_table=os.path.join(predgdb, 'predznpop10000highway_tab'))
    exportrastab(ExtractByMask(Raster(predcupop), Raster(hpms_ps_uhighbufras)),
                               out_table=os.path.join(predgdb, 'predcupop10000highway_tab'))
    exportrastab(ExtractByMask(Raster(predznexcessdipop), Raster(hpms_ps_uhighbufras)),
                               out_table=os.path.join(predgdb, 'predcudiratio10000highway_tab'))
    exportrastab(ExtractByMask(Raster(predcuexcessdipop), Raster(hpms_ps_uhighbufras)),
                               out_table=os.path.join(predgdb, 'predzndiratio10000highway_tab'))

#Create random sample to extract values for EJWA, predz, predcu, and highway masks
if not arcpy.Exists(samplepoints):
    arcpy.CreateRandomPoints_management(predgdb, os.path.split(samplepoints)[1],
                                        constraining_extent=arcpy.Describe(EJWApixeldipopras_sub).extent,
                                        number_of_points_or_field=1000000,
                                        minimum_allowed_distance=5)
    ExtractMultiValuesToPoints(samplepoints, [EJWApixeldipopras_sub, predzn36land, predcu19land, hpms_ps_uhighbufras],
                               bilinear_interpolate_values='NONE')