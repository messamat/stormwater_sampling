import arcpy
from arcpy.sa import *
import pandas as pd
import geopandas as gpd
import os
import re
import requests
import sys
import urllib
import urllib2
import io
import glob
from BeautifulSoup import BeautifulSoup
import zipfile
import numpy as np
import traceback
import us
import itertools
import shutil
import contextlib
import ftplib
import urlparse
from collections import defaultdict

from SpatialJoinLines_LargestOverlap import *
from GTFStoSHP import *
from explode_overlapping import *
from heatmap_custom import *
from Download_gist import *

pd.options.display.max_columns = 20
pd.options.display.max_rows = 200
arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = True
arcpy.CheckOutExtension("Spatial")

rootdir = "D:\Mathis\ICSL\stormwater"
USDOTdir = os.path.join(rootdir, "data\USDOT_0319")
NTMdir = os.path.join(rootdir, "data\NTM_0319")
resdir = os.path.join(rootdir, 'results/usdot')
AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')

NED19proj = os.path.join(rootdir, 'results/ned19_psproj')
NED13proj = os.path.join(rootdir, 'results/ned13_psproj')

XRFsites = os.path.join(rootdir, 'data/field_data/sampling_sites_edit.shp')
XRFsiteshull = os.path.join(rootdir, 'results/XRFsites_convexhull.shp')
XRFsiteshull_aea = os.path.join(rootdir, 'results/XRFsites_convexhull_aea.shp')
bingdir = os.path.join(rootdir, 'results/bing/')
binggdb = os.path.join(bingdir, 'postprocess.gdb')
bingeuc = os.path.join(binggdb, 'bingeuc1902')

NLCD_reclass = os.path.join(rootdir, 'results/LU.gdb/NLCD_reclass_final')
NLCD_imp = os.path.join(rootdir, 'data/NLCD_2016_Impervious_L48_20190405.img') #Zipped
cs_ref = arcpy.Describe(NLCD_imp).SpatialReference

STgtfs = os.path.join(rootdir, 'data\SoundTransit_201812\gtfs_puget_sound_consolidated.zip')
transitwiki_dir = os.path.join(rootdir, 'data/TransitWiki201812')

PSdissolve = os.path.join(rootdir, 'results/PSwtshd_dissolve.shp')

template_ras = os.path.join(rootdir, 'results/bing/heatbing1902log300proj.tif')
restemplate = arcpy.GetRasterProperties_management(template_ras, 'CELLSIZEX')
pollutgdb = os.path.join(rootdir,'results/pollution_variables.gdb')
if arcpy.Exists(pollutgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'pollution_variables.gdb')

pollutgdbPS = os.path.join(rootdir, 'results/pollution_variables_PS.gdb')
if arcpy.Exists(pollutgdbPS):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir, 'results'), 'pollution_variables_PS.gdb')

#Get buffer around AQI sites
AQIdir = os.path.join(rootdir, 'results/airdata')
AQIsites = os.path.join(AQIdir, 'airsites.shp')
AQIsites_bufunion = os.path.join(AQIdir, 'airsites_550bufunion.shp')
AQIgdb = os.path.join(AQIdir, 'AQI.gdb')
AQIaadttab = os.path.join(AQIdir, 'airsites_aadtlog100')
AQIbingtab = os.path.join(AQIdir, 'airsites_binglog300')
AQIimptab = os.path.join(AQIdir, 'airsites_impmean')
roadAQI = os.path.join(AQIgdb, 'AQI_hpmstigerinters')

#Create gdb for analysis
if arcpy.Exists(AQIgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results/airdata'), 'AQI.gdb')

#Output variables
usdotgdb = os.path.join(resdir, 'usdot.gdb')
hpms = os.path.join(USDOTdir, '2016.gdb/HPMS2016')
hmps_spdl = os.path.join(usdotgdb, 'hpms_spdl')

NED19smooth = os.path.join(rootdir, 'results/ned19_smooth')
NED13smooth = os.path.join(rootdir, 'results/ned13_smooth')
rangetab19 = os.path.join(pollutgdb, 'hpms_elv19range')
rangetab13 = os.path.join(pollutgdb, 'hpms_elv13range')
rangetab19_smooth = os.path.join(pollutgdb, 'hpms_elv19range_smooth')
rangetab13_smooth = os.path.join(pollutgdb, 'hpms_elv13range_smooth')

tiger16dir = os.path.join(rootdir, 'data/TIGER2016')
tigerroads = os.path.join(usdotgdb, 'tigerroads_merge')
tigerroads_sub = os.path.join(usdotgdb, 'tigerroads_sub')
tigerroads_UA = os.path.join(usdotgdb, 'tigerroads_UA')
tigerroads_format = tigerroads_UA + 'states'
states = os.path.join(tiger16dir, 'tl_2016_us_state/tl_2016_us_state.shp')
hpmstiger = os.path.join(usdotgdb, 'hpmstiger')
hpmstigerproj = os.path.join(usdotgdb, 'hpmstiger_proj')

hpms_sub = os.path.join(usdotgdb, 'hpmstiger_subset')
hpms_PS = os.path.join(usdotgdb, 'hpmstiger_PS')
hpms_ras19 = os.path.join(usdotgdb, 'hpmstiger_subsetras19')
hpms_ras13 = os.path.join(usdotgdb, 'hpmstiger_subsetras13')

NLCD_reclass_PS = os.path.join(rootdir, 'results/NLCD_reclass_final_PS.tif')
NLCD_imp_PS = os.path.join(rootdir, 'results/nlcd_imp_ps')
NLCD_imp_AQI_mean = os.path.join(rootdir, 'results/nlcd_imp_AQI_mean.tif')

NTMroutes = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes')
NTMsel = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes_sel')
NTMproj = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes_selproj')
NTMsplitdiss = NTMproj + '_splitv_diss'
NTMras = os.path.join(rootdir, 'results/transit.gdb/NTMras')

transitgdb = os.path.join(rootdir, 'results/transit.gdb')
soundtransit = os.path.join(rootdir, 'results/transit.gdb/SoundTransit')
PStransit = os.path.join(rootdir, 'results/transit.gdb/PStransit')
PStransitbus = PStransit + '_busroutes'
PStransitbus_proj = PStransit + '_busroutes_proj'
PStransitbus_splitdiss = PStransitbus_proj + '_splitv_diss'
PStransitras = os.path.join(rootdir, 'results/transit.gdb/PStransit_ras')

if not arcpy.Exists(usdotgdb):
    arcpy.CreateFileGDB_management(resdir, out_name = 'usdot')

XRFsites_aea = os.path.join(pollutgdb, 'XRFsites_aea')
XRFsites_aeaattri = os.path.join(pollutgdb, 'XRFsites_projallattri')

def Iter_ListRaster(workspaces, wildcard):
    "Build list of all rasters in list of workspaces that correspond to a wildcard"
    outlist = []
    for ws in workspaces:
        arcpy.env.workspace = ws
        rlist = arcpy.ListRasters(wildcard)
        if rlist is not None:
            outlist.extend([os.path.join(ws, r) for r in rlist if os.path.join(ws, r) not in outlist])
    arcpy.ClearEnvironment('Workspace')
    return outlist

########################################################################################################################
# COMPUTE FUNCTIONAL-CLASS BASED AADT AVERAGE FOR EVERY STATE
########################################################################################################################
#Get VM-2 estimates of Vehicle Miles Traveled for each state and average it for small roads at state level (idea fromhttps://www.pnas.org/content/pnas/suppl/2015/04/01/1421723112.DCSupplemental/pnas.1421723112.sapp.pdf)
#Highway Statistics Book 2017
# Public road lenth 2017 - miles by functional system, Table HM-20: https://www.fhwa.dot.gov/policyinformation/statistics/2017/hm20.cfm
hm20url = "https://www.fhwa.dot.gov/policyinformation/statistics/2016/hm20.cfm"

hm20_page = urllib2.urlopen(hm20url)
hm20_soup = BeautifulSoup(hm20_page)
hm20_regex = re.compile(".*hm20.*[.]xls")
hm20urltab = os.path.split(hm20url)[0] + '/' + \
             '/'.join([i for i in hm20_soup.find('a', attrs={'href': hm20_regex}).get('href').split('/')
                       if i not in os.path.split(hm20url)[0].split('/')])
hm20tab = os.path.join(USDOTdir, os.path.split(hm20urltab)[1])
if not arcpy.Exists(hm20tab):
    dlfile(hm20urltab, hm20tab)
else:
    print('{} already exists...'.format(hm20tab))

# Functional system travel 2017 - Annual vehicle-miles, Table VM-2: https://www.fhwa.dot.gov/policyinformation/statistics/2017/vm2.cfm
vm2url = "https://www.fhwa.dot.gov/policyinformation/statistics/2016/vm2.cfm"
vm2_page = urllib2.urlopen(vm2url)
vm2_soup = BeautifulSoup(vm2_page)
vm2_regex = re.compile(".*vm2.*[.]xls")
vm2urltab = os.path.split(vm2url)[0] + '/' + \
             '/'.join([i for i in vm2_soup.find('a', attrs={'href': vm2_regex}).get('href').split('/')
                       if i not in os.path.split(vm2url)[0].split('/')])
vm2tab = os.path.join(USDOTdir, os.path.split(vm2urltab)[1])
if not arcpy.Exists(vm2tab):
    dlfile(vm2urltab, vm2tab)
else:
    print('{} already exists...'.format(vm2tab))

#Format hm20 table
hm20pd = pd.ExcelFile(hm20tab)
hm20sheet = hm20pd.parse(hm20pd.sheet_names[1])
#Format header row
hm20sheet_tab = hm20sheet.loc[8:,:] #Remove title
hm20sheet_tab.loc[8,:][1:9] = hm20sheet_tab.loc[8,:][1] #Extend rural across columns
hm20sheet_tab.loc[8,:][9:17] = hm20sheet_tab.loc[8,:][9] #Extend urban across colums
hm20sheet_format = hm20sheet_tab.loc[12:,:].reset_index(drop=True)
hm20sheet_format.columns = hm20sheet_tab.loc[8:11,:].apply(lambda col: col.str.cat(), axis=0)
#Drop totals
hm20sheet_format = hm20sheet_format.drop(axis=1, labels = [col for col in hm20sheet_format.columns if re.compile('.*TOTAL$').match(col)])
hm20sheet_format = hm20sheet_format.loc[hm20sheet_format.STATE.str.contains('.*(?<!\sTotal)$'),:]

#Format vm2 table
vm2pd = pd.ExcelFile(vm2tab)
vm2sheet = vm2pd.parse(vm2pd.sheet_names[1])
#Format header row
vm2sheet_tab = vm2sheet.loc[8:65,:] #Remove title
vm2sheet_tab.loc[8,:][1:9] = vm2sheet_tab.loc[8,:][1] #Extend rural across columns
vm2sheet_tab.loc[8,:][9:17] = vm2sheet_tab.loc[8,:][9] #Extend urban across colums
vm2sheet_format = vm2sheet_tab.loc[12:,:].reset_index(drop=True)
vm2sheet_format.columns = vm2sheet_tab.loc[8:11,:].apply(lambda col: col.str.cat(), axis=0)
#Drop totals
vm2sheet_format = vm2sheet_format.drop(axis=1, labels = [col for col in vm2sheet_format.columns if re.compile('.*TOTAL$').match(col)])
vm2sheet_format = vm2sheet_format.loc[vm2sheet_format.STATE.str.contains('.*(?<!\sTotal)$'),:]

#Compute AADT for each road type as Yearly Vehicle Miles Traveled/(Total length of roads*365)
np.where((vm2sheet_format == 0) != (hm20sheet_format == 0)) #No mismatch in cells with 0 values (0 miles traveled is always matched by 0 miles of road)
USAADT = vm2sheet_format.copy()
USAADT.iloc[:,1:][USAADT.iloc[:,1:] != 0] = \
    ((10**6)*vm2sheet_format[vm2sheet_format!=0].iloc[:,1:])\
    .div(hm20sheet_format[hm20sheet_format!=0].iloc[:,1:])\
    .div(365)
USAADT.set_index("STATE", inplace=True)

#Join table to states FIPS code
statesfips = pd.DataFrame.from_dict(us.states.mapping('name', 'fips'), orient='index')
statesfips.rename(columns={0: 'state_code'}, inplace=True)
USAADTfips = USAADT.join(statesfips)

#Format pd of AADT value for merging with spatial dataset
USAADTmelt = pd.melt(USAADTfips, id_vars='state_code', value_vars=[f for f in USAADTfips.columns if f != 'state_code'],
        var_name='fclass_original', value_name='AADT')
USAADTmelt['urban'] = USAADTmelt.fclass_original.str[0:5]
USAADTmelt['fclass_original'] = USAADTmelt.fclass_original.str[5:]
USAADTmelt.set_index('fclass_original', inplace = True)
#Add F_SYSTEM code (see HPMS field manual: https://www.fhwa.dot.gov/policyinformation/hpms/fieldmanual/hpms_field_manual_dec2016.pdf)
f_ref = {'INTERSTATE':1,
         'OTHERFREEWAYS  ANDEXPRESSWAYS':2,
         'OTHERPRINCIPALARTERIAL':3,
         'MINORARTERIAL':4,
         'MAJORCOLLECTOR':5,
         'MINORCOLLECTOR':6,
         'LOCAL':7}
USAADTformat = USAADTmelt.join(pd.DataFrame.from_dict(f_ref, 'index'))\
    .rename(columns={0:'F_SYSTEM'})
USAADTformat['AADT'] = USAADTformat['AADT'].astype(int)
USAADTformat['state_code'] = USAADTformat['state_code'].astype(int)
USAADTformat.to_csv(os.path.join(resdir, 'checkAADT.csv'))

#Download USHPMS data (https://www.bts.gov/geography/geospatial-portal/NTAD-direct-download)
if not arcpy.Exists(hpms):
    dlfile(url="http://www.bts.gov/sites/bts.dot.gov/files/ntad/HPMS2016.gdb.zip",
           outpath=USDOTdir)

#Fill in AADT value for all roads based on state-wide averages
#[[f.name, f.type] for f in arcpy.ListFields(hpms)]
if 'aadt_filled' not in [f.name for f in arcpy.ListFields(hpms)]:
    arcpy.AddField_management(hpms, field_name = 'aadt_filled', field_type='LONG')
else:
    print('aadt_filled field already exists...')
arcpy.CalculateField_management(hpms, 'aadt_filled', expression = '!aadt!', expression_type='PYTHON')
with arcpy.da.UpdateCursor(hpms, ['aadt', 'urban_code', 'aadt_filled', 'state_code', 'f_system', 'OBJECTID']) as cursor:
    for row in cursor:
        print(row[5])
        if row[0] == 0: #If does not already have an AADT
            if row[4] == 0: #If no functional classification, assign local AADT
                if row[1] in [0, 99999]: #If no urban_code (aka is a rural road or from small urban section)
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'RURAL') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == 7)]['AADT'][0]
                else:
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'URBAN') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == 7)]['AADT'][0]
            else: #Otherwise, assign corresponding f-system aadt
                if row[1] in [0, 99999]: #If no urban_code (aka is a rural road or from small urban section)
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'RURAL') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == row[4])]['AADT'][0]
                else:
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'URBAN') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == row[4])]['AADT'][0]
            cursor.updateRow(row)

########################################################################################################################
# FILL IN LOCAL ROADS FOR STATES WHICH DID NOT SUBMIT TO THE HPMS
########################################################################################################################
#Identify states that did not include local roads
fsystemtab = os.path.join(usdotgdb, 'hpms_fsystemstats')
if not arcpy.Exists(fsystemtab):
    arcpy.Statistics_analysis(hpms, statistics_fields=[['f_system','COUNT']], case_field=['state_code', 'f_system'])
fsystemstats = pd.DataFrame(arcpy.da.TableToNumPyArray(fsystemtab, field_names=('state_code','f_system','frequency')))
fsystemstats['proportion'] = fsystemstats['frequency']/\
                             fsystemstats['frequency'].groupby(fsystemstats['state_code']).transform('sum')
missingstatefips = list(fsystemstats[(fsystemstats['f_system']==7) & (fsystemstats['proportion']<0.05)]['state_code'])

#Download list of all FIPS codes in the US (see https://www.census.gov/geo/reference/codes/cou.html for metadata)
dlfile(url = 'https://www2.census.gov/geo/docs/reference/codes/files/national_county.txt', outpath = tiger16dir,
       fieldnames = ['STATE', 'STATEFP', 'COUNTYFP', 'COUNTYNAME', 'CLASSFP'])
fipslist = pd.read_csv(os.path.join(tiger16dir, 'national_county.txt'))

missingcountyfips = fipslist[fipslist['STATEFP'].isin(missingstatefips)]
missingcountyfips_list = list(missingcountyfips['STATEFP'].astype(str).str.zfill(2) +
                              missingcountyfips['COUNTYFP'].astype(str).str.zfill(3))

#Download shapefile data for all counties with missing local roads in the US
downloadroads(countyfipslist=missingcountyfips_list, year=2016, outdir=tiger16dir)

#Unzip shapefiles
for county_code in missingcountyfips_list:
    outbase = 'tl_2016_{}_roads'.format(county_code)
    if not arcpy.Exists(os.path.join(tiger16dir, '{}.shp'.format(outbase))):
        inzip = os.path.join(tiger16dir, '{}.zip'.format(outbase))
        if os.path.exists(inzip):
            unzip(inzip)
        else:
            print('{} has not yet been downloaded...'.format(outbase))
    else:
        print('{} already exists...'.format('{}.shp'.format(outbase)))

#Merge all these layers
tigerlist = [os.path.join(dirpath, file)
              for (dirpath, dirnames, filenames) in os.walk(tiger16dir)
              for file in filenames if re.search('tl_2016_[0-9]{5}_roads.shp$', file)]
arcpy.Merge_management(tigerlist, tigerroads)

#Subselect to only keep local roads — MTFCC IN ('S1400','S1640) https://www2.census.gov/geo/pdfs/reference/mtfccs2018.pdf
arcpy.MakeFeatureLayer_management(tigerroads, 'tigerroads_lyr', where_clause = "MTFCC IN ('S1400','S1640')")
arcpy.GetCount_management('tigerroads_lyr')
arcpy.CopyFeatures_management('tigerroads_lyr', tigerroads_sub)

#Intersect with urbanized area boundaries to add urban code (but can't find adjusted urban areas)
#Download tiger 2016 urban area data (not adjusted)
UA2016_url = 'http://www2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_ua10_500k.zip'
dlfile(UA2016_url, outpath=tiger16dir)
UA2016 = '{}.shp'.format(os.path.join(tiger16dir,
                                      os.path.splitext(os.path.split(UA2016_url)[1])[0]))

#Download USDOT urban area boundary file (not adjusted)
UA_API = 'https://geo.dot.gov/server/rest/services/NTAD/Urbanized_Areas/MapServer/0/query'
#print("/server/rest" in requests.get(UA_API).headers.get('Path')) could make sure that it's a rest api
APIdownload(baseURL=UA_API,
            workspace=USDOTdir,
            basename='USDOT_UA',
            itersize=100,
            IDlist=[0,10000])
arcpy.Merge_management(arcpy.ListFeatureClasses('{}_*'.format('USDOT_UA')),output='{}'.format('USDOT_UA'))
for fc in arcpy.ListFeatureClasses('{}_*'.format('USDOT_UA')):
    arcpy.Delete_management(fc)
    print('Deleted {}...'.format(fc))

#It seems that the tiger 2016 boundary file is the one that fits the HPMS urban codes best, so use that
arcpy.SpatialJoin_analysis(tigerroads_sub, UA2016, tigerroads_UA, 'JOIN_ONE_TO_ONE',
                           join_type='KEEP_ALL', match_option='INTERSECT')

#Add state_code field by intersecting with state boundaries (https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2016&layergroup=States+%28and+equivalent%29)
arcpy.SpatialJoin_analysis(tigerroads_UA, states, tigerroads_format, 'JOIN_ONE_TO_ONE',
                           join_type='KEEP_ALL', match_option='INTERSECT')
arcpy.AddField_management(tigerroads_format, 'state_code', 'SHORT')
arcpy.CalculateField_management(tigerroads_format, 'state_code', 'int(!STATEFP!)', 'PYTHON')

#Assign mean AADT from USAADT format — add 'local road' (f_system = 7) to aadt_filled field
if 'aadt_filled' not in [f.name for f in arcpy.ListFields(tigerroads_format)]:
    arcpy.AddField_management(tigerroads_format, field_name = 'aadt_filled', field_type='LONG')
else:
    print('aadt_filled field already exists...')

with arcpy.da.UpdateCursor(tigerroads_format, ['UATYP10', 'state_code', 'aadt_filled', 'OBJECTID']) as cursor:
    for row in cursor:
        print(row[3])
        if row[0] in ['C', 'U']: #If in urban cluster or urban area
            row[2] = USAADTformat[(USAADTformat['urban'] == 'URBAN') &
                                  (USAADTformat['state_code'] == row[1]) &
                                  (USAADTformat['F_SYSTEM'] == 7)]['AADT'][0]
        else:
            row[2] = USAADTformat[(USAADTformat['urban'] == 'RURAL') &
                                  (USAADTformat['state_code'] == row[1]) &
                                  (USAADTformat['F_SYSTEM'] == 7)]['AADT'][0]
        cursor.updateRow(row)

#Add f_system field
if 'f_system' not in [f.name for f in arcpy.ListFields(tigerroads_format)]:
    arcpy.AddField_management(tigerroads_format, field_name = 'f_system', field_type='LONG')
else:
    print('f_system field already exists...')
arcpy.CalculateField_management(tigerroads_format, 'f_system', 7, 'PYTHON')

#Merge with hpms
arcpy.Merge_management([hpms, tigerroads_format], hpmstiger)

########################################################################################################################
# COMPUTE FUNCTIONAL-CLASS BASED SPEED LIMIT AVERAGE FOR EVERY STATE
########################################################################################################################
#Compute % of records that already have AADT and speed limit data
totalroadlength = sum([row[0] for row in arcpy.da.SearchCursor(hpmstiger, ['Shape_Length'])])
arcpy.MakeFeatureLayer_management(hpms, 'hpms_datlyr')
arcpy.SelectLayerByAttribute_management('hpms_datlyr', 'NEW_SELECTION', 'aadt > 0')
aadtroadlength = sum([row[0] for row in arcpy.da.SearchCursor('hpms_datlyr', ['Shape_Length'])])
print(aadtroadlength/totalroadlength) #% of road length with aadt attribute

arcpy.SelectLayerByAttribute_management('hpms_datlyr', 'NEW_SELECTION', 'speed_limit > 0')
spdlroadlength = sum([row[0] for row in arcpy.da.SearchCursor('hpms_datlyr', ['Shape_Length'])])
print(spdlroadlength/totalroadlength) #% of road length with spd limit attribute

#Get average and median speed limit for each functional class in every state
arcpy.MakeFeatureLayer_management(hpms, 'hpmslyr')
query = "speed_limit > 0"
arcpy.SelectLayerByAttribute_management('hpmslyr', 'NEW_SELECTION', where_clause=query)
arcpy.GetCount_management('hpmslyr')
spdltab = pd.DataFrame(
    arcpy.da.TableToNumPyArray('hpmslyr', field_names=('state_code', 'f_system', 'urban_code','speed_limit')))
spdltab['urban'] = np.where(spdltab['urban_code'].isin([0, 99999]), 'rural', 'urban')

#Get average and median speed limit for each functional class in every state
spdlstats = spdltab.groupby(['state_code', 'f_system', 'urban']). \
    aggregate({'speed_limit':['mean','median','count']}). \
    unstack(level=[1,2]).\
    stack(level=[2,3], dropna=False). \
    reset_index()
spdlstats.columns = [''.join(col).strip() for col in spdlstats.columns.values]

#Fill in values for all f_class in case simply no record of speed for that f_class in that state
def wavg(group):
    return(np.sum(group.speed_limitmedian*group.speed_limitcount)/group.speed_limitcount.sum())
spdlstats_US = (5*np.floor(spdlstats.groupby(['f_system','urban']).apply(wavg)/5)).reset_index() #Compute weighted average of median at US level floored to the closest 5 mph
spdlstats_US.rename(columns={0: 'meanUSmedian'}, inplace=True)
spdlstats_filled = spdlstats.merge(spdlstats_US, on=['f_system','urban'])
spdsub = (spdlstats_filled['speed_limitmedian'].isna()) | (spdlstats_filled['speed_limitcount']<50)
spdlstats_filled.loc[spdsub, 'speed_limitmedian'] = spdlstats_filled.loc[spdsub, 'meanUSmedian']

#Assign median state-specific 'local roads' speed limit to those road segments without an f_system attribute
spdlstats_filled.loc[spdlstats['f_system']==0, 'speed_limitmedian'].replace(
    spdlstats_filled.loc[spdlstats_filled['f_system']==1, 'speed_limitmedian'])

#Compute spdl_filled for hpms
if 'spdl_filled' not in [f.name for f in arcpy.ListFields(hpmstiger)]:
    arcpy.AddField_management(hpmstiger, field_name = 'spdl_filled', field_type='SHORT')
else:
    print('spdl_filled field already exists...')

with arcpy.da.UpdateCursor(hpmstiger, ['speed_limit', 'spdl_filled', 'OBJECTID', 'urban_code', 'UATYP10',
                                       'state_code', 'f_system']) as cursor:
    for row in cursor:
        if row[0] > 0: #If already speed limit data
            row[1] = row[0] #Use that value for spdl_filled
        else: #If no speed limit data
            print(row[2])
            if row[3] not in [None, '0', '99999'] or row[4] in ['U', 'C']: #If urban
                row[1] = spdlstats_filled.loc[(spdlstats_filled['urban'] == 'urban') &
                                          (spdlstats_filled['state_code'] == row[5]) &
                                          (spdlstats_filled['f_system'] == row[6]), 'speed_limitmedian'].values[0]
            else: #if rural
                row[1] = spdlstats_filled.loc[(spdlstats_filled['urban'] == 'rural') &
                                          (spdlstats_filled['state_code'] == row[5]) &
                                          (spdlstats_filled['f_system'] == row[6]), 'speed_limitmedian'].values[0]
        cursor.updateRow(row)

#Project hpmstiger data
arcpy.Project_management(hpmstiger, hpmstigerproj, out_coor_system= cs_ref)

########################################################################################################################
# PREPARE US-WIDE TRANSIT DATA FROM NTM TO CREATE HEATMAPS
# Note: The National Transit Map has many glitches and is not a comprehensive datasets. Many shapes do not have
# corresponding records in the 'trips' table (linked through AgencyName and shape_id) because either shape_id is null,
# missing altogether, or in the wrong format (e.g. 200062.0 -- formatted as a numeric -- rather than 0200062) leading to
# failures in the joining process. The number of issues is so high that we did not attempt to fix it. This concerns
#  ~3780/63000 records
########################################################################################################################
#Download tables and shapefile directly from API
basename='NTM_shapes'
if not arcpy.Exists(os.path.join(NTMdir, '{}.shp'.format(basename))):
    APIdownload(baseURL="https://geo.dot.gov/server/rest/services/NTAD/GTFS_NTM/MapServer/1/query",
                workspace = NTMdir,
                basename = basename,
                itersize = 1000,
                IDlist = [0,1000000],
                geometry=True)
    arcpy.Merge_management(arcpy.ListFeatureClasses('{}_*'.format(basename)),output='{}'.format(basename))
    for fc in arcpy.ListFeatureClasses('{}_*'.format(basename)):
        arcpy.Delete_management(fc)
        print('Deleted {}...'.format(fc))

#Calendar dates
outcaldates = os.path.join(NTMdir, 'NTMAPI_calendar_dates.csv')
if not arcpy.Exists(outcaldates):
    APIdownload(baseURL="https://geo.dot.gov/server/rest/services/NTAD/GTFS_NTM/MapServer/2/query",
                workspace = NTMdir, basename = 'NTM_calendar_dates', itersize = 1000, IDlist = [0,100000], geometry=False)
    mergedel(NTMdir, 'NTM_calendar_dates_[0-9]+_[0-9]+[.]csv$', outcaldates, verbose=True)
#Trips
outtrips = os.path.join(NTMdir, 'NTMAPI_trips.csv')
if not arcpy.Exists(outtrips):
    APIdownload(baseURL="https://geo.dot.gov/server/rest/services/NTAD/GTFS_NTM/MapServer/7/query",
                workspace = NTMdir, basename = 'NTM_trips', itersize = 1000, IDlist = [0,1700000], geometry=False)
    mergedel(NTMdir, 'NTM_trips_[0-9]+_[0-9]+[.]csv$', os.path.join(NTMdir, 'NTMAPI_trips.csv'), verbose=True)
#Calendar
outcalendar = os.path.join(NTMdir, 'NTMAPI_calendar.csv')
if not arcpy.Exists(outcalendar):
    APIdownload(baseURL="https://geo.dot.gov/server/rest/services/NTAD/GTFS_NTM/MapServer/9/query",
                workspace = NTMdir, basename = 'NTM_calendar', itersize = 500, IDlist = [0,10000], geometry=False)
    mergedel(NTMdir, 'NTM_calendar_[0-9]+_[0-9]+[.]csv$', outcalendar, verbose=True)
#Routes
outroutes = os.path.join(NTMdir, 'NTMAPI_routes.csv')
if not arcpy.Exists(outroutes):
    APIdownload(baseURL="https://geo.dot.gov/server/rest/services/NTAD/GTFS_NTM/MapServer/11/query",
                workspace = NTMdir, basename = 'NTM_routes', itersize = 1000, IDlist = [0,15000], geometry=False)
    mergedel(NTMdir, 'NTM_routes_[0-9]+_[0-9]+[.]csv$', outroutes, verbose=True)

#Format data
GTFStoSHPweeklynumber(gtfs_dir= NTMdir, out_gdb=os.path.join(rootdir, 'results/NTM.gdb'), out_fc = 'NTM',
                      current=False, keep = True)

#Identify routes for which the shape does not correspond to the actual trajectory but just stops
arcpy.AddField_management(NTMroutes, 'vrtx_count', 'LONG')
arcpy.AddField_management(NTMroutes, 'vrtxlength_ratio', 'FLOAT')
fc = int(arcpy.GetCount_management(NTMroutes).getOutput(0))
with arcpy.da.UpdateCursor(NTMroutes, ['SHAPE@', 'vrtx_count', 'SHAPE@LENGTH', 'vrtxlength_ratio']) as cursor:
    x=0
    for row in cursor:
        if x%100 == 0:
            print('{}% of records processed'.format(100*x/fc))
        vcount = 0
        if row[0]: #Make sure that shape is valid (some 'None' shapes sometimes)
            for part in row[0]: #Iterate through parts for multipart features
                for pnt in part: #Iterate through part's points
                    if pnt:
                        vcount += 1
            row[1] = vcount
            row[3] = vcount/row[2] #vrtexlength_ratio = vertex count/shape_length

            cursor.updateRow(row)
        x += 1

#Only keep overground transports (although some subways might have some overground sections)
routeSQL = "(route_type IN ('0','2','3','5')) AND" \
           "(service_len_MIN > 30) AND" \
           "(adjustnum_SUM > 0) AND" \
           "(((NOT AgencyName = 'NJTRANSITBUS_20080_277_1046') AND (vrtxlength_ratio > 0.0003)) OR" \
           " ((AgencyName = 'NJTRANSITBUS_20080_277_1046') AND (vrtxlength_ratio > 0.001)))"
arcpy.MakeFeatureLayer_management(NTMroutes, 'routes_lyr', where_clause= routeSQL)
arcpy.CopyFeatures_management('routes_lyr', NTMsel)
NTMproj = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes_selproj')
arcpy.Project_management(in_dataset=NTMsel, out_dataset=NTMproj, out_coor_system=cs_ref)

#Create raster of weekly number of buses at the same resolution as bing data
# Convert weekly number of buses to integer
arcpy.AddField_management(NTMproj, 'adjustnum_int', 'SHORT')
with arcpy.da.UpdateCursor(NTMproj, ['adjustnum_SUM', 'adjustnum_int']) as cursor:
    for row in cursor:
        if row[0]:
            row[1] = int(10*row[0]+0.5)
            cursor.updateRow(row)

#Split lines at all intersections so that small identical overlapping segments can be dissolved
arcpy.SplitLine_management(NTMproj, NTMproj + '_split') #Split at intersection
arcpy.FindIdentical_management(NTMproj + '_split', "explFindID", "Shape") #Find overlapping segments and make them part of a group (FEAT_SEQ)
arcpy.MakeFeatureLayer_management(NTMproj + '_split', "intlyr")
arcpy.AddJoin_management("intlyr", arcpy.Describe("intlyr").OIDfieldName, "explFindID", "IN_FID", "KEEP_ALL")
arcpy.Dissolve_management("intlyr", NTMsplitdiss, dissolve_field='explFindID.FEAT_SEQ',
                          statistics_fields=[[os.path.split(NTMproj)[1] + '_split.adjustnum_int', 'SUM']]) #Dissolve overlapping segments
arcpy.RepairGeometry_management(NTMsplitdiss, delete_null = 'DELETE_NULL') #sometimes creates empty geom
#Get the length of a half pixel diagonal to create buffers for
#guaranteeing that segments potentially falling within the same pixel are rasterized separately
tolerance = (2.0**0.5)*float(restemplate.getOutput(0))/2
ExplodeOverlappingLines(NTMsplitdiss, tolerance)

#For each set of non-overlapping lines, create its own raster
tilef = 'expl'
tilelist = list(set([row[0] for row in arcpy.da.SearchCursor(NTMsplitdiss, [tilef])]))
outras_base = os.path.join(rootdir, 'results/transitnum')
arcpy.env.snapRaster = template_ras

for tile in tilelist:
    outras = outras_base + str(tile)
    if not arcpy.Exists(outras):
        try:
            selexpr = '{0} = {1}'.format(tilef, tile)
            #Delete intermediate stuff from previous iteration
            try:
                arcpy.Delete_management('transit_lyr')
                arcpy.ClearEnvironment('scratchWorkspace')
            except:
                pass

            #Create scratch workspace
            tmpdir = os.path.join(os.path.dirname(os.path.dirname(outras)),
                                  'tmp_{}'.format(str(os.path.basename(outras))))
            os.mkdir(tmpdir)
            arcpy.env.scratchWorkspace = tmpdir

            #Rasterize
            print(selexpr)
            arcpy.MakeFeatureLayer_management(NTMsplitdiss, 'transit_lyr', where_clause= selexpr)
            tmplyr = os.path.join(tmpdir, 'transit{}.shp'.format(tile))
            arcpy.CopyFeatures_management('transit_lyr', tmplyr)
            arcpy.PolylineToRaster_conversion(tmplyr, value_field='SUM_NTM_routes_selproj_split_adjustnum_int'[0:10],
                                              out_rasterdataset=outras, cellsize=restemplate)
        except:
            traceback.print_exc()
            try:
                if arcpy.Exists(tmplyr):
                    arcpy.Delete_management(tmplyr)
                os.rmdir(tmpdir)

            except:
                pass

        #Remove intermediate products
        print('Deleting scratch workspace...')
        arcpy.Delete_management(tmplyr)
        os.rmdir(tmpdir)

    else:
        print('{} already exists...'.format(outras))

#Mosaic to new raster
arcpy.env.workspace = os.path.split(outras_base)[0]
transitras_tiles = arcpy.ListRasters('transitnum_*')
NTMras = os.path.join(rootdir, 'results/transit.gdb/NTMras')
arcpy.MosaicToNewRaster_management(transitras_tiles, os.path.split(NTMras)[0], os.path.split(NTMras)[1],
                                   pixel_type='32_BIT_UNSIGNED', number_of_bands= 1, mosaic_method = 'SUM')
for tile in transitras_tiles:
    print('Deleting {}...'.format(tile))
    arcpy.Delete_management(tile)
arcpy.ClearEnvironment('Workspace')

########################################################################################################################
# PREPARE DATA ON ROAD GRADIENTS
########################################################################################################################
#Subset road dataset based on convex hull around sites
arcpy.Project_management(XRFsiteshull, XRFsiteshull_aea, out_coor_system=cs_ref)
arcpy.Clip_analysis(hpmstigerproj, XRFsiteshull_aea, hpms_sub)

#Compute statistics
arcpy.env.snapRaster = NED19proj
arcpy.PolylineToRaster_conversion(hpms_sub, 'OBJECTID', hpms_ras19, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'aadt_filled', cellsize = NED19proj)
ZonalStatisticsAsTable(hpms_ras19, 'Value', NED19proj, out_table = rangetab19,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(hpms_ras19, 'Value', NED19smooth, out_table = rangetab19_smooth,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.env.snapRaster = NED13proj
arcpy.PolylineToRaster_conversion(hpms_sub, 'OBJECTID', hpms_ras13, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'aadt_filled', cellsize = NED13proj)
ZonalStatisticsAsTable(hpms_ras13, 'Value', NED13proj, out_table = rangetab13,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')
ZonalStatisticsAsTable(hpms_ras13, 'Value', NED13smooth, out_table = rangetab13_smooth,
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.ClearEnvironment('snapRaster')

#Get all range values for zonal statistics tables
tablist =  [rangetab19 + '_smooth', rangetab13 + '_smooth']
elvdic = defaultdict(lambda: [-99.0] * len(tablist))
for i in range(0, len(tablist)):
    print(tablist[i])
    for row in arcpy.da.SearchCursor(tablist[i], ['Value', 'RANGE']):
        elvdic[row[0]][i] = row[1]

#Compute slope
arcpy.AddField_management(hpms_sub, 'gradient', 'FLOAT')

#Assign composite value to road dataset
missinglist = []
with arcpy.da.UpdateCursor(hpms_sub, ['OBJECTID', 'Shape_Length', 'gradient']) as cursor:
    for row in cursor:
        if elvdic[row[0]][0] != -99.0:
            row[2] = elvdic[row[0]][0]/row[1] # gradient19_smooth = RANGE_19smooth/Shape_Length
            cursor.updateRow(row)
        elif elvdic[row[0]][1] != -99.0:
            row[2] = elvdic[row[0]][1]/row[1] # gradient13_smooth = RANGE_13smooth/Shape_Length
            cursor.updateRow(row)
        else:
            print('There are no data for OBJECTID {}'.format(row[0]))
            row[2] = 0
            missinglist.append(row[0])
        cursor.updateRow(row)
if len(missinglist)>0:
    print('{} roads did not have any elevation range value'.format(len(missinglist)))

########################################################################################################################
# CREATE HEATMAPS FOR SEATTLE MODEL TRAINING (for all variables aside from Bing data)
########################################################################################################################
#Get heatmaps
arcpy.env.workspace = pollutgdb
arcpy.env.snapRaster = template_ras
#AADT
arcpy.PolylineToRaster_conversion(hpms_sub, value_field='aadt_filled', out_rasterdataset='hpmssubAADT',
                                  priority_field='aadt_filled', cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(pollutgdb, 'hpmssubAADT'),
              out_gdb=pollutgdb, out_var='subAADT', divnum=100, keyw='(pow|log)[1235]00(_[123])*', verbose=True)

#Speed limit
arcpy.PolylineToRaster_conversion(hpms_sub, value_field='spdl_filled', out_rasterdataset='hpmssubspdl',
                                  priority_field='spdl_filled', cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(pollutgdb, 'hpmssubspdl'),
              out_gdb=pollutgdb, out_var='subspdl', divnum=100, keyw='(pow|log)(([1235]00)|50)(_[123])*', verbose=True)

#Public transit
arcpy.env.extent = hpms_sub
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=PStransitras,
              out_gdb=pollutgdb, out_var='bustransit', divnum=100,
              keyw='(pow|log)[1235]00(_[123])*', verbose=True)

#Road gradient
arcpy.PolylineToRaster_conversion(hpms_sub, value_field='gradient', out_rasterdataset='hpmssubslope',
                                  priority_field='aadt_filled', cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(pollutgdb, 'hpmssubslope'),
              out_gdb=pollutgdb, out_var='subslope', divnum=0.01, keyw='(pow|log)(([1235]00)|50)(_[123])*', verbose=True)

#Bing
#See Bing_format_PS.py (for Sound) and Bing_format_parallel.py (for AQI stations)

########################################################################################################################
# GET TREES HEATMAP VALUES
# Select candidate species of trees from the City of Seattle's street-tree dataset and extract heatmap values at their
# location
########################################################################################################################
#Get heat values for all trees
heatlist = Iter_ListRaster([pollutgdb, transitgdb], 'heat*') + glob.glob(os.path.join(bingdir, 'heat*bing*proj.tif')) +\
                       [NLCD_reclass_PS, NLCD_imp_PS, NLCD_imp_PS + '_mean.tif']

#Make sure of the projection of all layers
# for heatras in heatlist:
#     print(heatras)
#     print(arcpy.Describe(heatras).SpatialReference.name)

#Project
arcpy.Project_management(XRFsites, XRFsites_aea, cs_ref)

#Extract values
arcpy.env.qualifiedFieldNames = False
ExtractMultiValuesToPoints(XRFsites_aea, heatlist, bilinear_interpolate_values='NONE')
arcpy.AddGeometryAttributes_management(XRFsites_aea, 'POINT_X_Y_Z_M', Coordinate_System= arcpy.SpatialReference(4326))

arcpy.AddField_management(XRFsites_aea, 'SiteIDPair', field_type='TEXT')
arcpy.CalculateField_management(XRFsites_aea, 'SiteIDPair', expression='!F_!+!Pair!', expression_type='PYTHON')

########################################################################################################################
# CREATE HEATMAPS FOR PUGET SOUND MODEL PREDICTIONS (subset of variables based on selected model)
########################################################################################################################
#Subset hpms tiger road dataset based on Puget Sound
arcpy.Clip_analysis(hpmstigerproj, PSdissolve, hpms_PS)

arcpy.env.workspace = pollutgdbPS

#AADT
arcpy.env.snapRaster = template_ras
arcpy.PolylineToRaster_conversion(hpms_PS, value_field='aadt_filled', out_rasterdataset='hpmsPSAADT',
                                  priority_field='aadt_filled', cellsize=restemplate)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(pollutgdbPS, 'hpmsPSAADT'),
              out_gdb=pollutgdbPS, out_var='PSAADT', divnum=100, keyw='log100*', verbose=True)
arcpy.ClearEnvironment("snapRaster")

########################################################################################################################
# SUBSET HPMS TIGER AND GET AADT HEATMAP VALUES AROUND AQI STATIONS
########################################################################################################################
#Intersect roads with buffers
arcpy.Intersect_analysis([AQIsites_bufunion, hpmstigerproj], roadAQI, join_attributes='ALL')

############################## TO RUN ##################################################################################
##########Run AQI_AADTheatmap.py##################

#For each station, extract aadt_filled
heataadt_list = [os.path.join(dirpath, file)
              for (dirpath, dirnames, filenames) in os.walk(AQIdir)
              for file in filenames if re.search('heataadt_filled.*[.]tif$', file)]

#Get extent of every tile
heataadt_dict = {}
for tile in heataadt_list:
    tilext = arcpy.Describe(tile).Extent
    heataadt_dict[tile] = [tilext.XMin, tilext.XMax, tilext.YMin, tilext.YMax]

#For each AQI station, find the corresponding aadt heatmap tile and extract value to table
tablist = []
if 'aadtlog100' not in [f.name for f in arcpy.ListFields(AQIsites)]:
    arcpy.AddField_management(AQIsites, 'aadtlog100', 'DOUBLE')

with arcpy.da.UpdateCursor(AQIsites, ['OID@','SHAPE@XY', 'aadtlog100']) as cursor:
    for row in cursor:
        print(row[0])
        for k,v in heataadt_dict.iteritems(): #For each aadt heatmap
            if row[1][0] > v[0] and row[1][0] < v[1]: #If station X > tile XMin and station X < tile XMax
                if row[1][1] > v[2] and row[1][1] < v[3]: #If station Y > tile YMin and station Y < tile YMax
                    # outab = os.path.join(AQIdir, 'heataadt_stationID{0}'.format(row[0]))
                    # print(outab)
                    # tablist.append(outab)
                    outsamp = arcpy.GetCellValue_management(
                        in_raster=k,
                        location_point = "{0} {1}".format(row[1][0], row[1][1])).getOutput(0)

                    if outsamp == 'NoData':
                        row[2] = 0.0
                    else:
                        row[2] = outsamp
                    cursor.updateRow(row)
                    # outsamp = arcpy.sa.Sample(in_rasters=k, in_location_data=row[2], out_table=outab,
                    #                 resampling_type= 'BILINEAR')


arcpy.Merge_management(tablist, output = AQIaadttab)

########################################################################################################################
# PREPARE LAND USE DATA FOR AQI STATIONS
########################################################################################################################
#Compute focal stats within AQI buffers
arcpy.env.mask = AQIsites_bufunion
arcpy.env.snapRaster = NLCD_imp
imp_mean = arcpy.sa.FocalStatistics(NLCD_imp, neighborhood = NbrCircle(3, "CELL"), statistics_type= 'MEAN')
imp_mean.save(NLCD_imp_AQI_mean)

#STILL TO DO - 20191119
#Extract smoothed imperviousness values at AQI stations
arcpy.sa.ExtractMultiValuesToPoints(AQIsites, NLCD_imp_AQI_mean, bilinear_interpolate_values = 'NONE')