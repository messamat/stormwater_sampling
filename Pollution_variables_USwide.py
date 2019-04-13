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
PSgdb=os.path.join(rootdir,'results/PSOSM.gdb')
AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')

NED19proj = os.path.join(rootdir, 'results/ned19_psproj')
NED13proj = os.path.join(rootdir, 'results/ned13_psproj')

XRFsites = os.path.join(rootdir, 'data/field_data/sampling_sites_edit_select.shp')

NLCD_reclass = os.path.join(rootdir, 'results/LU.gdb/NLCD_reclass_final')
NLCD_imp = os.path.join(rootdir, 'data/nlcd_2011_impervious_2011_edition_2014_10_10/nlcd_2011_impervious_2011_edition_2014_10_10.img') #Zipped
UTM10 = arcpy.SpatialReference(26910)

STgtfs = os.path.join(rootdir, 'data\SoundTransit_201812\gtfs_puget_sound_consolidated.zip')
transitwiki_dir = os.path.join(rootdir, 'data/TransitWiki201812')

template_ras = os.path.join(rootdir,'results/bing/181204_02_00_class_mlc.tif')
res = arcpy.GetRasterProperties_management(template_ras, 'CELLSIZEX')
pollutgdb = os.path.join(rootdir,'results/pollution_variables.gdb')
if arcpy.Exists(pollutgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results'), 'pollution_variables.gdb')

#usa contiguous albers equal area conic - esri projection 102003
#Spatial reference: http://spatialreference.org/ref/esri/usa-contiguous-albers-equal-area-conic/
albers = arcpy.SpatialReference(102003) #"+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Output variables
usdotgdb = os.path.join(resdir, 'usdot.gdb')
hpms = os.path.join(USDOTdir, '2016.gdb/HPMS2016')
hmps_spdl = os.path.join(usdotgdb, 'hpms_spdl')

NED19smooth = os.path.join(rootdir, 'results/ned19_smooth')
NED13smooth = os.path.join(rootdir, 'results/ned13_smooth')
rangetab19 = os.path.join(PSgdb, 'hpms_elv19range')
rangetab13 = os.path.join(PSgdb, 'hpms_elv13range')
rangetab19_smooth = os.path.join(PSgdb, 'hpms_elv19range_smooth')
rangetab13_smooth = os.path.join(PSgdb, 'hpms_elv13range_smooth')

tiger16dir = os.path.join(rootdir, 'data/TIGER2016')
tigerroads = os.path.join(usdotgdb, 'tigerroads_merge')
tigerroads_sub = os.path.join(usdotgdb, 'tigerroads_sub')
tigerroads_UA = os.path.join(usdotgdb, 'tigerroads_UA')
tigerroads_format = tigerroads_UA + 'states'
states = os.path.join(tiger16dir, 'tl_2016_us_state/tl_2016_us_state.shp')
hpmstiger = os.path.join(usdotgdb, 'hpmstiger')

hpms_sub = os.path.join(usdotgdb, 'hpmstiger_subset')
hpms_ras19 = os.path.join(usdotgdb, 'hpmstiger_subsetras19')
hpms_ras13 = os.path.join(usdotgdb, 'hpmstiger_subsetras13')
siteshull = os.path.join(rootdir, 'results/XRFsites_convexhull.shp')

NLCD_reclass_PS = os.path.join(rootdir, 'results/NLCD_reclass_final_PS.tif')
NLCD_imp_PS = os.path.join(rootdir, 'results/nlcd_imp_ps')

NTMroutes = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes')
NTMsel = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes_sel')
NTMproj = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes_selproj')
NTMsplitdiss = NTMproj + '_splitv_diss'
NTMras = os.path.join(rootdir, 'results/transit.gdb/NTMras')

soundtransit = os.path.join(rootdir, 'results/transit.gdb/SoundTransit')
PStransit = os.path.join(rootdir, 'results/transit.gdb/PStransit')
PStransitbus = PStransit + '_busroutes'
PStransitbus_proj = PStransit + '_busroutes_proj'
PStransitbus_splitdiss = PStransitbus_proj + '_splitv_diss'
PStransitras = os.path.join(rootdir, 'results/transit.gdb/PStransit_ras')

if not arcpy.Exists(usdotgdb):
    arcpy.CreateFileGDB_management(resdir, out_name = 'usdot')

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

########################################################################################################################
# PREPARE TRANSIT DATA FROM NTM TO CREATE HEATMAPs
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
albers = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
NTMproj = os.path.join(rootdir, 'results/NTM.gdb/NTM_routes_selproj')
arcpy.Project_management(in_dataset=NTMsel, out_dataset=NTMproj, out_coor_system=albers)

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
##########RUN###########
arcpy.FindIdentical_management(NTMproj + '_split', "explFindID", "Shape") #Find overlapping segments and make them part of a group (FEAT_SEQ)
arcpy.MakeFeatureLayer_management(NTMproj + '_split', "intlyr")
arcpy.AddJoin_management("intlyr", arcpy.Describe("intlyr").OIDfieldName, "explFindID", "IN_FID", "KEEP_ALL")
arcpy.Dissolve_management("intlyr", NTMsplitdiss, dissolve_field='explFindID.FEAT_SEQ',
                          statistics_fields=[[os.path.split(NTMproj)[1] + '_split.adjustnum_int', 'SUM']]) #Dissolve overlapping segments
arcpy.RepairGeometry_management(NTMsplitdiss, delete_null = 'DELETE_NULL') #sometimes creates empty geom
#Get the length of a half pixel diagonal to create buffers for
#guaranteeing that segments potentially falling within the same pixel are rasterized separately
tolerance = (2.0**0.5)*float(res.getOutput(0))/2
ExplodeOverlappingLines(NTMsplitdiss, tolerance)

#For each set of non-overlapping lines, create its own raster
tilef = 'expl'
tilelist = list(set([row[0] for row in arcpy.da.SearchCursor(NTMsplitdiss, [tilef])]))
outras_base = os.path.join(rootdir, 'results/NTM.gdb/transitnum_')
arcpy.env.snapRaster = template_ras
for tile in tilelist:
    outras = outras_base + str(tile)
    if not arcpy.Exists(outras):
        selexpr = '{0} = {1}'.format(tilef, tile)
        print(selexpr)
        arcpy.MakeFeatureLayer_management(NTMsplitdiss, 'transit_lyr', where_clause= selexpr)
        arcpy.PolylineToRaster_conversion('transit_lyr', value_field='adjustnum_int', out_rasterdataset=outras, cellsize=res)

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
            logger.addHandler(fh) # add the handler to the logger
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
arcpy.Project_management(PStransitbus, PStransitbus_proj, UTM10)

#Create raster of weekly number of buses at the same resolution as bing data
# Convert weekly number of buses to integer
arcpy.AddField_management(PStransitbus_proj, 'adjustnum_int', 'SHORT')
arcpy.CalculateField_management(PStransitbus_proj, 'adjustnum_int',
                                expression='int(10*!SUM_adjustnum!+0.5)', expression_type='PYTHON')

#Split lines at all intersections so that small identical overlapping segments can be dissolved
arcpy.SplitLine_management(PStransitbus_proj, PStransitbus_proj + '_split') #Split at intersection
arcpy.FindIdentical_management(PStransitbus_proj + '_split', "explFindID", "Shape") #Find overlapping segments and make them part of a group (FEAT_SEQ)
arcpy.MakeFeatureLayer_management(PStransitbus_proj + '_split', "intlyr")
arcpy.AddJoin_management("intlyr", arcpy.Describe("intlyr").OIDfieldName, "explFindID", "IN_FID", "KEEP_ALL")
arcpy.Dissolve_management("intlyr", PStransitbus_splitdiss, dissolve_field='explFindID.FEAT_SEQ',
                          statistics_fields=[[os.path.split(PStransitbus_proj)[1] + '_split.adjustnum_int', 'SUM']]) #Dissolve overlapping segments
arcpy.RepairGeometry_management(PStransitbus_splitdiss, delete_null = 'DELETE_NULL') #sometimes creates empty geom
#Get the length of a half pixel diagonal to create buffers for
#guaranteeing that segments potentially falling within the same pixel are rasterized separately
tolerance = (2.0**0.5)*float(res.getOutput(0))/2
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


########################################################################################################################
# PREPARE DATA ON ROAD GRADIENTS
########################################################################################################################
#Subset road dataset based on convex hull around sites
arcpy.Buffer_analysis(XRFsites, os.path.join(rootdir, 'results/XRFsites_3kmbuf.shp'),
                      buffer_distance_or_field='3000 meters', method='GEODESIC')
arcpy.MinimumBoundingGeometry_management(os.path.join(rootdir, 'results/XRFsites_3kmbuf.shp'),
                                         siteshull, geometry_type='CONVEX_HULL', group_option='ALL')
arcpy.Clip_analysis(hpmstiger, siteshull, hpms_sub)

#Compute statistics
[f.name for f in arcpy.ListFields(hpms_sub)]
arcpy.PolylineToRaster_conversion(hpms_sub, 'OBJECTID', hpms_ras19, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'aadt_filled', cellsize = NED19proj)
ZonalStatisticsAsTable(hpms_ras19, 'Value', NED19smooth, out_table = rangetab19 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

arcpy.PolylineToRaster_conversion(hpms_sub, 'OBJECTID', hpms_ras13, cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'aadt_filled', cellsize = NED13proj)
ZonalStatisticsAsTable(hpms_ras13, 'Value', NED13smooth, out_table = rangetab13 + '_smooth',
                       statistics_type= 'RANGE', ignore_nodata='NODATA')

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
# PREPARE LAND USE DATA
########################################################################################################################
#Export NLCD data to Puget Sound scale
arcpy.env.extent = siteshull
arcpy.ProjectRaster_management(NLCD_reclass, NLCD_reclass_PS, UTM10, resampling_type='NEAREST')
#Export NLCD impervious data
arcpy.ProjectRaster_management(NLCD_imp, NLCD_imp_PS, UTM10, resampling_type='BILINEAR')
#Compute focal stats
imp_mean = arcpy.sa.FocalStatistics(NLCD_imp_PS, neighborhood = NbrCircle(3, "CELL"), statistics_type= 'MEAN')
imp_mean.save(NLCD_imp_PS + '_mean.tif')

########################################################################################################################
# CREATE HEATMAPS FOR SEATTLE AREA MODEL TRAINING (for all variables aside from Bing data)
########################################################################################################################
arcpy.env.workspace = pollutgdb
#AADT
arcpy.PolylineToRaster_conversion(hpms_sub, value_field='aadt_filled', out_rasterdataset='hpmssubAADT',
                                  priority_field='aadt_filled', cellsize=res)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(pollutgdb, 'hpmssubAADT'),
              out_gdb=pollutgdb, out_var='subAADT', divnum=100, keyw='(pow|log)[1235]00(_[123])*', verbose=True)

#Speed limit
[f.name for f in arcpy.ListFields(hpms_sub)]
arcpy.PolylineToRaster_conversion(hpms_sub, value_field='spdl_filled', out_rasterdataset='hpmssubspdl',
                                  priority_field='spdl_filled', cellsize=res)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(pollutgdb, 'hpmssubspdl'),
              out_gdb=pollutgdb, out_var='subspdl', divnum=100, keyw='(pow|log)(([1235]00)|50)(_[123])*', verbose=True)

#Public transit (

#Road gradient
arcpy.PolylineToRaster_conversion(hpms_sub, value_field='gradient', out_rasterdataset='hpmssubslope',
                                  priority_field='aadt_filled', cellsize=res)
customheatmap(kernel_dir=os.path.join(rootdir, 'results/bing'), in_raster=os.path.join(pollutgdb, 'hpmssubslope'),
              out_gdb=pollutgdb, out_var='subslope', divnum=0.01, keyw='(pow|log)(([1235]00)|50)(_[123])*', verbose=True)

#Bing
