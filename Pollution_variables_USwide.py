import arcpy
import pandas as pd
import geopandas as gpd
import os
import re
import itertools
import shapely
from BeautifulSoup import BeautifulSoup
import urllib2
import zipfile
import numpy as np
from SpatialJoinLines_LargestOverlap import *

pd.options.display.max_columns = 20
arcpy.env.overwriteOutput = True

rootdir = "D:\Mathis\ICSL\stormwater"
USDOTdir = os.path.join(rootdir, "data\USDOT_0319")
resdir = os.path.join(rootdir, 'results/usdot')
USOSMdir = os.path.join(rootdir, 'data/OSM_20190220/')
AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')
UA = os.path.join(USDOTdir, '2010CensusUrbanizedArea/UrbanizedArea2010.shp')

#Output variables
usdotgdb = os.path.join(resdir, 'usdot.gdb')
ushpms = os.path.join(usdotgdb, 'ushpms')
USOSM = os.path.join(AQIgdb, 'USOSM')
USOSM_UA = os.path.join(AQIgdb, 'USOSM_UA')
OSMAQI =  os.path.join(AQIgdb, 'AQI_OSMinters')
OSMAQIproj = os.path.join(AQIgdb, 'OSMAQIproj')
OSMhpms = os.path.join(usdotgdb, 'OSMhpms_join')


if not arcpy.Exists(usdotgdb):
    arcpy.CreateFileGDB_management(resdir, out_name = 'usdot')
    
########################################################################################################################
# FORMAT HIGHWAY AADT SPATIAL DATA
########################################################################################################################
#Function to download file
def dlfile(url, outpath):
    # Open the url
    try:
        f = urllib2.urlopen(url)
        print "downloading " + url

        # Open our local file for writing
        with open(outpath, "wb") as local_file:
            local_file.write(f.read())

    #handle errors
    except urllib2.HTTPError, e:
        print "HTTP Error:", e.code, url
    except urllib2.URLError, e:
        print "URL Error:", e.reason, url

    #Unzip folder
    if zipfile.is_zipfile(outpath):
        print('Unzipping {}...'.format(os.path.split(outpath)[1]))
        with zipfile.ZipFile(outpath) as zipf:
            zipfilelist = [info.filename for info in zipf.infolist()]
            listcheck = [f for f in zipfilelist if os.path.exists(os.path.join(outpath, f))]
            if len(listcheck) > 0:
                print('Overwriting {}...'.format(', '.join(listcheck)))
            zipf.extractall(os.path.split(outpath)[0])
        del zipf



########################################################################################################################
# COMPUTE FUNCTIONAL-CLASS BASED AADT AVERAGE FOR EVERY STATE
########################################################################################################################
#Get VM-2 estimates of Vehicle Miles Traveled for each state and average it for small roads at state level (idea fromhttps://www.pnas.org/content/pnas/suppl/2015/04/01/1421723112.DCSupplemental/pnas.1421723112.sapp.pdf)
#Highway Statistics Book 2017
# Public road lenth 2017 - miles by functional system, Table HM-20: https://www.fhwa.dot.gov/policyinformation/statistics/2017/hm20.cfm
hm20url = "https://www.fhwa.dot.gov/policyinformation/statistics/2017/hm20.cfm"
hm20_page = urllib2.urlopen(hm20url)
hm20_soup = BeautifulSoup(hm20_page)
hm20_regex = re.compile(".*hm20.*[.]xls")
hm20urltab = os.path.split(hm20url)[0] + '/' + \
             '/'.join([i for i in hm20_soup.find('a', attrs={'href': hm20_regex}).get('href').split('/')
                       if i not in os.path.split(hm20url)[0].split('/')])
hm20tab = os.path.join(USDOTdir, os.path.split(hm20urltab)[1])
dlfile(hm20urltab, hm20tab)

# Functional system travel 2017 - Annual vehicle-miles, Table VM-2: https://www.fhwa.dot.gov/policyinformation/statistics/2017/vm2.cfm
vm2url = "https://www.fhwa.dot.gov/policyinformation/statistics/2017/vm2.cfm"
vm2_page = urllib2.urlopen(vm2url)
vm2_soup = BeautifulSoup(vm2_page)
vm2_regex = re.compile(".*vm2.*[.]xls")
vm2urltab = os.path.split(vm2url)[0] + '/' + \
             '/'.join([i for i in vm2_soup.find('a', attrs={'href': vm2_regex}).get('href').split('/')
                       if i not in os.path.split(vm2url)[0].split('/')])
vm2tab = os.path.join(USDOTdir, os.path.split(vm2urltab)[1])
dlfile(vm2urltab, vm2tab)

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
USAADT[USAADT != 0] = \
    ((10**6)*vm2sheet_format[vm2sheet_format!=0].iloc[:,1:])\
    .div(hm20sheet_format[hm20sheet_format!=0].iloc[:,1:])\
    .div(365)

########################################################################################################################
# MERGE OSM WITH HPMS
########################################################################################################################
#-----------------------------------------------------------------------------------------------------------------------
# FORMAT OSM DATA FOR ANALYSIS
#-----------------------------------------------------------------------------------------------------------------------
#Merge OSM data for US
USOSMlist = []
for (dirname, dirs, files) in os.walk(USOSMdir):
   for filename in files:
       if filename == 'gis_osm_roads_free_1.shp':
           print(os.path.split(dirname)[1])
           USOSMlist.append(os.path.join(dirname,filename))
arcpy.Merge_management(USOSMlist, USOSM)

#Only select road types with car traffic
#Select OSM roads, but this map with service roads, as enough through traffic to potentially have some impact
arcpy.MakeFeatureLayer_management(USOSM, 'OSMroads_lyr')
np.unique([row[0] for row in arcpy.da.SearchCursor('OSMroads_lyr', ['fclass'])])
sel = "{} IN ('motorway','motorway_link','living_street','primary','primary_link','residential','secondary','secondary_link'," \
      "'tertiary','tertiary_link','trunk','trunk_link','unclassified','unknown', 'service')".format('"fclass"')
arcpy.SelectLayerByAttribute_management('OSMroads_lyr', 'NEW_SELECTION', sel)
arcpy.CopyFeatures_management('OSMroads_lyr', USOSM + 'sub')
arcpy.Delete_management('OSMroads_lyr')

#Overlay OSM with urbanized areas
arcpy.SpatialJoin_analysis(USOSM+'sub', UA, USOSM_UA, join_operation="JOIN_ONE_TO_ONE", join_type='KEEP_ALL',
                           match_option='INTERSECT')

#Merge with OSM by spatial overlap
SpatialJoinLines_LargestOverlap(target_features= USOSM, join_features= ushpms,
                                out_fc= os.path.split(OSMhpms)[1], outgdb=os.path.split(OSMhpms)[0],
                                bufsize='10 meters', keep_all=True,
                                fields_select=[f.name for f in arcpy.ListFields(ushpms)])


# sites_outbufunion = os.path.join(outdir, 'airsites_600bufunion.shp')
#Project OSM data
# arcpy.Project_management(OSMAQI, OSMAQIproj, out_coor_system= sites_outbufunion
# #Intersect roads with buffers
# arcpy.Intersect_analysis([sites_outbufunion, USOSM], OSMAQI, join_attributes='ALL')
#
# #Get SPD and AADT medians for each fclass
# fclass_SPDADTmedian = pd.read_csv(os.path.join(rootdir, 'results/OSM_SPDADTmedian.csv'))
#
# #OSM AADT AND SPD (apply Puget Sound data for now  get US-wide data in the future)
# #Convert OSM functional categories to numbers
# arcpy.AddField_management(OSMAQI, 'fclassADT', 'LONG')
# arcpy.AddField_management(OSMAQI, 'fclassSPD', 'LONG')
# with arcpy.da.UpdateCursor(OSMAQI, ['fclass','fclassADT', 'fclassSPD']) as cursor:
#     for row in cursor:
#         if any(fclass_SPDADTmedian['first'].isin([row[0]])):
#             if row[0] in ['service','unclassified','unknown', 'living_street']: #Value for service and unclassified seem overestimated
#                 row[1] = int(fclass_SPDADTmedian[fclass_SPDADTmedian['first'] == 'residential'].mean_ADT)
#             else:
#                 row[1] = int(fclass_SPDADTmedian.loc[fclass_SPDADTmedian['first'] == row[0]].mean_ADT)
#
#              row[2] = int(fclass_SPDADTmedian.loc[fclass_SPDADTmedian['first'] == row[0]].mean_SPD)
#         else:
#             row[1]=0
#             row[2]=0
#         cursor.updateRow(row)

#Donwload National Transit Map data