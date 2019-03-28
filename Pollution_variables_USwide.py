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
import us
import itertools
from SpatialJoinLines_LargestOverlap import *
from GTFStoSHP import *
import shutil
import contextlib
import ftplib
import urlparse


pd.options.display.max_columns = 20
arcpy.env.overwriteOutput = True

rootdir = "D:\Mathis\ICSL\stormwater"
USDOTdir = os.path.join(rootdir, "data\USDOT_0319")
NTMdir = os.path.join(rootdir, "data\NTM_0319")
resdir = os.path.join(rootdir, 'results/usdot')
PSgdb=os.path.join(rootdir,'results/PSOSM.gdb')
AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')
UA = os.path.join(USDOTdir, '2010CensusUrbanizedArea/UrbanizedArea2010.shp')

NED19proj = os.path.join(rootdir, 'results/ned19_psproj')
NED13proj = os.path.join(rootdir, 'results/ned13_psproj')
NED19smooth = os.path.join(rootdir, 'results/ned19_smooth')
NED13smooth = os.path.join(rootdir, 'results/ned13_smooth')

#Output variables
usdotgdb = os.path.join(resdir, 'usdot.gdb')
hpms = os.path.join(USDOTdir, '2016.gdb/HPMS2016')
hmps_spdl = os.path.join(usdotgdb, 'hpms_spdl')
tiger16dir = os.path.join(rootdir, 'data/TIGER2016')
rangetab19 = os.path.join(PSgdb, 'OSM_elv19range')
rangetab13 = os.path.join(PSgdb, 'OSM_elv13range')
rangetab19_smooth = os.path.join(PSgdb, 'OSM_elv19range_smooth')
rangetab13_smooth = os.path.join(PSgdb, 'OSM_elv13range_smooth')
sroadsras = os.path.join(PSgdb, 'Seattle_roadras')
srangetab19 = os.path.join(PSgdb, 'Seattle_elv19range')
srangetab13 = os.path.join(PSgdb, 'Seattle_elv13range')

if not arcpy.Exists(usdotgdb):
    arcpy.CreateFileGDB_management(resdir, out_name = 'usdot')

#Function to download and unzip miscellaneous types of files
#Partly inspired from https://www.codementor.io/aviaryan/downloading-files-from-urls-in-python-77q3bs0un
def is_downloadable(url):
    """
    Does the url contain a downloadable resource
    """
    h = requests.head(url, allow_redirects=True)
    header = h.headers
    content_type = header.get('content-type')
    if 'html' in content_type.lower():
        return False
    return True

def get_filename_from_cd(url):
    """
    Get filename from content-disposition
    """
    r = requests.get(url, allow_redirects=True)
    cd = r.headers.get('content-disposition')
    if not cd:
        return None
    fname = re.findall('filename=(.+)', cd)
    if len(fname) == 0:
        return None
    return fname[0]

def dlfile(url, outpath, outfile=None, fieldnames=None):
    """Function to download file from URL path and unzip it.
    URL (required): URL of file to download
    outpath (required): the full path including
    outfile (optional): the output name without file extension, otherwise gets it from URL
    fieldnames (optional): fieldnames in output table if downloading plain text"""

    try:
        if is_downloadable(url): #check that url is not just html
            # Get output file name
            if outfile is None:
                outfile = get_filename_from_cd(url)
                if outfile is not None:
                    out = os.path.join(outpath, outfile)
                else:
                    out = os.path.join(outpath, os.path.split(url)[1])
            else:
                out = os.path.join(outpath, outfile + os.path.splitext(url)[1])
            del outfile

            #http request
            f = requests.get(url, allow_redirects=True)
            print "downloading " + url

            # Open local file for writing
            if not os.path.exists(out):
                if 'csv' in f.headers.get('content-type').lower(): #If csv file
                    df = pd.read_csv(io.StringIO(f.text))
                    df.to_csv(out, index=False)

                elif f.headers.get('content-type').lower() == 'text/plain': #If plain text
                    dialect = csv.Sniffer().sniff(f.text)
                    txtF = csv.DictReader(f.text.split('\n'),
                                          delimiter=dialect.delimiter,
                                          fieldnames=fieldnames)
                    with open(out, "wb") as local_file:
                        writer = csv.DictWriter(local_file, fieldnames=fieldnames)
                        writer.writeheader()
                        for row in txtF:
                            writer.writerow(row)

                else: #Otherwise, just try reading
                    with open(out, "wb") as local_file:
                        local_file.write(f.read())
            else:
                print('{} already exists...'.format(out))

    #handle errors
    except requests.exceptions.HTTPError, e:
        print "HTTP Error:", e.code, url
    except Exception as e:
        print e
        if os.path.exists(out):
            os.remove(out)

    #Unzip folder
    if zipfile.is_zipfile(out):
        print('Unzipping {}...'.format(os.path.split(out)[1]))
        with zipfile.ZipFile(out) as zipf:
            zipfilelist = [info.filename for info in zipf.infolist()]
            listcheck = [f for f in zipfilelist if os.path.exists(os.path.join(out, f))]
            if len(listcheck) > 0:
                print('Overwriting {}...'.format(', '.join(listcheck)))
            zipf.extractall(os.path.split(out)[0])
        del zipf


# url = "ftp://ftp2.census.gov/geo/tiger/TIGER2016/ROADS/tl_2016_04013_roads.zip"
# out = "D:\Mathis\ICSL\stormwater\data\TIGER2016/tl_2016_04013_roads.zip"

def downloadroads(countyfipslist, year=None, outdir=None):
    if year is None:
        year=2018
    if year < 2008:
        raise ValueError("Roads data are not currently available via FTP for years prior to 2008.")
    if outdir==None:
        print('Downloading to {}...'.format(os.getcwd()))
        outdir = os.getcwd()
    if not os.path.exists(outdir):
        print('Creating {}...'.format(outdir))
        os.mkdir(outdir)

    #Open ftp connection
    try:
        rooturl = "ftp://ftp2.census.gov/geo/tiger/TIGER{0}/ROADS".format(year)
        urlp = urlparse.urlparse(rooturl)
        ftp = ftplib.FTP(urlp.netloc)
        ftp.login()
        ftp.cwd(urlp.path)

        #Iterate over county fips
        x=0
        N = len(countyfipslist)

        for county_code in countyfipslist:
            outfile = "tl_{0}_{1}_roads.zip".format(year, county_code)
            out = os.path.join(outdir, outfile)
            if not os.path.exists(out):
                # ftp download
                print "downloading " + outfile
                with open(out, 'w') as fobj:
                    ftp.retrbinary('RETR {}'.format(outfile), fobj.write)

                ######DID NOT WORK
                # urllib.urlretrieve(url, out)
                ######KEPT HANGING
                # try:
                #     r = urllib2.urlopen(url)
                #     print "downloading " + url
                #     with open(out, 'wb') as f:
                #         shutil.copyfileobj(r, f)
                # finally:
                #     if r:
                #         r.close()

                ######KEPT HANGING
                # with contextlib.closing(urllib2.urlopen(url)) as ftprequest:
                #     print "downloading " + url
                #     with open(out, 'wb') as local_file:
                #         shutil.copyfileobj(ftprequest, local_file)

            else:
                print('{} already exists... skipping'.format(outfile))
            x += 1
            print('{}% of county-level data downloaded'.format(100 * round(x / N, 2)))
    finally:
        if ftp:
            ftp.quit()

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
if not arcpy.Exists(hm20tab):
    dlfile(hm20urltab, hm20tab)
else:
    print('{} already exists...'.format(hm20tab))

# Functional system travel 2017 - Annual vehicle-miles, Table VM-2: https://www.fhwa.dot.gov/policyinformation/statistics/2017/vm2.cfm
vm2url = "https://www.fhwa.dot.gov/policyinformation/statistics/2017/vm2.cfm"
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
#Download USHPMS data
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
                if row[1] in [0, 99998, 99999]: #If no urban_code (aka is a rural road or from small urban section)
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'RURAL') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == 7)]['AADT'][0]
                else:
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'URBAN') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == 7)]['AADT'][0]
            else: #Otherwise, assign corresponding f-system aadt
                if row[1] in [0, 99998, 99999]: #If no urban_code (aka is a rural road or from small urban section)
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'RURAL') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == row[4])]['AADT'][0]
                else:
                    row[2] = USAADTformat[(USAADTformat['urban'] == 'URBAN') &
                                          (USAADTformat['state_code'] == row[3]) &
                                          (USAADTformat['F_SYSTEM'] == row[4])]['AADT'][0]
            cursor.updateRow(row)

########################################################################################################################
# FILL IN LOCAL ROADS FOR STATE WHICH DID NOT SUBMIT TO THE HPMS
########################################################################################################################
#Identify states that did not include local roads
fsystemtab = os.path.join(usdotgdb, 'hpms_fsystemstats')
if not arcpy.Exists(fsystemtab):
    arcpy.Statistics_analysis(hpms, statistics_fields=[['f_system','COUNT']], case_field=['state_code', 'f_system'])
fsystemstats = pd.DataFrame(arcpy.da.TableToNumPyArray(fsystemtab, field_names=('state_code','f_system','frequency')))
fsystemstats['proportion'] = fsystemstats['frequency']/\
                             fsystemstats['frequency'].groupby(fsystemstats['state_code']).transform('sum')
missingstatefips = list(fsystemstats[(fsystemstats['f_system']==7) & (fsystemstats['proportion']<0.1)]['state_code'])

#Download list of all FIPS codes in the US (see https://www.census.gov/geo/reference/codes/cou.html for metadata)
dlfile(url = 'https://www2.census.gov/geo/docs/reference/codes/files/national_county.txt', outpath = tiger16dir,
       fieldnames = ['STATE', 'STATEFP', 'COUNTYFP', 'COUNTYNAME', 'CLASSFP'])
fipslist = pd.read_csv(os.path.join(tiger16dir, 'national_county.txt'))

missingcountyfips = fipslist[fipslist['STATEFP'].isin(missingstatefips)]
missingcountyfips_list = list(missingcountyfips['STATEFP'].astype(str).str.zfill(2) +
                              missingcountyfips['COUNTYFP'].astype(str).str.zfill(3))

#Download shapefile data for all counties with missing local roads in the US
downloadroads(countyfipslist=missingcountyfips_list, year=2016, outdir=tiger16dir)

#Merge all these layers

#Subselect to only keep MTFCC IN ('S1400','S1640) https://www2.census.gov/geo/pdfs/reference/mtfccs2018.pdf

#Intersect with urbanized area boundaries to add urban code

#Add aadt_filled field

#Add state_code field

#Add f_system field

#Assign mean AADT from USAADT format

#Merge with hpms

########################################################################################################################
# COMPUTE FUNCTIONAL-CLASS BASED SPEED LIMIT AVERAGE FOR EVERY STATE
########################################################################################################################
#Make sure to replace 0s by NULLs
arcpy.MakeFeatureLayer_management(hpms, 'hpmslyr')
query = "speed_limit > 0"
arcpy.SelectLayerByAttribute_management('hpmslyr', 'NEW_SELECTION', where_clause=query)
arcpy.GetCount_management('hpmslyr')
arcpy.Statistics_analysis(hpms, hmps_spdl, statistics_fields=[['Speed_Limit','MEAN'],['Speed_Limit','MEDIAN']],
                          case_field=['state_code', 'f_system', 'urban_code'])


########################################################################################################################
# PREPARE DATA ON ROAD GRADIENTS
########################################################################################################################
hpms_sub = os.path.join(usdotgdb, 'hpms_subset')
hpms_dens = os.path.join(usdotgdb, 'hpms_dens')
arcpy.CopyFeatures_management(hpms, hpms_dens)
#Densify roads
arcpy.Describe(hpms_sub).SpatialReference.linearUnitName #Check that crs is in meters
arcpy.Densify_edit(hpms_sub, densification_method='DISTANCE', distance='10', max_deviation='1.5')
#Split at vertices
hpms_split = os.path.join(usdotgdb, 'hpms_split')
arcpy.SplitLine_management(hpms_sub,hpms_split)
#Zonal statistics by OBJECTID
arcpy.PolylineToRaster_conversion(hpms_split, 'OBJECTID', hpms_sub + '19', cell_assignment='MAXIMUM_COMBINED_LENGTH',
                                  priority_field= 'aadt_filled', cellsize = NED19proj)
#Rejoin to main
#Average by streetID

# #Run on NED19 max-min for line itself
# arcpy.PolylineToRaster_conversion(PSOSM_allproj, 'osm_id', PSOSMras + '19', cell_assignment='MAXIMUM_COMBINED_LENGTH',
#                                   priority_field= 'fclassADT', cellsize = NED19proj)
# ZonalStatisticsAsTable(PSOSMras + '19', 'osm_id', NED19proj, out_table = rangetab19, statistics_type= 'RANGE', ignore_nodata='NODATA')
#
# #Run on NED13 max-min for line itself
# arcpy.PolylineToRaster_conversion(PSOSM_allproj, 'osm_id', PSOSMras + '13', cell_assignment='MAXIMUM_COMBINED_LENGTH',
#                                   priority_field= 'fclassADT', cellsize = NED13proj)
# ZonalStatisticsAsTable(PSOSMras + '13', 'osm_id', NED13proj, out_table = rangetab13, statistics_type= 'RANGE', ignore_nodata='NODATA')
#
# #Run on smoothed NED 19 max-min for line itself
# ned19_smooth = FocalStatistics(NED19proj, NbrRectangle(3, 3, 'cells'), statistics_type='MEDIAN', ignore_nodata='DATA')
# ned19_smooth.save(NED19smooth)
# ZonalStatisticsAsTable(PSOSMras + '19', 'osm_id', NED19smooth, out_table = rangetab19_smooth, statistics_type= 'RANGE', ignore_nodata='NODATA')
#
# #Run on smoothed NED 13 max-min for line itself
# ned13_smooth = FocalStatistics(NED13proj, NbrRectangle(3, 3, 'cells'), statistics_type='MEDIAN', ignore_nodata='DATA')
# ned13_smooth.save(NED13smooth)
# ZonalStatisticsAsTable(PSOSMras + '13', 'osm_id', NED13smooth, out_table = rangetab13_smooth, statistics_type= 'RANGE', ignore_nodata='NODATA')
#
# #Join all data to road vector
# arcpy.MakeFeatureLayer_management(PSOSM_allproj, 'osmroads')
# arcpy.AddJoin_management('osmroads', 'osm_id', rangetab19, 'osm_id')
# arcpy.AddJoin_management('osmroads', 'osm_id', rangetab13, 'osm_id')
# arcpy.AddJoin_management('osmroads', 'osm_id', rangetab19_smooth, 'osm_id')
# arcpy.AddJoin_management('osmroads', 'osm_id', rangetab13_smooth, 'osm_id')
# arcpy.CopyFeatures_management('osmroads', PSOSM_elv) #Often stays stuck in Python
#
# arcpy.AlterField_management(PSOSM_elv, 'RANGE', 'RANGE19', 'RANGE19')
# arcpy.AlterField_management(PSOSM_elv, 'RANGE_1', 'RANGE13', 'RANGE13')
# arcpy.AlterField_management(PSOSM_elv, 'RANGE_12', 'RANGE19smooth', 'RANGE19smooth')
# arcpy.AlterField_management(PSOSM_elv, 'RANGE_12_13', 'RANGE13smooth', 'RANGE13smooth')
#
# arcpy.AddField_management(PSOSM_elv, 'gradient19', 'FLOAT')
# arcpy.AddField_management(PSOSM_elv, 'gradient13', 'FLOAT')
# arcpy.AddField_management(PSOSM_elv, 'gradient19_smooth', 'FLOAT')
# arcpy.AddField_management(PSOSM_elv, 'gradient13_smooth', 'FLOAT')
# arcpy.CalculateField_management(PSOSM_elv, 'gradient19', '!RANGE19!/!Shape_Length!', 'PYTHON')
# arcpy.CalculateField_management(PSOSM_elv, 'gradient13', '!RANGE13!/!Shape_Length!', 'PYTHON')
# arcpy.CalculateField_management(PSOSM_elv, 'gradient19_smooth', '!RANGE19smooth!/!Shape_Length!', 'PYTHON')
# arcpy.CalculateField_management(PSOSM_elv, 'gradient13_smooth', '!RANGE13smooth!/!Shape_Length!', 'PYTHON')
#
# #Compare to Seattle roads slope values. Apply same method to Seattle road dataset
# roadproj = arcpy.Project_management(roads, os.path.join(PSgdb, 'Seattle_roadproj'), UTM10)
# arcpy.PolylineToRaster_conversion(roadproj, 'OBJECTID', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
#                                   priority_field= 'SURFACEWID', cellsize = NED19proj)
# ZonalStatisticsAsTable(sroadsras, 'Value', NED19proj, out_table = srangetab19,
#                        statistics_type= 'RANGE', ignore_nodata='NODATA')
# ZonalStatisticsAsTable(sroadsras, 'Value', NED19smooth, out_table = srangetab19 + '_smooth',
#                        statistics_type= 'RANGE', ignore_nodata='NODATA')
#
# arcpy.PolylineToRaster_conversion(roadproj, 'OBJECTID', sroadsras, cell_assignment='MAXIMUM_COMBINED_LENGTH',
#                                   priority_field= 'SURFACEWID', cellsize = NED13proj)
# ZonalStatisticsAsTable(sroadsras, 'Value', NED13proj, out_table = srangetab13,
#                        statistics_type= 'RANGE', ignore_nodata='NODATA')
# ZonalStatisticsAsTable(sroadsras, 'Value', NED13smooth, out_table = srangetab13 + '_smooth',
#                        statistics_type= 'RANGE', ignore_nodata='NODATA')
#
# # Fill in and adjust values for those roads outside of NED 1/9 extent
# arcpy.AddField_management(PSOSM_elv, 'gradient_composite', 'FLOAT')
# with arcpy.da.UpdateCursor(PSOSM_elv, ['gradient_composite','gradient19_smooth', 'gradient13_smooth']) as cursor:
#     for row in cursor:
#         if row[1] is not None:
#             row[0] = min(row[1], 0.5)
#         elif row[2] is not None:
#             row[0] = min(row[2], 0.5)
#         else:
#             pass
#         cursor.updateRow(row)

########################################################################################################################
# PREPARE TRANSIT DATA TO CREATE HEATMAP BASED ON BUS ROUTES
########################################################################################################################
#Donwload National Transit Map data (from http://osav-usdot.opendata.arcgis.com/, type NTM)
NTMurls = {'trips': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_7.csv",
           'locations':"http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_11.csv",
           'shapes': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_1.csv",
           'agency': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_6.csv",
           'transfers': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_3.csv",
           'stops': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_0.csv",
           'calendar': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_9.csv",
           'frequencies': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_8.csv",
           'calendar dates': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_2.csv",
           'fare rules': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_10.csv",
           'feed info': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_5.csv",
           'routes': "http://osav-usdot.opendata.arcgis.com/datasets/816b89974fad4fcf9c7356a3b704ffd1_11.csv"}
for tab in NTMurls:
    print(tab)
    try:
        dlfile(url = NTMurls[tab], outpath = NTMdir)
    except Exception as e:
        print e
        pass

#Download shapes as opendata website doesn't work
baseURL = "https://geo.dot.gov/server/rest/services/NTAD/GTFS_NTM/MapServer/1/query"
arcpy.env.workspace = NTMdir
basename='NTM_shapes'
itersize=1000
IDrange = range(0,100000, itersize)
for i, j in itertools.izip(IDrange, IDrange[1:]):
    itertry = itersize
    downlist = []
    # It seems like REST API server also has a limitation on the size of the downloads so sometimes won't allow
    # Therefore, if fails to download, try smaller increments until reaches increments of 2 if still fails at increments
    # of 2, then throw a proper error and break
    while True:
        try:
            IDrangetry = list(
                sorted(set(range(i, j + 1, itertry) + [j])))  # To make sure that the list goes until the maximum
            # Loop with smaller increment within range
            for k, l in itertools.izip(IDrangetry, IDrangetry[1:]):
                print('From {0} to {1}'.format(k, l))
                where = "OBJECTID>={0} AND OBJECTID<{1}".format(k, l)
                # &geometryType=esriGeometryPoint
                query = "?where={}&returnGeometry=true&f=json&outFields=*".format(where)
                fsURL = baseURL + query
                fs = arcpy.FeatureSet()
                fs.load(fsURL)
                if long(arcpy.GetCount_management(fs)[0]) > 0:
                    outname = '{0}_{1}_{2}'.format(basename, k, l)
                    downlist.append(outname)
                    arcpy.CopyFeatures_management(fs, outname)
                    print(outname)
                else:
                    print('No data from OBJECTID {0} to OBJECTID {1}'.format(k, l))
            break
        except:
            if itertry > 5:
                print('Count not download, delete previous {0} datasets, try downloading in smaller increments'.format(
                    len(downlist)))
                if len(downlist) > 0:
                    for fc in downlist:
                        arcpy.Delete_management(fc)
                    downlist = []
                itertry = itertry / 2
            else:
                e = sys.exc_info()[1]
                print('Exit with error: ' + e.args[0])
                # sys.exit(1)
                break
arcpy.Merge_management(arcpy.ListFeatureClasses('{}_*'.format(basename)),output='{}'.format(basename))
for fc in arcpy.ListFeatureClasses('{}_*'.format(basename)):
    arcpy.Delete_management(fc)
    print('Deleted {}...'.format(fc))

#Format data
GTFStoSHPweeklynumber(gtfs_dir= NTMdir, out_gdb=os.path.join(rootdir, 'results/NTM.gdb'), out_fc = 'NTM',
                      keep = True)

# #Only keep transports overground (although some subways might have some overground sections)
# arcpy.MakeFeatureLayer_management(PStransit, 'PStransit_lyr',
#                                   where_clause= '(route_type = 3) AND (MIN_service_len > 1) AND (SUM_adjustnum > 0)')
# arcpy.CopyFeatures_management('PStransit_lyr', PStransitbus)
# arcpy.Project_management(PStransitbus, PStransitbus_proj, UTM10)
#
# #Create raster of weekly number of buses at the same resolution as bing data
# # Convert weekly number of buses to integer
# arcpy.AddField_management(PStransitbus_proj, 'adjustnum_int', 'SHORT')
# arcpy.CalculateField_management(PStransitbus_proj, 'adjustnum_int',
#                                 expression='int(10*!SUM_adjustnum!+0.5)', expression_type='PYTHON')
#
# #Split lines at all intersections so that small identical overlapping segments can be dissolved
# arcpy.SplitLine_management(PStransitbus_proj, PStransitbus_proj + '_split') #Split at intersection
# arcpy.FindIdentical_management(PStransitbus_proj + '_split', "explFindID", "Shape") #Find overlapping segments and make them part of a group (FEAT_SEQ)
# arcpy.MakeFeatureLayer_management(PStransitbus_proj + '_split', "intlyr")
# arcpy.AddJoin_management("intlyr", arcpy.Describe("intlyr").OIDfieldName, "explFindID", "IN_FID", "KEEP_ALL")
# arcpy.Dissolve_management("intlyr", PStransitbus_splitdiss, dissolve_field='explFindID.FEAT_SEQ',
#                           statistics_fields=[[os.path.split(PStransitbus_proj)[1] + '_split.adjustnum_int', 'SUM']]) #Dissolve overlapping segments
# arcpy.RepairGeometry_management(PStransitbus_splitdiss, delete_null = 'DELETE_NULL') #sometimes creates empty geom
# #Get the length of a half pixel diagonal to create buffers for
# #guaranteeing that segments potentially falling within the same pixel are rasterized separately
# tolerance = (2.0**0.5)*float(res.getOutput(0))/2
# ExplodeOverlappingLines(PStransitbus_splitdiss, tolerance)
#
# #For each set of non-overlapping lines, create its own raster
# tilef = 'expl'
# tilelist = list(set([row[0] for row in arcpy.da.SearchCursor(PStransitbus_splitdiss, [tilef])]))
# outras_base = os.path.join(rootdir, 'results/transit.gdb/busnum_')
# arcpy.env.snapRaster = template_ras
# for tile in tilelist:
#     outras = outras_base + str(tile)
#     if not arcpy.Exists(outras):
#         selexpr = '{0} = {1}'.format(tilef, tile)
#         print(selexpr)
#         arcpy.MakeFeatureLayer_management(PStransitbus_splitdiss, 'bus_lyr', where_clause= selexpr)
#         arcpy.PolylineToRaster_conversion('bus_lyr', value_field='adjustnum_int', out_rasterdataset=outras, cellsize=res)
#
# #Mosaic to new raster
# arcpy.env.workspace = os.path.split(outras_base)[0]
# transitras_tiles = arcpy.ListRasters('busnum_*')
# arcpy.MosaicToNewRaster_management(transitras_tiles, arcpy.env.workspace, os.path.split(PStransitras)[1],
#                                    pixel_type='32_BIT_UNSIGNED', number_of_bands= 1, mosaic_method = 'SUM')
# for tile in transitras_tiles:
#     print('Deleting {}...'.format(tile))
#     arcpy.Delete_management(tile)
# arcpy.ClearEnvironment('Workspace')




########################################################################################################################
# PREPARE LAND USE DATA
########################################################################################################################

