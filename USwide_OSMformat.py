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

pd.options.display.max_columns = 20

rootdir = "D:\Mathis\ICSL\stormwater"
USDOTdir = os.path.join(rootdir, "data\USDOT_0319")
resdir = os.path.join(rootdir, 'results/usdot')

#Output variables
usdotgdb = os.path.join(resdir, 'usdot.gdb')
ushpms = os.path.join(usdotgdb, 'ushpms')

if not arcpy.Exists(usdotgdb):
    arcpy.CreateFileGDB_management(resdir, out_name = 'usdot')

########################################################################################################################
# FORMAT TMAS DATA
########################################################################################################################
exampletab = os.path.join(USDOTdir, "dec_2017_ccs_data\\12. DEC\\AK_DEC_2017 (TMAS).VOL")

STAlist = list(itertools.chain.from_iterable( #To unnest list
    [filter(re.compile('.*[.]STA$').search, [os.path.join(dirpath, file) for file in filenames])
     for (dirpath, dirnames, filenames) in os.walk(USDOTdir)]))

#Import tabs
dectab = list(itertools.chain.from_iterable( #To unnest list
    [filter(re.compile('.*[.]VOL$').search, [os.path.join(dirpath, file) for file in filenames])
     for (dirpath, dirnames, filenames) in os.walk(USDOTdir)]))

#Set crs
wgs84 = {'init' :'epsg:4326'}

#Format .STA files
stafields = {
    'type':         [0,1],      #Record Type
    'FIPS_state':   [1, 3],     #FIPS State Code
    'stationID':    [3, 9],     #Station ID
    'direction':    [9, 10],    #Direction of Travel Code
    'lane':         [10, 11],   #Lane of travel
    'year':         [11, 15],   #Year of data
    'fclass':       [15, 17],   #Functional classification code
    'laneN':        [17, 18],   #Number of lanes in direction indicated
    'TMAStype':     [18, 19],   #Sample type for TMAS
    'laneNvol':     [19, 20],   #Number of lanes monitored for traffic volume
    'methodvol':    [20, 21],   #Method of traffic volume counting
    'laneNclass':   [21, 22],   #Number of lanes monitored for vehicle class
    'mechaclass':   [22, 23],   #Mechanism of vehicle classification
    'methodclass':  [23, 24],   #Method for vehicle classification
    'classgroups':  [24, 26],   #Vehicle classification groupings
    'laneNweight':  [26, 27],   #Number of lanes monitored for truck weight
    'methodweight': [27, 28],   #Method of truck weighing
    'calibweight':  [28, 29],   #Calibration of weighing system
    'methoddata':   [29, 30],   #Method of data retrieval
    'sensortype':   [30, 31],   #Type of sensor
    'sensortype2':  [31, 32],   #Second type of sensor
    'purpose':      [32, 33],   #Primary purpose
    'LRSID':        [33, 93],   #LRS Route ID
    'LRSloc':       [93, 101],  #LRS location point
    'Latitude':     [101, 109], #Latitude
    'Longitude':    [109, 118], #Longitude
    'LTTPID':       [118, 122], #LTPP Site identification
    'stationIDpre': [122, 128], #Previous station ID
    'yearstart':    [128, 132], #Year station established
    'yearend':      [132, 136], #Year station discontinued
    'FIPS_county':  [136, 139], #FIPS county code
    'HPMStype':     [139, 140], #HPMS sample type
    'HPMSID':       [140, 152], #HPMS sample identifier
    'NHS':          [152, 153], #National Highway System
    'routesign':    [153, 155], #Posted route signing
    'routenumber':  [155, 163], #Posted signed routed number
    'stationloc':   [163, 213], #Station location
}

#Define old field system based on 2001_TMG_Chapter_6
staTMASfields = {
    'type':         [0,1],      #Record Type
    'FIPS_state':   [1, 3],     #FIPS State Code
    'stationID':    [3, 9],     #Station ID
    'direction':    [9, 10],    #Direction of Travel Code
    'lane':         [10, 11],   #Lane of travel
    'year':         [11, 13],   #Year of data
    'fclass':       [13, 15],   #Functional classification code
    'laneN':        [15, 16],   #Number of lanes in direction indicated
    'voltype':      [16, 17],   #Sample type for traffic volume
    'laneNvol':     [17, 18],   #Number of lanes monitored for traffic volume
    'methodvol':    [18, 19],   #Method of traffic volume counting
    'classtype':    [19, 20],   #Sample type for vehicle classification
    'laneNclass':   [20, 21],   #Number of lanes monitored for vehicle class
    'methodclass':  [21, 22],   #Method for vehicle classification
    'algoclass':    [22, 23],   #Algorithm of vehicle classification
    'classgroups':  [23, 25],   #Classification System for Vehicle Classification
    'typeweight':   [25, 26],   #Sample type for truck weight
    'laneNweight':  [26, 27],   #Number of lanes monitored for truck weight
    'methodweight': [27, 28],   #Method of truck weighing
    'calibweight':  [28, 29],   #Calibration of weighing system
    'methoddata':   [29, 30],   #Method of data retrieval
    'sensortype':   [30, 31],   #Type of sensor
    'sensortype2':  [31, 32],   #Second type of sensor
    'purpose':      [32, 33],   #Primary purpose
    'LRSID':        [33, 45],   #LRS Route ID
    'LRSloc':       [45, 51],  #LRS location point
    'Latitude':     [51, 59], #Latitude
    'Longitude':    [59, 68], #Longitude
    'SHRPID':       [68, 72], #SHRP Site identification
    'stationIDpre': [72, 78], #Previous station ID
    'yearstart':    [78, 80], #Year station established
    'yearend':      [80, 82], #Year station discontinued
    'FIPS_county':  [82, 85], #FIPS county code
    'HPMStype':     [85, 86], #HPMS sample type
    'HPMSID':       [86, 98], #HPMS sample identifier
    'NHS':          [87, 99], #National Highway System
    'routesign':    [99, 100], #Posted route signing
    'routenumber':  [100, 108], #Posted signed routed number
    'stationloc':   [108, 167], #Station location
}

#Format .VOL files
volfields = {
    'type':         [0,1],      #Record type
    'FIPS_state':   [1,3],      #FIPS State Code
    'fclass':       [3,5],      #Functional classification
    'stationID':    [5,11],     #Station identification
    'direction':    [11,12],    #Direction of travel
    'lane':         [12, 13],   #Lane of travel
    'year':         [13, 17],   #Year of data
    'month':        [17, 19],   #Month of data
    'day':          [19, 21],   #Day of data
    'dayweek':      [21, 22],   #Day of week
    'count00_01':   [22, 27],   #Traffic volume counted, after 00:00- to 01:00
    'count01_02':   [27, 32],
    'count02_03':   [32, 37],
    'count03_04':   [37, 42],
    'count04_05':   [42, 47],
    'count05_06':   [47, 52],
    'count06_07':   [52, 57],
    'count07_08':   [57, 62],
    'count08_09':   [62, 67],
    'count09_10':   [67, 72],
    'count10_11':   [72, 77],
    'count11_12':   [77, 82],
    'count12_13':   [82, 87],
    'count13_14':   [87, 92],
    'count14_15':   [92, 97],
    'count15_16':   [97, 102],
    'count16_17':   [102, 107],
    'count17_18':   [107, 112],
    'count18_19':   [112, 117],
    'count19_20':   [117, 122],
    'count20_21':   [122, 127],
    'count21_22':   [127, 132],
    'count22_23':   [132, 137],
    'count23_24':   [137, 142],
    'restrictions': [142, 143]  #Restrictions
}

#Format station info table
sta2017 = {}
for tabpath in STAlist:
    print(tabpath)
    sta2017[tabpath] = pd.read_table(tabpath, sep='\n', header=None)
    sta2017[tabpath] = sta2017[tabpath].rename({0: 'rawfield'}, axis='columns')
    for f in staTMASfields:
        sta2017[tabpath][f] = sta2017[tabpath]['rawfield'].str[staTMASfields[f][0]:staTMASfields[f][1]]

#Concatenate pds
sta2017_all = pd.concat(sta2017)
sta2017_all.columns

#Format coordinates (make up for weird entries)
sta2017_all['Longitude'] = sta2017_all.Longitude.astype(float)
sta2017_all['Latitude'] = sta2017_all.Latitude.astype(float)


sta2017_all.loc[sta2017_all.Latitude < 15000000, 'Latitude'] = \
    sta2017_all.loc[sta2017_all.Latitude < 15000000, 'Latitude']/10**5
sta2017_all.loc[sta2017_all.Latitude >= 15000000, 'Latitude'] = \
    sta2017_all.loc[sta2017_all.Latitude >= 15000000, 'Latitude']/10**6

sta2017_all.loc[sta2017_all.Longitude < 18000000, 'Longitude'] = \
    sta2017_all.loc[sta2017_all.Longitude < 18000000, 'Longitude']/(-10**5)
sta2017_all.loc[(18000000 < sta2017_all.Longitude) & (sta2017_all.Longitude <= 180000000), 'Longitude'] = \
    sta2017_all.loc[(18000000 < sta2017_all.Longitude) & (sta2017_all.Longitude <= 180000000), 'Longitude']/(-10**6)
sta2017_all.loc[sta2017_all.Longitude > 180000000, 'Longitude'] = \
    sta2017_all.loc[sta2017_all.Longitude > 180000000, 'Longitude']/(-10**7)

#Write out geodataframe
sta2017_all['Coordinates'] = list(zip(sta2017_all.Longitude, sta2017_all.Latitude))
sta2017_all['Coordinates'] = sta2017_all['Coordinates'].apply(shapely.geometry.Point)
sta2017_gdf = gpd.GeoDataFrame(sta2017_all, geometry='Coordinates')
sta2017_gdf.crs = wgs84
sta2017_gdf.to_file(driver = 'ESRI Shapefile', filename= os.path.join(resdir, "sta2017.shp"))

#Format december data table
dec2017 = {}
for tabpath in dectab:
    print(tabpath)
    dec2017[tabpath] = pd.read_table(tabpath, sep='\n', header=None)
    dec2017[tabpath] = dec2017[tabpath].rename({0: 'rawfield'}, axis='columns')
    for f in volfields:
        dec2017[tabpath][f] = dec2017[tabpath]['rawfield'].str[volfields[f][0]:volfields[f][1]]
dec2017_all = pd.concat(dec2017)
dec2017_all.columns

len(list(set(dec2017_all['stationID'])))

########################################################################################################################
# FORMAT HIGHWAY AADT DATA
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

#Get HPMS data for every state
hpmsurl = "https://www.fhwa.dot.gov/policyinformation/hpms/shapefiles.cfm"
hpms_page = urllib2.urlopen(hpmsurl)
hpms_soup = BeautifulSoup(hpms_page)
hpmsshp = []
hpms_regex = re.compile(".*hpms[/]shapefiles[/](?!nationalarterial)[a-z]+2017[.]zip")
for link in hpms_soup.findAll('a', attrs={'href': hpms_regex}):
    zipurl = os.path.split(hpmsurl)[0] + '/' + \
             '/'.join([i for i in link.get('href').split('/') if i not in os.path.split(hpmsurl)[0].split('/')])
    outpath = os.path.join(USDOTdir, 'HPMS', os.path.split(zipurl)[1])

    if outpath not in hpmsshp:
        hpmsshp.append('{}.shp'.format(os.path.splitext(outpath)[0]))

    if not os.path.exists(outpath):
        dlfile(zipurl, outpath)
    else:
        print('{} already exists...'.format(outpath))

#Merge all HPMS data
print('{} missing shapefiles...'.format(len([file for file in hpmsshp if not os.path.exists(file)])))
arcpy.Merge_management([file for file in hpmsshp if os.path.exists(file)], ushpms)

#Merge with OSM by spatial overlap

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