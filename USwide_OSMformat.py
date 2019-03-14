import arcpy
import pandas as pd
import geopandas as gpd
import os
import re
import itertools
import shapely
from BeautifulSoup import BeautifulSoup
import urllib2

rootdir = "D:\Mathis\ICSL\stormwater"
USDOTdir = os.path.join(rootdir, "data\USDOT_0319")
resdir = os.path.join(rootdir, 'results/usdot')

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
#Get HPMS data for every state
"https://www.fhwa.dot.gov/policyinformation/hpms/shapefiles.cfm"


#Merge it all
#Merge by spatial overlap
#Get VM-2 estimates of Vehicle Miles Traveled for each state and average it for small roads at state level (idea fromhttps://www.pnas.org/content/pnas/suppl/2015/04/01/1421723112.DCSupplemental/pnas.1421723112.sapp.pdf)
#Highway Statistics Book 2017
# Public road lenth 2017 - miles by functional system, Table HM-20: https://www.fhwa.dot.gov/policyinformation/statistics/2017/hm20.cfm
# Functional system travel 2017 - Annual vehicle-miles, Table VM-2: https://www.fhwa.dot.gov/policyinformation/statistics/2017/vm2.cfm