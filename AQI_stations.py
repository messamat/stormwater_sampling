import arcpy
import mendeleev
import pandas as pd
import itertools
from arcpy.sa import *
import os
import re
import geopandas as gpd
from shapely.geometry import Point
import us
import numpy as np
import netCDF4
from datetime import datetime
import timeit
import glob

#Custom functions
from Download_gist import *

#Set options
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'
AQIdir = os.path.join(rootdir, 'data/EPA_AirData_201902')
NARRdir = os.path.join(rootdir, 'data/NARR_201904')
if not os.path.isdir(NARRdir):
    os.mkdir(NARRdir)

#Import variables
monitors = pd.read_csv(os.path.join(AQIdir, 'aqs_monitors/aqs_monitors.csv'))
sites = pd.read_csv(os.path.join(AQIdir, 'aqs_sites/aqs_sites.csv'))

#Output variables
outdir = os.path.join(rootdir, 'results/airdata')
if not os.path.isdir(outdir):
    os.mkdir(outdir)
sites_out = os.path.join(outdir, 'airsites.shp')
sites_outbuf = os.path.join(outdir, 'airsites_550buf.shp')
sites_outbufdis = os.path.join(outdir, 'airsites_550bufdis.shp')
sites_outbufunion = os.path.join(outdir, 'airsites_550bufunion.shp')

AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')
#Create gdb for analysis
if arcpy.Exists(AQIgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results/airdata'), 'AQI.gdb')

#Functions
def getclosest_ij(lats, lons, latpt, lonpt):
    # find squared distance of every point on grid
    dist_sq = (lats - latpt) ** 2 + (lons - lonpt) ** 2
    # 1D index of minimum dist_sq element
    minindex_flattened = dist_sq.argmin()
    # Get 2D index for latvals and lonvals arrays from 1D index
    return np.unravel_index(minindex_flattened, lats.shape)

def getclosest_ij_df(lats, lons, latlonpt):
    """# A function to find the index of the point closest pt(in squared distance) to give lat/lon value.
    latlonpt may be a 2-column dataframe (lat-lon) or a list/tuple of lat=lon for a point.
    If data.frame, returns a data.frame of unique coordinates from the original data.frame together with index columns"""

    if isinstance(latlonpt, pd.DataFrame):
        print('Latitude and longitude to locate are panda series, proceeding in pd framework...')
        if latlonpt.shape[1] == 2:
            #Get unique set of coordinates
            latlonpt_unique = latlonpt.drop_duplicates()
            #Get closest grid for each unique coordinate and return data.frame to join to original one
            return latlonpt_unique.apply(
                lambda x: pd.Series([x[0], x[1]] + list(getclosest_ij(lats, lons, float(x[0]), float(x[1]))),
                                    index = list(latlonpt.columns.values) + ['iy_min', 'ix_min']), axis=1)
        else:
            raise ValueError('Wrong number of columns in latlonpt')

def extractCDFtoDF(dfcsv, pattern, indir, varname, level=None, outcsv=None):
    """Warning: will edit source df with extra columns"""
    pathpattern = os.path.join(indir, pattern)
    cdflist = glob.glob(pathpattern)
    print('Extracting {0}'.format(cdflist))

    #Get netCDF data
    try:
        sourcef = netCDF4.MFDataset(cdflist) #Source it as an MFD, otherwise
    except Exception as e:
        traceback.print_exc()
        if isinstance(e, RuntimeError):
            for fpath in glob.glob(cdflist):
                try:
                    netCDF4.Dataset(fpath)
                except:
                    print('Error stems from {}...'.format(fpath))
                    pass

    #Print dimensions and their length
    print('Dimensions:')
    for dim in sourcef.dimensions:
        dimobj = sourcef.dimensions[dim]
        print('{0}, length:{1}, unlimited:{2}'.format(dim, len(dimobj), dimobj.isunlimited()))

    #Extract dimensional variables
    timev = sourcef.variables['time']
    latv, lonv = sourcef.variables['lat'], sourcef.variables['lon']

    #If supplied path to csv, read csv as df
    if isinstance(dfcsv, str) and os.path.exists(dfcsv):
        df = pd.read_csv(dfcsv)
    elif isinstance(dfcsv, pd.DataFrame):
        df = dfcsv
    else:
        raise ValueError('Input df argument neither csv nor dataframe ')

    #Make sure that right columns are in df
    if any([c not in dfcsv.columns for c in ['Latitude', 'Longitude', 'date']]):
        raise ValueError('Column missing in df; either Latitude, Longitude or date')

    #Get index of the pixels closest to each point in df
    print('Getting index of pixels closest to points in df...')
    df_ij = df.merge(getclosest_ij_df(latv[:], lonv[:], df[['Latitude', 'Longitude']]),
                     how='left', on=['Latitude', 'Longitude'])
    #Convert dates in df to index
    print('Converting dates to index...')
    df_ij['date_index'] = \
        df_ij['date'].apply(lambda x: np.where(timev[:] == netCDF4.date2num(datetime.strptime(x, '%Y-%m-%d'),
                                                                            timev.units))[0][0])
    #Extract variables
    print('Extracting variable...')
    if level is not None:
        if 'level' in sourcef.dimensions:
            level_i = np.where(sourcef.variables['level'][:] == level)[0][0]
            varname = '{0}{1}'.format(varname, level)
            df_ij[varname] = \
                df_ij.apply(lambda x: sourcef[x.date_index, level_i, int(x['iy_min']), int(x['ix_min'])], axis=1)
        else:
            raise ValueError("A 'level' argument was provided but netCDF does not have a level dimension")
    else:
        df_ij[varname] = \
            df_ij.apply(lambda x: sourcef[x.date_index, int(x['iy_min']), int(x['ix_min'])], axis=1)

    #Delete intermediate columns
    print('Deleting indices columns')
    df.drop(['iymin','ix_min','date_index'], axis=1)

    #Write df out
    if (isinstance(dfcsv, str) and os.path.exists(dfcsv)) or (outcsv is not None):
        print('Writing df out...')
        df_ij.to_csv(dfcsv if outcsv is None else outcsv)

#-----------------------------------------------------------------------------------------------------------------------
# SELECT SITES THAT RECORD CHEMICAL CONCENTRATIONS
#-----------------------------------------------------------------------------------------------------------------------
#Create list of elements
elems = list(itertools.chain.from_iterable(mendeleev.get_table('elements')[['name', 'symbol']].values))
elems_out = ['monoxide', 'dioxide', 'H', 'C', 'O', 'N', 'S', 'Hydrogen', 'Carbon', 'Oxygen', 'Nitrogen', 'Sulfur',
             '.*ene']
elems_regex = re.compile('|'.join(['\\b{}\\b'.format(e) for e in elems]), re.IGNORECASE)
elemsout_regex = re.compile('|'.join(['\\b{}\\b'.format(e) for e in elems_out]), re.IGNORECASE)
#Create list of parameter codes based on https://aqs.epa.gov/aqsweb/documents/codetables/methods_all.html
######################### RECHECK PARAMETER CODES ###############
paramsel = ['11','12','14','22','65','82','85','86','89']
param_regex = re.compile('|'.join(['(^{}.*)'.format(p) for p in paramsel]))
monitors_chem = monitors.loc[(monitors['Parameter Code'].astype(str).str.contains(param_regex)) &
                             (monitors['Parameter Name'].str.contains(elems_regex)) &
                             ~(monitors['Parameter Name'].str.contains(elemsout_regex)),:]
#pd.unique(monitors_chem['Parameter Name'])

#Pad monitors' codes with 0s to match sites' codes
pd.unique(monitors_chem['County Code'])
monitors_chem['State Code'] = monitors_chem['State Code'].astype(str).str.pad(2, 'left', '0')

#Remove sites outside of the conterminous US (to not exceed 50,000 tile limit)
fips_excluded = [us.states.lookup(abr).fips for abr in ['PR','AK','HI','VI','UM', 'MP', 'GU'] if us.states.lookup(abr)]

#Compute Unique Site Identifier and select sites that record chemical concentrations
monitors_chem['UID'] = monitors_chem['State Code'] + \
                       monitors_chem['County Code'].astype(str) + \
                       monitors_chem['Site Number'].astype(str)
sites['UID'] = sites['State Code'].astype(str).astype(str) + \
               sites['County Code'].astype(str) + \
               sites['Site Number'].astype(str)
sites_chem = sites[(sites['UID'].isin(pd.unique(monitors_chem['UID']))) &
                   ~(sites['State Code'].isin(fips_excluded))]
len(sites_chem)
len(pd.unique(sites_chem['UID'])) #Most sites are unique

#-----------------------------------------------------------------------------------------------------------------------
# CREATE POINT SHAPEFILE AND BUFFER
#-----------------------------------------------------------------------------------------------------------------------
sites_chem.groupby('Datum')['UID'].nunique() #Check what datums data are in

#Create geodataframes
def pdtogpd_datum(df, datum):
    "Subset df based on datum and create Geodataframe with the right assigned datum"
    crs = {'WGS84': {'init': 'epsg:4326'}, 'NAD27': {'init': 'epsg:4267'},'NAD83': {'init': 'epsg:6318'}} #'NAD83': {'init': 'epsg:1116'}}
    subdat = df[df['Datum']==datum]
    geom = [Point(xy) for xy in zip(subdat.Longitude, subdat.Latitude)]
    return(gpd.GeoDataFrame(subdat, crs=crs[datum], geometry=geom))

sites_wgs84 = pdtogpd_datum(sites_chem, 'WGS84')
sites_nad27 = pdtogpd_datum(sites_chem, 'NAD27')
sites_nad83 = pdtogpd_datum(sites_chem, 'NAD83')

sites_chem['Longitude'].describe()

#Project them all to esri projection 102003 - usa contiguous albers equal area conic
#Spatial reference: http://spatialreference.org/ref/esri/usa-contiguous-albers-equal-area-conic/
albers = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
sites_gpd = pd.concat([sites_wgs84.to_crs(albers),sites_nad27.to_crs(albers), sites_nad83.to_crs(albers)])
#Remove outliers with lon=0
sites_gpd['Longitude'].describe()
sites_gpd = sites_gpd[sites_gpd['Longitude'] != 0]

#Create 550 m radius buffers and dissolve them
sites_buf = sites_gpd['geometry'].buffer(distance=550)
sites_bufgpd = gpd.GeoDataFrame(sites_buf, crs=albers).rename(columns={0:'geometry'}).set_geometry('geometry')

#Dissolve all
sites_bufgpd['diss'] = 1
sites_bufdis = sites_bufgpd.dissolve(by='diss')

#Dissolve only overlapping buffers
sites_bufunion = gpd.GeoDataFrame([polygon for polygon in sites_bufgpd.unary_union], crs=albers). \
    rename(columns={0:'geometry'}).set_geometry('geometry')

#Output to shapefile
sites_gpd.to_file(sites_out, driver = 'ESRI Shapefile')
sites_bufgpd.to_file(sites_outbuf, driver = 'ESRI Shapefile')
sites_bufdis.to_file(sites_outbufdis, driver = 'ESRI Shapefile')
sites_bufunion.to_file(sites_outbufunion, driver = 'ESRI Shapefile')

#-----------------------------------------------------------------------------------------------------------------------
# DOWNLOAD DAILY METEOROLOGICAL DATA FOR AQI STATIONS
# Get selected covariates from Porter et al. 2015:
# Investigating the observed sensitivities of air-quality extremes to meteorological drivers via quantile regression
#-----------------------------------------------------------------------------------------------------------------------

# All text in Table 2, second panel of Porter et al. 2015
porter2015_covarlist_raw = """rhum.2m_mean dswrf_mean air.2m_max hpbl_mean
air.2m_max wspd.10m_mean air.sfc_9x9_nightmin.6daymax vwnddir.10m_mean
dswrf_daymin.3daymean dswrf_daymean.diff crain_9x9_nightmax wspd.10m_daymax.3daymax
vwnddir.10m_mean vwnddir.10m_mean wspd.10m_max.3daymean crain_9x9_nightmax
crain_9x9_daymean lts_daymin vwnddir.10m_mean wspd.10m_nightmax
fire shum.2m_min lftx4_mean rhum.2m_mean
tke.hl1_9x9_daymax uwnddir.10m_mean lts_daymin uwnddir.10m_mean
uwnddir.10m_daymean.3daymean crain_9x9_daymax uwnddir.10m_daymean.3daymean wspd.10m_max.3daymin
air.sfc_9x9_daymin.3daymean dswrf_min.3daymin shum.2m_daymean.diff rpi_max
rpi_max fire crain_9x9_max.6daymean uwnddir.10m_nightmean.3daymean
lts_mean air.sfc_9x9_mean.6daymean rpi_max dswrf_daymin.6daymax
dswrf_min.6daymin hpbl_daymax vwnd.500_min lftx4_nightmin.3daymean
vwnd.500_min hcdc_9x9_daymax vwnd.500_daymax.6daymax shum.2m_nightmin.6daymean
hpbl_nightmean.3daymin pres.sfc_nightmin.6daymean pres.sfc_max fire
vvel.hl1_mean.6daymean rpi_nightmax.6daymean hgt.850_max.6daymax"""
porter2015_covarlist = re.split('\s|\n', porter2015_covarlist_raw)
porter2015_covarlist.append(None)

#Convert to data frame with columns
covar_df = pd.DataFrame(
    np.array(porter2015_covarlist). \
        reshape(int(np.ceil(float(len(porter2015_covarlist))/4)),4),
    columns = ['O3summer', 'O3winter', 'PM25summer', 'PM25winter'])

#Select PM2.5 subsetted covariates and separate into main variable and subvariable names
covar_pm25 = [var.split('_') for var in
              list(set(itertools.chain.from_iterable(covar_df[['PM25summer', 'PM25winter']].values.tolist())))
              if var is not None]
#Get variable base name for downloading NARR data
covar_pm25_base = list(set([var[0] for var in covar_pm25]))
covar_pm25_base

#Modify list to only keep raw variables and match NARR abbreviations
covar_sub = [var for var in covar_pm25_base if var not in ['lts', 'fire', 'uwnddir.10m', 'vwnddir.10m', 'wspd.10m', 'rpi']] +\
            ['uwnd.10m', 'vwnd.10m']

#Assign variable to either monolevel or pressure direction in the ftp server (ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/)
varfolder_dic = {}
for var in covar_sub:
    if re.search('.*[.][0-9]{1,4}$', var):
        varfolder_dic[var] = 'pressure'
    else:
        varfolder_dic[var] = 'monolevel'

#Download all NARR data from 2010 to 2018 for each variable
yearlist = range(2010, 2019)
for var in varfolder_dic:
    if varfolder_dic[var] == 'pressure':
        varbase = var.split('.')[0]
        downloadNARR(folder=varfolder_dic[var], variable=varbase,
                     years=['{0}{1}'.format(year,month) for year in yearlist
                            for month in [str(i).zfill(2) for i in range(1,13)]],
                     outdir=NARRdir)
    else:
        downloadNARR(folder=varfolder_dic[var], variable=var, years=yearlist, outdir=NARRdir)

#-----------------------------------------------------------------------------------------------------------------------
# DOWNLOAD SPECIATION DATA FOR ALL MONITORING SITES
#-----------------------------------------------------------------------------------------------------------------------
yearlist = range(2010, 2019)
epadl_url = "https://aqs.epa.gov/aqsweb/airdata/"
spec25_list = [os.path.join(epadl_url, "daily_SPEC_{}.zip".format(year)) for year in yearlist]
spec10_list = [os.path.join(epadl_url, "daily_PM10SPEC_{}.zip".format(year)) for year in yearlist]
for file in spec25_list + spec10_list:
    dlfile(file, outpath=os.path.join(AQIdir))

#Collate all data
airdat_uniquetab = os.path.join(AQIdir, 'daily_SPEC_unique.csv')
if not os.path.exists(airdat_uniquetab):
    airdatall = os.path.join(AQIdir, 'daily_SPEC_collate.csv')
    airdatlist = mergedel(AQIdir, 'daily_.*SPEC.*[.]csv$', airdatall, verbose=True)
    airdat_df = pd.read_csv(airdatall)
    airdat_df.shape
    airdat_df.dtypes

    #Compute unique ID for each monitoring site
    airdat_df['UID'] = airdat_df['State Code'].astype(str).astype(str) + \
                       airdat_df['County Code'].astype(str) + \
                       airdat_df['Site Num'].astype(str)

    #Compute df of unique UID-time records to extract meteorological variables with
    subcols = ['UID', 'Latitude', 'Longitude', 'Date Local']
    airdat_unique = pd.unique(
        airdat_df[subcols].
            apply(lambda x: '._.'.join(x.astype(str)), axis=1)) #Concatenate columns

    airdat_uniquedf = pd.DataFrame(
        [i.split('._.') for i in list(airdat_unique)],
        columns = subcols) #Resplit and format to df
    airdat_uniquedf.to_csv(os.path.join(AQIdir, 'daily_SPEC_unique.csv'))
else:
    airdat_uniquedf = pd.read_csv(airdat_uniquetab)

airdat_uniquedf['date'] = airdat_uniquedf['Date Local']

#-----------------------------------------------------------------------------------------------------------------------
# EXTRACT ALL METEOROLOGICAL VARIABLES FOR EACH STATION-DATE COMBINATION
#-----------------------------------------------------------------------------------------------------------------------
#Read all years of crain data
extractCDFtoDF(airdat_uniquedf, pattern='crain.*.nc', indir=NARRdir,
               varname='crain', level=None, outcsv=os.path.join(rootdir, 'results/daily_SPEC_crain.csv'))




# pattern = 'hgt.*.nc'
# varname = 'hgt'
# df = airdat_uniquedf
# indir = NARRdir
# level = 850



#Check variables
for a in air2mf.variables:
    print(a)
    varobj = air2mf.variables[a]
    print(varobj)

levelv = air2mf.variables['level']
levelv[:]
#Get variables
crainv = crainf.variables['crain']
timev = crainf.variables['time']
lat, lon = crainf.variables['lat'], crainf.variables['lon']
timevals = timev[:]
tarray = netCDF4.num2date(timevals, timev.units) #Get all dates formated in date.time

#Add x and y grid index for each unique air monitoring station record
airdat_uniquedf_ij = airdat_uniquedf.merge(getclosest_ij_df(lat[:],lon[:],airdat_uniquedf[['Latitude', 'Longitude']]),
                                           how='left', on=['Latitude', 'Longitude'])

#Get netcdf value for every unique air monitoring station record (x, y, time) into new column
airdat_uniquedf_ij['date_index'] = \
    airdat_uniquedf_ij['date'].apply(lambda x:
                                     np.where(timevals == netCDF4.date2num(datetime.strptime(x, '%Y-%m-%d'), timev.units))[0][0])

airdat_uniquedf_ij['crain'] =\
    airdat_uniquedf_ij.apply(lambda x: crainv[x.date_index, int(x['iy_min']), int(x['ix_min'])], axis=1)