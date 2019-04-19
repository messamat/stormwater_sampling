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
from datetime import datetime, timedelta
import timeit
import glob
import dateparser
import xarray as xr
#Custom functions
from Download_gist import *
import cPickle as pickle

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
monitortab = os.path.join(AQIdir, 'aqs_monitors.csv')
sitetab = os.path.join(AQIdir, 'aqs_sites.csv')
if not os.path.exists(monitortab):
    print('Monitor tab does not exist...')
    dlfile('https://aqs.epa.gov/aqsweb/airdata/aqs_monitors.zip', AQIdir)
monitors = pd.read_csv(monitortab)
if not os.path.exists(sitetab):
    print('Site tab does not exist...')
    dlfile('https://aqs.epa.gov/aqsweb/airdata/aqs_sites.zip', AQIdir)
sites = pd.read_csv(sitetab)

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

airdatall = os.path.join(AQIdir, 'daily_SPEC_collate.csv')
#airdat_uniquetab = os.path.join(AQIdir, 'daily_SPEC_unique.csv')
airdat_uniquedfexp_pickle =  os.path.join(rootdir, 'results/airdat_uniquedfexp.p')

#Functions
def pdtogpd_datum(df, datum):
    "Subset df based on datum and create Geodataframe with the right assigned datum"
    crs = {'WGS84': {'init': 'epsg:4326'}, 'NAD27': {'init': 'epsg:4267'},'NAD83': {'init': 'epsg:6318'}} #'NAD83': {'init': 'epsg:1116'}}
    subdat = df[df['Datum']==datum]
    geom = [Point(xy) for xy in zip(subdat.Longitude, subdat.Latitude)]
    return(gpd.GeoDataFrame(subdat, crs=crs[datum], geometry=geom))

def getclosest_ij(X, Y, Xpt, Ypt):
    # find squared distance of every point on grid
    dist_sq = np.square(X - Xpt) + np.square(Y - Ypt)
    # 1D index of minimum dist_sq element
    minindex_flattened = dist_sq.argmin()
    # Get 2D index for latvals and lonvals arrays from 1D index
    return np.unravel_index(minindex_flattened, X.shape)

def getclosest_ij_df(X, Y, XYpt):
    """# A function to find the index of the point closest pt(in squared distance) to give lat/lon value.
    latlonpt may be a 2-column dataframe (lat-lon) or a list/tuple of lat=lon for a point.
    If data.frame, returns a data.frame of unique coordinates from the original data.frame together with index columns"""

    if isinstance(XYpt, pd.DataFrame):
        print('x and y to locate are panda series, proceeding in pd framework...')
        if XYpt.shape[1] == 2:
            #Get unique set of coordinates
            xypt_unique = XYpt.drop_duplicates()
            #Get closest grid for each unique coordinate and return data.frame to join to original one
            return xypt_unique.apply(
                lambda x: pd.Series([x[0], x[1]] +
                                    list(getclosest_ij(X, Y, float(x[0]), float(x[1]))),
                                    index = list(XYpt.columns.values) + ['ix_min', 'iy_min']), axis=1)
        else:
            raise ValueError('Wrong number of columns in XYpt')

def extractCDFtoDF(indf, pattern, indir, varname, keepcols=None, datecol = None, level=None, outfile=None, overwrite=False):
    if not (((outfile != None and os.path.exists(outfile)) or
             (outfile == None and os.path.exists(indf))) and
            overwrite == False):
        """Warning: will edit source df with extra columns"""
        pathpattern = os.path.join(indir, pattern)
        cdflist = glob.glob(pathpattern)
        print('Extracting {0}'.format(cdflist))

        #Get netCDF data
        try:
            sourcef = xr.open_mfdataset(cdflist)
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
        for dim in sourcef.dims:
            print('{0}, length: {1}'.format(dim, sourcef.dims[dim]))

        #If supplied path to csv, read csv as df
        if isinstance(indf, str) and os.path.exists(indf):
            if os.path.splitext(indf)[1] == '.csv':
                df = pd.read_csv(indf)
            elif os.path.splitext(indf)[1] == '.p':
                df = pickle.load(open(indf, "rb"))
            #Could use .feather too
            # elif os.path.splitext(indf)[1] == '.feather':
            #     df = pd.read_feather(indf)
            else:
                raise ValueError('Input df argument unknown format, only pd dfs, .csv and .p are supported')
        elif isinstance(indf, pd.DataFrame):
            df = indf
        else:
            raise ValueError('Input df argument neither csv nor dataframe ')

        #Make sure that right columns are in df
        if datecol is None:
            datecol = 'date'
        if any([c not in df.columns for c in ['x', 'y', datecol]]):
            raise ValueError('Column missing in df; either x, y or date')

        #Get index of the pixels closest to each point in df
        print('Getting index of pixels closest to points in df...')

        df_ij = df.merge(getclosest_ij_df(
            X=np.tile(sourcef['x'].values, [sourcef['y'].shape[0], 1]),
            Y=np.tile(np.transpose(np.asmatrix(sourcef['y'].values[::-1])), [1, sourcef['x'].shape[0]]),
            XYpt=df[['x', 'y']]),
            how='left', on=['x', 'y'])

        #Convert dates in df to index
        print('Converting dates to index...')
        if df_ij[datecol].dtype.name != 'datetime64[ns]':
            print('First converting {} to date datetime64 format'.format(datecol))
            df_ij[datecol] = dateparser.parse(datecol)
        df_ij = df_ij.merge(pd.DataFrame({'date_index': range(sourcef.dims['time']),
                                          datecol: sourcef['time'].values}),
                            on=datecol)

        #Extract variables
        print('Extracting variable...')
        if level is not None:
            if 'level' in sourcef.dimensions:
                level_i = np.where(sourcef.variables['level'][:] == level)[0][0]
                varname = '{0}{1}'.format(varname, level)
                df_ij[varname] = (sourcef[varname].isel(time=xr.DataArray(df_ij.date_index, dims='z'),
                                                        x=xr.DataArray(df_ij.ix_min.astype(int), dims='z'),
                                                        y=xr.DataArray(df_ij.iy_min.astype(int), dims='z'),
                                                        level=level_i)).values
            else:
                raise ValueError("A 'level' argument was provided but netCDF does not have a level dimension")
        else:
            df_ij[varname] = (sourcef[varname].isel(time=xr.DataArray(df_ij.date_index, dims='z'),
                                                    x=xr.DataArray(df_ij.ix_min.astype(int), dims='z'),
                                                    y=xr.DataArray(df_ij.iy_min.astype(int), dims='z'))).values

        #Delete intermediate columns
        print('Deleting indices columns')
        df_ij.drop(['ix_min','iy_min','date_index'], axis=1, inplace=True)
        if keepcols is not None:
            df_ij.drop([c for c in df_ij.columns if c not in keepcols + list(varname)], axis=1, inplace=True)

        #Write df out
        if (isinstance(indf, str) and os.path.exists(indf)) or (outfile is not None):
            print('Writing df out...')
            outf = indf if outcsv is None else outfile
            if os.path.splitext(outf)[1] == '.csv':
                df_ij.to_csv(outf, sep='\t', encoding = 'utf-8')
            elif os.path.splitext(outf)[1] == '.p':
                pickle.dump(df_ij, open(outf, "wb" ))
            # elif os.path.splitext(outfile)[1] == '.feather':
            #     df_ij.to_feather(outfile)
    else:
        raise ValueError('{} already exists and overwrite==False, '
                         'either set overwrite==True or change outfile'.format(outfile))

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
#paramsel = ['11','12','14','22','65','82','85','86','89']
#param_regex = re.compile('|'.join(['(^{}.*)'.format(p) for p in paramsel]))
#(monitors['Parameter Code'].astype(str).str.contains(param_regex)) &
monitors_chem = monitors.loc[(monitors['Parameter Name'].str.contains(elems_regex)) &
                             ~(monitors['Parameter Name'].str.contains(elemsout_regex)),:]


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
# DOWNLOAD DAILY METEOROLOGICAL DATA
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
yearlist = range(2014, 2020)
for var in varfolder_dic:
    if varfolder_dic[var] == 'pressure':
        varbase = var.split('.')[0]
        downloadNARR(folder=varfolder_dic[var], variable=varbase,
                     years=['{0}{1}'.format(year,month) for year in yearlist
                            for month in [str(i).zfill(2) for i in range(1,13)]],
                     outdir=NARRdir)
    else:
        downloadNARR(folder=varfolder_dic[var], variable=var, years=yearlist, outdir=NARRdir)

#Get proj4 from netcdf
templatef = xr.open_mfdataset(glob.glob(os.path.join(NARRdir, varfolder_dic.keys()[0]+'*')))
crsatt = templatef.variables['Lambert_Conformal'].attrs
lcc_proj4 = ("+proj=lcc +lat_1={0} +lat_2={1} +lat_0={2} +lon_0={3} +x_0={4} +y_0={5} +units=m +no_def".
             format(crsatt['standard_parallel'][0],
                    crsatt['standard_parallel'][1],
                    crsatt['latitude_of_projection_origin'],
                    crsatt['longitude_of_central_meridian'],
                    crsatt['false_easting'],
                    crsatt['false_northing']))

#-----------------------------------------------------------------------------------------------------------------------
# CREATE POINT SHAPEFILE AND BUFFER
#-----------------------------------------------------------------------------------------------------------------------
sites_chem.groupby('Datum')['UID'].nunique() #Check what datums data are in

#Create geodataframes
sites_wgs84 = pdtogpd_datum(sites_chem, 'WGS84')
sites_nad27 = pdtogpd_datum(sites_chem, 'NAD27')
sites_nad83 = pdtogpd_datum(sites_chem, 'NAD83')

#Project them all to same projection as meteorological data (Lambert Comformal Conic)
sites_gpd = pd.concat([sites_wgs84.to_crs(crs=lcc_proj4), sites_nad27.to_crs(crs=lcc_proj4), sites_nad83.to_crs(crs=lcc_proj4)])
sites_gpd['x'], sites_gpd['y'] = sites_gpd.centroid.x, sites_gpd.centroid.y
#Remove outliers with lon=0
sites_gpd['Longitude'].describe()
sites_gpd = sites_gpd[sites_gpd['Longitude'] != 0]

#Create 550 m radius buffers and dissolve them
sites_buf = sites_gpd['geometry'].buffer(distance=550)
sites_bufgpd = gpd.GeoDataFrame(sites_buf, crs=lcc_proj4).rename(columns={0:'geometry'}).set_geometry('geometry')

#Dissolve all
sites_bufgpd['diss'] = 1
sites_bufdis = sites_bufgpd.dissolve(by='diss')

#Dissolve only overlapping buffers
sites_bufunion = gpd.GeoDataFrame([polygon for polygon in sites_bufgpd.unary_union], crs=lcc_proj4). \
    rename(columns={0:'geometry'}).set_geometry('geometry')

#Output to shapefile
sites_gpd.to_file(sites_out, driver = 'ESRI Shapefile')
sites_bufgpd.to_file(sites_outbuf, driver = 'ESRI Shapefile')
sites_bufdis.to_file(sites_outbufdis, driver = 'ESRI Shapefile')
sites_bufunion.to_file(sites_outbufunion, driver = 'ESRI Shapefile')

#-----------------------------------------------------------------------------------------------------------------------
# DOWNLOAD SPECIATION DATA FOR ALL MONITORING SITES
#-----------------------------------------------------------------------------------------------------------------------
#Download data
yearlist = range(2014, 2020)
epadl_url = "https://aqs.epa.gov/aqsweb/airdata/"
spec25_list = [os.path.join(epadl_url, "daily_SPEC_{}.zip".format(year)) for year in yearlist]
spec10_list = [os.path.join(epadl_url, "daily_PM10SPEC_{}.zip".format(year)) for year in yearlist]
for file in spec25_list + spec10_list:
    dlfile(file, outpath=os.path.join(AQIdir))

#Collate all data
if not os.path.exists(airdatall):
    mergedel(AQIdir, 'daily_.*SPEC.*[0-9]{4,6}[.]csv$', airdatall, verbose=True)
airdat_df = pd.read_csv(airdatall)

#-----------------------------------------------------------------------------------------------------------------------
# FORMAT AQI DATA
#-----------------------------------------------------------------------------------------------------------------------
if not os.path.exists(airdat_uniquedfexp_pickle):
    print('Compute unique ID for each monitoring site...')
    airdat_df['UID'] = airdat_df['State Code'].astype(str).str.zfill(2) + \
                       airdat_df['County Code'].astype(str) + \
                       airdat_df['Site Num'].astype(str)

    airdat_df[airdat_df['UID']=='01101102']
    print('Compute df of unique UID-time records to extract meteorological variables with...')
    subcols = ['UID', 'Latitude', 'Longitude', 'Date Local']
    airdat_unique = pd.unique(airdat_df[subcols[0]].str.cat(airdat_df[subcols[1:]].astype(str), sep='._.')) #Concatenate columns

    airdat_uniquedf = pd.DataFrame(
        [i.split('._.') for i in list(airdat_unique)],
        columns = subcols) #Resplit and format to df
    #Faster than airdat_unique.str.split("._.", n = 4, expand = True) ?
    airdat_uniquedf['date'] = airdat_uniquedf['Date Local']

    # Get projected x and y coordinates for each site
    airdat_uniquedfproj = airdat_uniquedf.merge(pd.DataFrame(sites_gpd.drop(columns='geometry')), on='UID', how='inner')

    # Expand df to extract 3-hourly meteorological averages
    airdat_uniquedfexp = airdat_uniquedfproj.reindex(airdat_uniquedfproj.index.repeat(8))
    airdat_uniquedfexp['datetime_dupli'] = pd.to_datetime(airdat_uniquedfexp['date'], format='%Y-%m-%d') + \
                                           pd.to_timedelta(airdat_uniquedfexp.groupby(level=0).cumcount() * 3, unit='h')

    pickle.dump(airdat_uniquedfexp, open(airdat_uniquedfexp_pickle, "wb"))
else:
    airdat_uniquedfexp = pickle.load(open(airdat_uniquedfexp_pickle, "rb"))

#-----------------------------------------------------------------------------------------------------------------------
# EXTRACT ALL METEOROLOGICAL VARIABLES FOR EACH STATION-DATE COMBINATION
#-----------------------------------------------------------------------------------------------------------------------
#Dictionary of variables to extract
cdfvardict = {
    'air.2m': ['air.2m*.nc', None],
    'air.sfc': ['air.sfc*.nc', None],
    'crain': ['crain.*.nc', None],
    'dswrf': ['dswrf.*.nc', None],
    'hgt.850': ['hgt.*.nc', 850],
    'hpbl': ['hpbl.*.nc', None],
    'lftx4': ['lftx4.*.nc', None],
    'pres.sfc': ['pres.sfc.*.nc', None],
    'rhum.2m': ['rhum.2m.*.nc', None],
    'shum.2m': ['shum.2m.*.nc', None],
    'uwnd.10m': ['uwnd.10m.*.nc', None],
    'vwnd.10m': ['vwnd.10m.*.nc', None],
    'vwnd.500': ['vwnd.2*.nc', 500]
}

#Iterate through variables to extract
for var in cdfvardict.keys():
    try:
        print('Processing {}...'.format(var))
        extractCDFtoDF(airdat_uniquedfexp_pickle, pattern= cdfvardict[var][0], indir=NARRdir, varname= var,
                       datecol = 'datetime_dupli', level=cdfvardict[var][1],
                       keepcols = ['UID', 'datetime_dupli', 'x', 'y'], overwrite=False,
                       outfile= os.path.join(rootdir, 'results/daily_SPEC_{}.p'.format(var.replace('.', '_'))))
    except Exception:
        traceback.print_exc()
        print('Skipping...')
        pass

#Extract data



#------------------------------------------------------------------------------------------------------
# COMPUTE DERIVED VARIABLES
#-----------------------------------------------------------------------------------------------------------------------
cdfvardict.keys()
covar_pm25


##############################################################################################################
####################################################################################################################################
# EXTRA STUFF
air2m = netCDF4.Dataset(os.path.join(NARRdir, 'air.2m.2016.nc'))
# air2m = netCDF4.Dataset(os.path.join(NARRdir, 'air.2m.2016.nc'))
# air2m.close()
# #Check variables
for a in air2m.variables:
    print(a)
    varobj = air2m.variables[a]
    print(varobj)

air2m.variables['Lambert_Conformal']