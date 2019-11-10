import arcpy
import mendeleev
import pandas as pd
import itertools
import os
import re
import geopandas as gpd
from shapely.geometry import Point
import us
import numpy as np
import netCDF4
from datetime import datetime, timedelta
import time
import glob
import dateparser
import xarray as xr
import cPickle as pickle
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

#Custom functions
from xarrayDataAssembly import *
from Download_gist import *

#Set options
#arcpy.CheckOutExtension("Spatial")
#arcpy.env.overwriteOutput=True
#arcpy.env.qualifiedFieldNames = False

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'
AQIdir = os.path.join(rootdir, 'data/EPA_AirData_201902')
NARRdir = os.path.join(rootdir, 'data/NARR_201904')
if not os.path.isdir(NARRdir):
    os.mkdir(NARRdir)
NARRoutdir = os.path.join(rootdir, 'results/NARR')
if not os.path.isdir(NARRoutdir):
    os.mkdir(NARRoutdir)


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
NLCD_imp = os.path.join(rootdir, 'data/NLCD_2016_Impervious_L48_20190405.img') #Based on 2016 dara
ref_cs = arcpy.Describe(NLCD_imp).SpatialReference

#Output variables
outdir = os.path.join(rootdir, 'results/airdata')
if not os.path.isdir(outdir):
    os.mkdir(outdir)
sites_out = os.path.join(outdir, 'airsites.shp')
sites_out_lambers = os.path.join(outdir, 'airsites_lambers.shp')
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
def save_rdata_file(df, filename):
    #Save panda df to rdata file
    #https://stackoverflow.com/questions/15081858/can-i-use-rpy2-to-save-a-pandas-dataframe-to-an-rdata-file
    r_data = pandas2ri.py2ri(df)
    robjects.r.assign("my_df", r_data)
    robjects.r("save(my_df, file='{}')".format(filename))

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

def subsetNARRlevel(indir, pattern, sel_level, outnc):
    pathpattern = os.path.join(indir, pattern)
    cdflist = glob.glob(pathpattern)
    print('Extracting {0}'.format(cdflist))

    try:
        sourcef = xr.open_mfdataset(cdflist)
    except Exception as e:
        traceback.print_exc()
        if isinstance(e, RuntimeError):
            for fpath in cdflist:
                try:
                    netCDF4.Dataset(fpath)
                except:
                    print('Error stems from {}...'.format(fpath))
                    pass

    #level_i = np.where(sourcef['level'] == level)[0][0]
    sourcef.sel(level=sel_level).to_netcdf(outnc)

def extractCDFtoDF(indf, pattern, indir, varname, keepcols=None, datecol = None, level=None, outfile=None, overwrite=False):
    if not (((outfile != None and os.path.exists(outfile)) or
             (outfile == None and os.path.exists(indf))) and
            overwrite == False):
        pathpattern = os.path.join(indir, pattern)
        cdflist = glob.glob(pathpattern)
        print('Extracting {0}'.format(cdflist))

        #Get netCDF data
        try:
            sourcef = xr.open_mfdataset(cdflist)
        except Exception as e:
            traceback.print_exc()
            if isinstance(e, RuntimeError):
                for fpath in cdflist:
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
            if 'level' in sourcef.dims:
                level_i = np.where(sourcef['level'] == level)[0][0]
                df_ij[varname] = (sourcef[varname.split('_')[0]].isel(
                    time=xr.DataArray(df_ij.date_index, dims='z'),
                    x=xr.DataArray(df_ij.ix_min.astype(int), dims='z'),
                    y=xr.DataArray(df_ij.iy_min.astype(int), dims='z'),
                    level=level_i)).values
            else:
                raise ValueError("A 'level' argument was provided but netCDF does not have a level dimension")
        else:
            df_ij[varname] = (sourcef[varname.split('_')[0]].isel(
                time=xr.DataArray(df_ij.date_index, dims='z'),
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
            outf = indf if outfile is None else outfile
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
sites['UID'] = sites['State Code'].astype(str) + \
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

#Project them all to same projection as meteorological data for extraction (Lambert Comformal Conic)
sites_gpd_lambers = pd.concat([sites_wgs84.to_crs(crs=lcc_proj4), sites_nad27.to_crs(crs=lcc_proj4), sites_nad83.to_crs(crs=lcc_proj4)])
sites_gpd_lambers['x'], sites_gpd_lambers['y'] = sites_gpd_lambers.geometry.x, sites_gpd_lambers.geometry.y
#Remove outliers with lon=0
sites_gpd_lambers['Longitude'].describe()
sites_gpd_lambers = sites_gpd_lambers[sites_gpd_lambers['Longitude'] != 0]

#Project them all to same projection as NLCD data (and all other data)
aea_proj4 = ("+proj=aea +lat_1={0} +lat_2={1} +lat_0={2} +lon_0={3} +x_0={4} +y_0={5} "
             "+ellps=GRS80 +datum=NAD83 +units=m +no_def".
             format(ref_cs.standardParallel1,
                    ref_cs.standardParallel2,
                    ref_cs.latitudeOfOrigin,
                    ref_cs.centralMeridian,
                    ref_cs.falseEasting,
                    ref_cs.falseNorthing))
sites_gpd = pd.concat([sites_wgs84.to_crs(crs=lcc_proj4), sites_nad27.to_crs(crs=lcc_proj4), sites_nad83.to_crs(crs=lcc_proj4)])
sites_gpd['x'], sites_gpd['y'] = sites_gpd.geometry.x, sites_gpd.geometry.y
#Remove outliers with lon=0
sites_gpd['Longitude'].describe()
sites_gpd = sites_gpd_lambers[sites_gpd_lambers['Longitude'] != 0]

#Create 550 m radius buffers and dissolve them
sites_buf = sites_gpd['geometry'].buffer(distance=550)
sites_bufgpd = gpd.GeoDataFrame(sites_buf, crs=aea_proj4).rename(columns={0:'geometry'}).set_geometry('geometry')

#Dissolve all
sites_bufgpd['diss'] = 1
sites_bufdis = sites_bufgpd.dissolve(by='diss')

#Dissolve only overlapping buffers
sites_bufunion = gpd.GeoDataFrame([polygon for polygon in sites_bufgpd.unary_union], crs=lcc_proj4). \
    rename(columns={0:'geometry'}).set_geometry('geometry')

#Output to shapefile
sites_gpd.to_file(sites_out, driver = 'ESRI Shapefile')
sites_gpd_lambers.to_file(sites_out_lambers, driver = 'ESRI Shapefile')
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
    airdat_uniquedfproj = airdat_uniquedf.merge(pd.DataFrame(sites_gpd_lambers.drop(columns='geometry')), on='UID', how='inner')

    # Expand df to extract 3-hourly meteorological averages
    airdat_uniquedfexp = airdat_uniquedfproj.reindex(airdat_uniquedfproj.index.repeat(8))
    airdat_uniquedfexp['datetime_dupli'] = pd.to_datetime(airdat_uniquedfexp['date'], format='%Y-%m-%d') + \
                                           pd.to_timedelta(airdat_uniquedfexp.groupby(level=0).cumcount() * 3, unit='h')

    pickle.dump(airdat_uniquedfexp, open(airdat_uniquedfexp_pickle, "wb"))
else:
    airdat_uniquedfexp = pickle.load(open(airdat_uniquedfexp_pickle, "rb"))


#------------------------------------------------------------------------------------------------------
# COMPUTE DERIVED VARIABLES
#-----------------------------------------------------------------------------------------------------------------------
covar_pm25

#Get fired data
#Porter et al. 2015 used: MODIS Global Monthly Fire Location Product
#https://www.ospo.noaa.gov/Products/land/hms.html
#https://satepsanone.nesdis.noaa.gov/pub/FIRE/HMS/GIS/ARCHIVE/
#https://www.ospo.noaa.gov/Products/land/fire.html: try ASDTA Smoke-East AOD and ASDTA Smoke-West AOD

#Data request:
#https://www.ncdc.noaa.gov/has/has.dsselect
#GOES-13, SMOKEE_GRD 2014/01/01 - 2018/08/01;
#       Check status at: https://www.ncdc.noaa.gov/has/HAS.CheckOrderStatus?hasreqid=HAS011301876&emailadd=messamat@uw.edu
#       Delivery Location: http://www1.ncdc.noaa.gov/pub/has/HAS011301876/
#GOES-15, SMOKEW_GRD 2014/01/01 - 2019/04/17
#       Check status at: https://www.ncdc.noaa.gov/has/HAS.CheckOrderStatus?hasreqid=HAS011301886&emailadd=messamat@uw.edu
#       Delivery location: http://www1.ncdc.noaa.gov/pub/has/HAS011301886/

#-------------- subset datasets by pressure level
subsetNARRlevel(indir=NARRdir, pattern = 'hgt.*.nc', sel_level = 850, outnc = os.path.join(NARRoutdir, 'hgt.850.nc'))
subsetNARRlevel(indir=NARRdir, pattern = 'vwnd.2*.nc', sel_level = 500, outnc = os.path.join(NARRoutdir, 'vwnd.500.nc'))

#-------------- Generate 9x9 mean -----------------#
#crain_9x9
for yeardat in glob.glob(os.path.join(NARRdir, 'crain.*.nc')):
    print('Processing {}...'.format(yeardat))
    outdat = os.path.join(NARRoutdir, os.path.splitext(os.path.split(yeardat)[1])[0]+'_9x9.nc')
    if not os.path.exists(outdat):
        crainf = xr.open_dataset(yeardat, chunks={'time': 10})
        crainroll = crainf['crain'].rolling(y=9, center=True).mean().rolling(x=9, center=True).mean()
        crainroll.to_netcdf(outdat)
    else:
        print('{} already exists, skipping...'.format(outdat))

#air.sfc 9x9
for yeardat in glob.glob(os.path.join(NARRdir, 'air.sfc*.nc')):
    print('Processing {}...'.format(yeardat))
    outdat = os.path.join(NARRoutdir, os.path.splitext(os.path.split(yeardat)[1])[0]+'_9x9.nc')
    if not os.path.exists(outdat):
        crainf = xr.open_dataset(yeardat, chunks={'time': 10})
        crainroll = crainf['air'].rolling(y=9, center=True).mean().rolling(x=9, center=True).mean()
        crainroll.to_netcdf(outdat)
    else:
        print('{} already exists, skipping...'.format(outdat))

#-------------- Compute day and night extrema -----------------#
arr = xr.open_dataset(os.path.join(NARRdir, 'air.sfc.2015.nc'))
#Add 8h to datetime so that e.g. 02/01 16h-02/02 8h is 02/02 00*02/02 16h (preceding night is now consider part of that date
arr['shifted_datetime'] = arr['time'].values + pd.to_timedelta(timedelta(hours=8))
arr['shifted_date'] = arr.shifted_datetime.dt.floor('d')
#Assign night vs day
arr.coords['daynight'] = xr.where(arr.time.dt.hour.isin(range(0,17)), 'night', 'day')
arr.coords['dateperiod'] = pd.Series(np.datetime_as_string(arr.shifted_date.values)).str[0:10] + arr.daynight.values

#Add dateperiod dimension (https://stackoverflow.com/questions/39626402/add-dimension-to-an-xarray-dataarray)
#Stack over x, y, dateperiod
#Get groupby.mean
#unstack


stackedarr = arr.stack(allpoints=['x', 'y'])
# define a function to compute a linear trend of a timeseries
def meandateperiod(x):
    pf = np.polyfit(x.time, x, 1)
    # we need to return a dataarray or else xarray's groupby won't be happy
    return xr.DataArray(pf[0])
tic=time.time()
trend = stackedarr.groupby('allpoints').apply(linear_trend)
print(time.time() - tic)

#arr.expand_dims(dim='dateperiod')





tic=time.time()
arr_periodmean = arr.groupby('dateperiod').mean()
print(time.time() - tic)




tic=time.time()
results = []
for label1, group1 in arr.groupby('shifted_date'):
    for label2, group2, in group1.groupby('daynight'):
        results.append(group2.mean())
xr.concat(results, dim=['shifted_date','daynight'])
print(time.time() - tic)





arr = xr.DataArray(np.arange(0, 625, 1).reshape(25, 25),dims=('x', 'y'))
r= arr.rolling(y=3)
r.mean(skipna=True)
for label, arr_window in r:
    print(arr_window+1)



#get mean, min, max groupby('date', 'dperiod')




#air.sfc_9x9 - nightmin
#crain_9x9 - nightmax
#lts_daymin
#uwnddir.10m - daymean
#uwnddir.19m_nightmean
#vwnd.500 - daymax
#wspd.10m-daymax
#wspd.10m-nightmax
#dswrf-daymin
#lftx4-nightmin
#shum.2m-nightmin


#Compute 3- and 6-day maxima, minima, and means

#Compute 1-day delta variable



#Compute daily recirculation potential index (RPI)
"surface wind speeds based on the ratio between the vector sum magnitude (L) and scalar sum (S) of wind speeds over the" \
"previous 24 h (Levy et al., 2009)"

#Compute lower-tropospheric stability (LTS)
# difference between surface and 700 hPa potential temperatures


#-----------------------------------------------------------------------------------------------------------------------
# EXTRACT ALL METEOROLOGICAL VARIABLES FOR EACH STATION-DATE COMBINATION
#-----------------------------------------------------------------------------------------------------------------------
#Dictionary of variables to extract (do not rely on level argument for hgt and vwnd as netcdfs are too big and saturate memory_
cdfvardict = {
    'air_2m': ['air.2m*.nc', None],
    'air_sfc': ['air.sfc*.nc', None],
    'crain': ['crain.*.nc', None],
    'dswrf': ['dswrf.*.nc', None],
    'hgt_850': ['hgt.850*.nc', None],
    'hpbl': ['hpbl.*.nc', None],
    'lftx4': ['lftx4.*.nc', None],
    'pres_sfc': ['pres.sfc.*.nc', None],
    'rhum_2m': ['rhum.2m.*.nc', None],
    'shum_2m': ['shum.2m.*.nc', None],
    'uwnd_10m': ['uwnd.10m.*.nc', None],
    'vwnd_10m': ['vwnd.10m.*.nc', None],
    'vwnd_500': ['vwnd.500*.nc', None]
}

#Iterate through variables to extract
for var in ['pres_sfc', 'shum_2m']: #cdfvardict.keys():
    try:
        print('Processing {}...'.format(var))
        extractCDFtoDF(indf = airdat_uniquedfexp_pickle, pattern= cdfvardict[var][0], indir=NARRdir, varname=var,
                       datecol = 'datetime_dupli', level=cdfvardict[var][1],
                       keepcols = ['UID', 'datetime_dupli', 'x', 'y', var], overwrite=True,
                       outfile= os.path.join(rootdir, 'results/daily_SPEC_{}.p'.format(var.replace('.', '_'))))
    except Exception:
        traceback.print_exc()
        print('Skipping...')
        pass

for var in ['vwnd_500', 'hgt_850']:
    print('Processing {}...'.format(var))
    extractCDFtoDF(indf=airdat_uniquedfexp_pickle, pattern=cdfvardict[var][0], indir=NARRoutdir, varname=var,
                   datecol='datetime_dupli', level=cdfvardict[var][1],
                   keepcols=['UID', 'datetime_dupli', 'x', 'y', var], overwrite=True,
                   outfile=os.path.join(rootdir, 'results/daily_SPEC_{}.p'.format(var.replace('.', '_'))))

#Collate all air data
airdat_climmerge = pd.read_csv(airdatall)
airdat_climmerge['UID'] = airdat_climmerge ['State Code'].astype(str).str.zfill(2) + \
                   airdat_climmerge ['County Code'].astype(str) + \
                   airdat_climmerge ['Site Num'].astype(str)

for var in cdfvardict.keys():
    print(var)
    with open(os.path.join(rootdir, 'results/daily_SPEC_{}.p'.format(var.replace('.', '_'))), "rb") as input_file:
         indf = pickle.load(input_file)
    indf['Date Local'] = indf['datetime_dupli'].dt.date.astype(str)
    indf_stat = indf.groupby(['UID', 'Date Local'], as_index=False).agg({var: ["max", "mean"]})
    indf_stat.columns = ["".join(x) for x in indf_stat.columns.ravel()]
    airdat_climmerge = airdat_climmerge.merge(indf_stat, on=['UID', 'Date Local'])

    pickle.dump(airdat_climmerge, open(os.path.join(rootdir, "airdat_scratch"), "wb"))

#Write out data to table (should use feather but some module conflicts and don't want to deal with it)
airdat_climmerge.to_csv(os.path.join(rootdir, 'results/airdat_NARRjoin.csv'))



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