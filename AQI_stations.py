#Make sure to run this with 64bit Python. Even chunking won't help.

import mendeleev
import arcpy
import pandas as pd
import itertools
import os
import re
import geopandas as gpd
from shapely.geometry import Point
from collections import defaultdict
import us
import numpy as np
import netCDF4
from datetime import date, timedelta
import time
import glob
import traceback
import dask
import dateparser
import xarray as xr
import bottleneck # rolling window aggregations are faster and use less memory when bottleneck is installed http://xarray.pydata.org/en/stable/computation.html#rolling-window-operations
import cPickle as pickle
import rpy2
import feather
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

#Custom functions
from Download_gist import *

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'
AQIdir = os.path.join(rootdir, 'data/EPA_AirData_201902')
AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')
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
NLCD_imp = os.path.join(rootdir, 'data/NLCD_2016_Impervious_L48_20190405.img') #Based on 2016 data
cs_ref = arcpy.Describe(NLCD_imp).SpatialReference
smoke1419 = os.path.join(AQIgdb, 'smoke1419_aea')

#Output variables
outdir = os.path.join(rootdir, 'results/airdata')
if not os.path.isdir(outdir):
    os.mkdir(outdir)
sites_out = os.path.join(outdir, 'airsites.shp')
sites_out_lambers = os.path.join(outdir, 'airsites_lambers.shp')
sites_outbuf = os.path.join(outdir, 'airsites_550buf.shp')
sites_outbufdis = os.path.join(outdir, 'airsites_550bufdis.shp')
sites_outbufunion = os.path.join(outdir, 'airsites_550bufunion.shp')
sites_smokejoin = os.path.join(AQIgdb, 'sitessmokejoin')
sites_smokejointab = os.path.join(AQIgdb, 'sitessmokejoin_tab')
sites_smokejoin_subpd = os.path.join(AQIdir, 'sitessmokejoin_tabsub.p')

airdatall = os.path.join(AQIdir, 'daily_SPEC_collate.csv')
#airdat_uniquetab = os.path.join(AQIdir, 'daily_SPEC_unique.csv')
airdat_uniquedf_pickle =  os.path.join(rootdir, 'results/airdat_uniquedfproj.p')

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

def subsetNARRlevel(indir, pattern, sel_level, outnc, overwrite=False):
    if os.path.exists(outnc) and overwrite == False:
        print('{} already exists and overwrite=False, skipping...'.format(outnc))
    else:
        pathpattern = os.path.join(indir, pattern)
        cdflist = glob.glob(pathpattern)
        print('Extracting {0}'.format(cdflist))

        try:
            sourcef = xr.open_mfdataset(cdflist, chunks={'time': 10})
            sourcef.sel(level=sel_level).to_netcdf(outnc)
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

def narr_daynightstat(indir, regexpattern, outdir, dnlist = ['day', 'night'], statlist = ['mean', 'min', 'max']):
    for yeardat in getfilelist(indir, os.path.split(regexpattern)[1]):
        print('Processing {}...'.format(yeardat))
        outdat = os.path.join(outdir, '{}_'.format(os.path.splitext(os.path.split(yeardat)[1])[0]))

        if not glob.glob('{0}{1}*'.format(outdat, dnlist[0])):
            datname = os.path.splitext(os.path.split(yeardat)[1])[0]
            varname = datname.split('.')[0]

            xrd = xr.open_dataset(yeardat, chunks={'time': 10})
            # Add 8h to datetime so that e.g. 02/01 16h-02/02 8h is 02/02 00*02/02 16h (preceding night is now consider part of that date
            xrd.coords['shifted_datetime'] = ('time', xrd['time'].values + pd.to_timedelta(timedelta(hours=8)))
            xrd.coords['shifted_date'] = ('time', xrd.shifted_datetime.dt.floor('d'))
            # Assign night vs day
            xrd.coords['daynight'] = ('time', xr.where(xrd.time.dt.hour.isin(range(0, 17)), 'night', 'day'))

            #Drop pressure level coordinates
            if 'level' in xrd.coords._names:
                xrd = xrd.drop('level')

            #For each of 'day' and 'night' in dnlist
            for per in dnlist:
                #Output dataset name
                outnc = '{0}{1}'.format(outdat, per)

                #If xarray.DataArray was directly saved to netcdf before being converted manually to xarray.Dataset
                datavars = xrd.data_vars.keys()
                if len(datavars) == 1 and datavars[0] == '__xarray_dataarray_variable__':
                    varname = '__xarray_dataarray_variable__'

                if 'mean' in statlist:
                    xrd[varname].where(xrd.daynight == per, drop=True). \
                        groupby('shifted_date').mean(dim='time').to_netcdf('{}mean.nc'.format(outnc))
                if 'min' in statlist:
                    xrd[varname].where(xrd.daynight == per, drop=True). \
                        groupby('shifted_date').min(dim='time').to_netcdf('{}min.nc'.format(outnc))
                if 'max' in statlist:
                    xrd[varname].where(xrd.daynight == per, drop=True). \
                        groupby('shifted_date').max(dim='time').to_netcdf('{}max.nc'.format(outnc))
        else:
            print('{} already exists, skipping...'.format(outdat))

def narr_d36stat(indir, regexpattern, outdir, dstat = None, multidstat = None):
    xrlist = getfilelist(indir, os.path.split(regexpattern)[1])
    if len(xrlist) > 0:
        for yeardat in xrlist:
            print('Processing {}...'.format(yeardat))
            outdat = os.path.join(outdir, '{0}_{1}.nc'.format(os.path.splitext(os.path.split(yeardat)[1])[0], dstat))

            if multidstat is not None:
                outdat = '{0}_{1}day{2}.nc'.format(os.path.splitext(outdat)[0], multidstat.keys()[0], multidstat.values()[0])

            if not os.path.exists(outdat):
                datname = os.path.splitext(os.path.split(yeardat)[1])[0]
                varname = datname.split('.')[0]

                xrd = xr.open_dataset(yeardat)
                if 'shifted_date' in xrd.coords.keys():
                    xrd = xrd.rename({'shifted_date':'date'})
                    xrd = xrd.chunk({'date': 100, 'x': 40, 'y': 40})
                else:
                    xrd.coords['date'] = ('time', xrd.time.dt.floor('d'))
                    xrd = xrd.chunk({'time': 100, 'x': 40, 'y': 40})

                #Drop pressure level coordinates
                if 'level' in xrd.coords._names:
                    xrd = xrd.drop('level')

                #If xarray.DataArray was directly saved to netcdf before being converted manually to xarray.Dataset
                datavars = xrd.data_vars.keys()
                if len(datavars) == 1 and datavars[0] == '__xarray_dataarray_variable__':
                    varname = '__xarray_dataarray_variable__'

                #Compute daily stat
                if dstat is not None:
                    if dstat == 'mean':
                        xrd = xrd[varname].groupby('date').mean(dim='time')
                    elif dstat == 'min':
                        xrd = xrd[varname].groupby('date').min(dim='time')
                    elif dstat == 'max':
                        xrd = xrd[varname].groupby('date').max(dim='time')
                    elif dstat == 'diff':
                        xrd = xrd[varname].diff(dim='date', n=1)
                    else:
                        raise ValueError('Daily statistics is not mean, min, or max')
                    xrd = xrd.chunk({'date':100, 'x':100, 'y':100})

                if multidstat is not None:
                    if isinstance(xrd, xr.Dataset):
                        xrd = xrd[varname]

                    if multidstat.values()[0] == 'mean':
                        xrd = xrd.rolling(date=multidstat.keys()[0], center=False,
                                          min_periods=multidstat.keys()[0]).construct('window').mean('window')
                    elif multidstat.values()[0] == 'min':
                        xrd = xrd.rolling(date=multidstat.keys()[0], center=False,
                                          min_periods=multidstat.keys()[0]).construct('window').min('window')
                    elif multidstat.values()[0] == 'max':
                        xrd = xrd.rolling(date=multidstat.keys()[0], center=False,
                                          min_periods=multidstat.keys()[0]).construct('window').max('window')
                    else:
                        raise ValueError('Multiday statistics is not mean, min, or max over 3 or 6 days')

                print("Saving {} to netcdf. This might take a while...".format(outdat))
                xrd.to_dataset(name=datname.split('.')[0]).chunk({'date':100, 'x':100, 'y':100}).to_netcdf(outdat)
            else:
                print('{} already exists, skipping...'.format(outdat))
    else:
        raise ValueError('regexpattern does not correspond to any existing dataset')

def extractCDFtoDF(indf, incdf, indir, varname, keepcols=None, datecol = None, level=None, outfile=None, overwrite=False):
    """
    Function to extract netcdf values at a set of points and dates for a list of netcdf files.
    Two main inputs: a table of points with x, y, and a date column; a wildcard pattern to fetch netcdf files

    :param indf: input table containing dataframe with x and y coordinates of points to which netcdf values should be extracted.
                    can be .csv or .p file.
    :param incdf: either xarray loaded in memory or regex pattern to use to fetch netcdf (does not include full path)
    :param indir: full path of directory from which to fetch netcdf
    :param varname: name of data variable to extract to points
    :param keepcols: whether to drop (None) or keep other columns than the ones that were extracted from the netcdf
    :param datecol: date column in the dataframe
    :param level: if netcdf has multiple pressure levels, which level to extract
    :param outfile: if outfile is None, overwrites indf, otherwise writes to .csv or .p
    :param overwrite: False or True, whether to overwrite outfile
    """

    if not (((outfile != None and os.path.exists(outfile)) or
             (outfile == None and os.path.exists(indf))) and
            overwrite == False):


        if isinstance(incdf, xr.DataArray) or isinstance(incdf, xr.Dataset):
            sourcef = incdf
        else:
            if isinstance(incdf, str):
                cdflist = getfilelist(indir, incdf)
                print('Extracting {0}'.format(cdflist))
            elif isinstance(incdf, list):
                cdflist = incdf
            else:
                raise ValueError('incdf parameter is neither a DataArray, Dataset, regex pattern, or file list')

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

        if 'shifted_date' in sourcef.coords.keys():
            sourcef = sourcef.rename({'shifted_date': 'date'})

        # If xarray.DataArray was directly saved to netcdf before being converted manually to xarray.Dataset
        datavars = sourcef.data_vars.keys()
        if len(datavars) == 1 and datavars[0] == '__xarray_dataarray_variable__':
            sourcef = sourcef.rename({'__xarray_dataarray_variable__': varname})
        else:
            sourcef = sourcef.rename({varname.split('.')[0]: varname})

        #Print dimensions and their length
        # print('Dimensions:')
        # for dim in sourcef.dims:
        #     print('{0}, length: {1}'.format(dim, sourcef.dims[dim]))

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
            raise ValueError('Input df argument neither csv nor pickle or dataframe ')

        #Make sure that right columns are in df
        if datecol is None:
            datecol = 'date'
        if any([c not in df.columns for c in ['x', 'y', datecol]]):
            raise ValueError('Column missing in df; either x, y or date')

        #Get index of the pixels closest to each point in df
        #print('Getting index of pixels closest to points in df...')

        df_ij = df.merge(getclosest_ij_df(
            X=np.tile(sourcef['x'].values, [sourcef['y'].shape[0], 1]),
            Y=np.tile(np.transpose(np.asmatrix(sourcef['y'].values[::-1])), [1, sourcef['x'].shape[0]]),
            XYpt=df[['x', 'y']]),
            how='left', on=['x', 'y'])

        #Convert dates in df to index
        #print('Converting dates to index...')
        if df_ij[datecol].dtype.name != 'datetime64[ns]':
            print('First converting {} to date datetime64 format'.format(datecol))
            df_ij[datecol] = pd.to_datetime(df_ij[datecol], format='%Y-%m-%d')
        df_ij = df_ij.merge(pd.DataFrame({'date_index': range(sourcef.dims['date']),
                                          datecol: sourcef['date'].values}),
                            on=datecol)

        #Extract variables
        print('Extracting variable...')
        if level is not None:
            if 'level' in sourcef.dims:
                level_i = np.where(sourcef['level'] == level)[0][0]
                df_ij[varname] = (sourcef[varname].isel(
                    time=xr.DataArray(df_ij.date_index, dims='z'),
                    x=xr.DataArray(df_ij.ix_min.astype(int), dims='z'),
                    y=xr.DataArray(df_ij.iy_min.astype(int), dims='z'),
                    level=level_i)).values
            else:
                raise ValueError("A 'level' argument was provided but netCDF does not have a level dimension")
        else:
            df_ij[varname] = (sourcef[varname].isel(
                date=xr.DataArray(df_ij.date_index, dims='z'),
                x=xr.DataArray(df_ij.ix_min.astype(int), dims='z'),
                y=xr.DataArray(df_ij.iy_min.astype(int), dims='z'))).values

        #Delete intermediate columns
        #print('Deleting indices columns')
        df_ij.drop(['ix_min','iy_min','date_index'], axis=1, inplace=True)
        if keepcols is not None:
            df_ij.drop([c for c in df_ij.columns if c not in keepcols + list(varname)], axis=1, inplace=True)

        #Write df out
        print('Writing df out...')
        if os.path.splitext(outfile)[1] == '.csv':
            df_ij.to_csv(outfile, sep='\t', encoding = 'utf-8')
        elif os.path.splitext(outfile)[1] == '.p':
            pickle.dump(df_ij, open(outfile, "wb" ))
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
monitors_chem.loc[:, 'State Code'] = monitors_chem.loc[:, 'State Code'].astype(str).str.pad(2, 'left', '0')

#Remove sites outside of the conterminous US (to not exceed 50,000 tile limit)
fips_excluded = [us.states.lookup(abr).fips for abr in ['PR','AK','HI','VI','UM', 'MP', 'GU'] if us.states.lookup(abr)]

#Compute Unique Site Identifier and select sites that record chemical concentrations
monitors_chem.loc[:, 'UID'] = monitors_chem.loc[:, 'State Code'] + \
                              monitors_chem.loc[:, 'County Code'].astype(str) + \
                              monitors_chem.loc[:, 'Site Number'].astype(str) + \
                              monitors_chem.loc[:, "Latitude"].astype(str).str[-3:] + \
                              monitors_chem.loc[:, "Longitude"].astype(str).str[-3:]


sites.loc[:, 'UID'] = sites.loc[:, 'State Code'].astype(str) + \
                      sites.loc[:, 'County Code'].astype(str) + \
                      sites.loc[:, 'Site Number'].astype(str) + \
                      sites.loc[:, "Latitude"].astype(str).str[-3:] + \
                      sites.loc[:, "Longitude"].astype(str).str[-3:]

sites_chem = sites[(sites.loc[:, 'UID'].isin(pd.unique(monitors_chem.loc[:, 'UID']))) &
                   ~(sites.loc[:, 'State Code'].isin(fips_excluded))]
len(sites_chem)
len(pd.unique(sites_chem['UID'])) #Most sites are unique

#-----------------------------------------------------------------------------------------------------------------------
# DOWNLOAD DAILY NARR METEOROLOGICAL DATA
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
covar_sub = [var for var in covar_pm25_base if var not in ['lts', 'fire', 'uwnddir.10m', 'vwnddir.10m', 'wspd.10m', 'rpi']] + \
            ['uwnd.10m', 'vwnd.10m'] + ['pottmp.sfc', 'air.700'] #to compute lts

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

#-----------------------------------------------------------------------------------------------------------------------
# COMPUTE DERIVED VARIABLES
# Does not include wind direction
#-----------------------------------------------------------------------------------------------------------------------
#-------------- subset datasets by pressure level
subsetNARRlevel(indir=NARRdir, pattern = 'hgt.*.nc', sel_level = 850, outnc = os.path.join(NARRoutdir, 'hgt.850.nc'))
subsetNARRlevel(indir=NARRdir, pattern = 'vwnd.2*.nc', sel_level = 500, outnc = os.path.join(NARRoutdir, 'vwnd.500.nc'))
subsetNARRlevel(indir=NARRdir, pattern = 'air.20*.nc', sel_level = 700, outnc = os.path.join(NARRoutdir, 'air.700.nc'))

#-------------- compute derived variables -------------------
#lts
# (lower-tropospheric stability) = potential temperatures (700 hPa) - potential temperatures (surface)
#To compute potential temperature at 700 hPa: PT700 = Temperature @ 700 hPa * (standard pressure/700)^(gas constant/specific heat)
#with standard pressure = 1000 hPa and gas constant/specific heat =  0.286
#So lts = (temperature_700*(1000/700)^0.286)- potential temp_surface

air700 = xr.open_dataset(os.path.join(NARRoutdir, 'air.700.nc'), chunks={'time':20}).drop('level')
for yeardat in glob.glob(os.path.join(NARRdir, 'pottmp.sfc.*.nc')):
    print('Processing {}...'.format(yeardat))
    outdat = os.path.join(NARRoutdir,
                          'lts.{}.nc'.format(re.compile('[0-9]{4}').search(os.path.split(yeardat)[1]).group()))
    if not os.path.exists(outdat):
        with xr.open_dataset(yeardat, chunks={'time': 20}) as pottmp_sfc:
            airpottmp_join = xr.merge([air700, pottmp_sfc], join='inner')
            airpottmp_join['pottmp700'] = airpottmp_join.air*pow(1000.0/700.0, 0.286)
            lts = airpottmp_join.pottmp700 - airpottmp_join.pottmp
            ltsdt = lts.to_dataset(name  ='lts')
            print('Saving to netdcdf. This might take a while... {}'.format(outdat))
            ltsdt.to_netcdf(outdat)
            del lts
            del ltsdt
    else:
        print('{} already exists, skipping...'.format(outdat))

#wspd (windspeed at 10 m) and rpi (recirculation potential index) = 1 - L/S
#the vector sum magnitude (L) and scalar sum (S) of surface wind speeds over the previous 24 h (Allwine and Whitemane 1994)
#https://www.sciencedirect.com/science/article/abs/pii/1352231094900485
#L = sqrt((3*sumforpast24h(uwind))^2+(3*sumforpast24h(vwind))^2)
#S = 3*sumoverpast24h(uwind^2+vwind^2)^1/2
for yr in range(2014, 2020):
    print('Processing {}...'.format(yr))
    outwspd = os.path.join(NARRoutdir, 'wspd.10m.{}.nc'.format(yr))
    outrpi = os.path.join(NARRoutdir, 'rpi.{}.nc'.format(yr))

    if not (os.path.exists(outrpi) and os.path.exists(outwspd)):
        uvw = xr.merge([xr.open_dataset(os.path.join(NARRdir, 'uwnd.10m.{}.nc'.format(yr)),
                                        chunks={'x': 10, 'y': 10}),
                        xr.open_dataset(os.path.join(NARRdir, 'vwnd.10m.{}.nc'.format(yr)),
                                        chunks={'x': 10, 'y': 10}).vwnd])
        # Compute wspd
        if not os.path.exists(outwspd):
            wspd = pow(uvw.uwnd**2 + uvw.vwnd**2, 0.5)
            wspddt = wspd.to_dataset(name='wspd')
            print('Saving {} to netcdf. This might take hours...'.format(outwspd))
            wspddt.to_netcdf(outwspd)
            del wspd
            del wspddt

        # Compute rpi
        if not os.path.exists(outrpi):
            #Used .construct('window').sum('window') rather than simply .sum() due to MemoryError issues (see https://github.com/pydata/xarray/issues/3165)
            rpi = (pow((3*uvw.uwnd.rolling(time=8, center=False, min_periods=8).construct('window').sum('window')**2) +
                       (3*uvw.vwnd.rolling(time=8, center=False, min_periods=8).construct('window').sum('window')**2), 0.5)/
                   3*pow(uvw.uwnd**2 + uvw.vwnd**2, 0.5).rolling(time=8, center=False, min_periods=8).construct('window').sum('window'))

            # Can easily lead to memory error. Switching to single-threaded scheduler did not help
            # Need to chunk but will likely also lead to error if chunking across time — must mess with the rolling
            # Lazy computing pushes all calculations to .to_netcdf
            # If does not work again. Try doing a first rolling, saving. Then another rolling, saving. Then combining, etc.
            rpidt = rpi.to_dataset(name  ='rpi')
            print('Saving {} to netcdf. This might take hours...'.format(outrpi))
            rpidt.to_netcdf(outrpi)

            del rpi
            del rpidt
            del uvw

    else:
        print('{} already exists, skipping...'.format(outdat))

#crain_9x9
for yeardat in glob.glob(os.path.join(NARRdir, 'crain.*.nc')):
    print('Processing {}...'.format(yeardat))
    outdat = os.path.join(NARRoutdir, '{}_9x9.nc'.format(os.path.splitext(os.path.split(yeardat)[1])[0]))
    if not os.path.exists(outdat):
        with xr.open_dataset(yeardat, chunks={'time':292}) as crainf:
            crainroll = crainf['crain'].rolling(y=9, center=True).mean().rolling(x=9, center=True).mean()
            crainrolldt = crainroll.to_dataset(name='crain')
            crainrolldt.to_netcdf(outdat)
        del crainroll
        del crainrolldt
    else:
        print('{} already exists, skipping...'.format(outdat))

#air.sfc 9x9
for yeardat in glob.glob(os.path.join(NARRdir, 'air.sfc*.nc')):
    print('Processing {}...'.format(yeardat))
    outdat = os.path.join(NARRoutdir, '{}_9x9.nc'.format(os.path.splitext(os.path.split(yeardat)[1])[0]))
    if not os.path.exists(outdat):
        with xr.open_dataset(yeardat, chunks={'time':5}) as airsfcf:
            airsfcroll = airsfcf['air'].rolling(y=9, center=True).mean().rolling(x=9, center=True).mean()
            airsfcrolldt = airsfcroll.to_dataset(name='air')
            airsfcrolldt.to_netcdf(outdat)
        del airsfcroll
    else:
        print('{} already exists, skipping...'.format(outdat))

#Night/Day stats
narr_daynightstat(indir=NARRdir, regexpattern='shum[.]2m.*[.]nc',
                  outdir=NARRoutdir, statlist=['mean', 'min'])
narr_daynightstat(indir=NARRdir, regexpattern='dswrf.*[.]nc',
                  outdir=NARRoutdir, dnlist=['day'], statlist=['min'])
narr_daynightstat(indir=NARRdir, regexpattern='lftx4.*[.]nc',
                  outdir=NARRoutdir, dnlist=['night'], statlist=['min'])
narr_daynightstat(indir=NARRdir, regexpattern='lts.*[.]nc',
                  outdir=NARRoutdir, dnlist=['day'], statlist=['min'])
narr_daynightstat(indir=NARRoutdir, regexpattern='wspd[.]10m[.][0-9]{4}[.]nc',
                  outdir=NARRoutdir, statlist=['mean', 'min', 'max'])
narr_daynightstat(indir=NARRoutdir, regexpattern='crain[.][0-9]{4}[_]9x9[.]nc',
                  outdir=NARRoutdir, dnlist=['night'], statlist=['max'])
narr_daynightstat(indir=NARRoutdir, regexpattern='vwnd[.]500[.]nc',
                  outdir=NARRoutdir, dnlist=['day'], statlist=['max'])
narr_daynightstat(indir=NARRoutdir, regexpattern='air[.]sfc[.][0-9]{4}_9x9[.]nc',
                  outdir=NARRoutdir, dnlist=['night'], statlist=['min'])

#Compute 1-, 3- and 6-day maxima, minima, and means
narr_d36stat(indir=NARRoutdir, regexpattern='crain[.][0-9]{4}[_]9x9[.]nc',
             outdir=NARRoutdir, dstat = 'max', multidstat = {6: 'mean'})
narr_d36stat(indir=NARRoutdir, regexpattern='dswrf.*daymin.*[.]nc',
             outdir=NARRoutdir, multidstat = {6: 'max'})
narr_d36stat(indir=NARRdir, regexpattern='air[.]2m.*[.]nc',
             outdir=NARRoutdir, dstat = 'max')
narr_d36stat(indir=NARRoutdir, regexpattern='shum[.]2m[.][0-9]{4}_nightmin[.]nc',
             outdir=NARRoutdir, multidstat = {6: 'mean'})
narr_d36stat(indir=NARRoutdir, regexpattern='lftx4[.][0-9]{4}_nightmin[.]nc',
             outdir=NARRoutdir, multidstat = {3: 'mean'})
narr_d36stat(indir=NARRdir, regexpattern='rhum[.]2m[.][0-9]{4}[.]nc',
             outdir=NARRoutdir, dstat = 'mean')
narr_d36stat(indir=NARRoutdir, regexpattern='air[.]sfc.[0-9]{4}_9x9_nightmin[.]nc',
             outdir=NARRoutdir,multidstat= {6: 'max'})
narr_d36stat(indir=NARRoutdir, regexpattern='vwnd[.]500_daymax[.]nc',
             outdir=NARRoutdir,multidstat= {6: 'max'})
narr_d36stat(indir=NARRoutdir, regexpattern='shum[.]2m[.][0-9]{4}_daymean[.]nc',
             outdir=NARRoutdir, dstat='diff')
narr_d36stat(indir=NARRdir, regexpattern='pres[.]sfc[.][0-9]{4}[.]nc',
             outdir=NARRoutdir, dstat='max')
narr_d36stat(indir=NARRoutdir, regexpattern='rpi[.][0-9]{4}[.]nc',
             outdir=NARRoutdir, dstat='max')
narr_d36stat(indir=NARRoutdir, regexpattern='wspd[.]10m[.][0-9]{4}_daymax[.]nc',
             outdir=NARRoutdir, multidstat = {3: 'max'})
narr_d36stat(indir=NARRoutdir, regexpattern='wspd[.]10m[.][0-9]{4}[.]nc',
             outdir=NARRoutdir, dstat = 'max', multidstat = {3: 'min'})
narr_d36stat(indir=NARRoutdir, regexpattern='wspd[.]10m[.][0-9]{4}[.]nc',
             outdir=NARRoutdir, dstat = 'max', multidstat = {3: 'mean'})
narr_d36stat(indir=NARRdir, regexpattern='lftx4.*[.]nc',
             outdir=NARRoutdir, dstat = 'mean')
narr_d36stat(indir=NARRoutdir, regexpattern='vwnd[.]500[.]nc',
             outdir=NARRoutdir, dstat = 'min')
narr_d36stat(indir=NARRoutdir, regexpattern='hgt[.]850[.]nc',
             outdir=NARRoutdir, dstat='max', multidstat = {6: 'max'})

#-----------------------------------------------------------------------------------------------------------------------
# CREATE POINT SHAPEFILE AND BUFFER
#-----------------------------------------------------------------------------------------------------------------------
sites_chem.groupby('Datum')['UID'].nunique() #Check what datums data are in

#Create geodataframes
sites_wgs84 = pdtogpd_datum(sites_chem, 'WGS84')
sites_nad27 = pdtogpd_datum(sites_chem, 'NAD27')
sites_nad83 = pdtogpd_datum(sites_chem, 'NAD83')

#Get proj4 from NARR netcdfs
templatef = xr.open_mfdataset(glob.glob(os.path.join(NARRdir, varfolder_dic.keys()[0]+'*')))
crsatt = templatef.variables['Lambert_Conformal'].attrs
lcc_proj4 = ("+proj=lcc +lat_1={0} +lat_2={1} +lat_0={2} +lon_0={3} +x_0={4} +y_0={5} +units=m +no_def".
             format(crsatt['standard_parallel'][0],
                    crsatt['standard_parallel'][1],
                    crsatt['latitude_of_projection_origin'],
                    crsatt['longitude_of_central_meridian'],
                    crsatt['false_easting'],
                    crsatt['false_northing']))

#Project them all to same projection as meteorological data for extraction (Lambert Comformal Conic)
sites_gpd_lambers = pd.concat([sites_wgs84.to_crs(crs=lcc_proj4), sites_nad27.to_crs(crs=lcc_proj4), sites_nad83.to_crs(crs=lcc_proj4)])
sites_gpd_lambers['x'], sites_gpd_lambers['y'] = sites_gpd_lambers.geometry.x, sites_gpd_lambers.geometry.y
#Remove outliers with lon=0
sites_gpd_lambers['Longitude'].describe()
sites_gpd_lambers = sites_gpd_lambers[sites_gpd_lambers['Longitude'] != 0]

#Project them all to same projection as NLCD data (and all other data)
aea_proj4 = ("+proj=aea +lat_1={0} +lat_2={1} +lat_0={2} +lon_0={3} +x_0={4} +y_0={5} "
             "+datum=WGS84 +units=m +no_def".
             format(cs_ref.standardParallel1,
                    cs_ref.standardParallel2,
                    cs_ref.latitudeOfOrigin,
                    cs_ref.centralMeridian,
                    cs_ref.falseEasting,
                    cs_ref.falseNorthing))
sites_gpd = pd.concat([sites_wgs84.to_crs(crs=aea_proj4), sites_nad27.to_crs(crs=aea_proj4), sites_nad83.to_crs(crs=aea_proj4)])
sites_gpd['x'], sites_gpd['y'] = sites_gpd.geometry.x, sites_gpd.geometry.y
#Remove outliers with lon=0
sites_gpd['Longitude'].describe()
sites_gpd = sites_gpd[sites_gpd['Longitude'] != 0]

#Create 550 m radius buffers and dissolve them
sites_buf = sites_gpd['geometry'].buffer(distance=550)
sites_bufgpd = gpd.GeoDataFrame(sites_buf, crs=aea_proj4).rename(columns={0:'geometry'}).set_geometry('geometry')

#Dissolve all
sites_bufgpd['diss'] = 1
sites_bufdis = sites_bufgpd.dissolve(by='diss')

#Dissolve only overlapping buffers
sites_bufunion = gpd.GeoDataFrame([polygon for polygon in sites_bufgpd.unary_union], crs=aea_proj4). \
    rename(columns={0:'geometry'}).set_geometry('geometry')

#Output to shapefile
if not os.path.exists(sites_out):
    sites_gpd.to_file(sites_out, driver = 'ESRI Shapefile')
else:
    sites_gpd = gpd.read_file(sites_out)
if not os.path.exists(sites_out_lambers):
    sites_gpd_lambers.to_file(sites_out_lambers, driver = 'ESRI Shapefile')
else:
    sites_gpd_lambers = gpd.read_file(sites_out_lambers)
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
airdat_df_dtypes = {'State Code':np.object, 'County Code':np.object, 'Site Num':np.object, 'Parameter Code': np.object,
                    'Latitude': np.float64, 'Longitude':np.float64, 'Datum':np.object, 'Parameter Name':np.object,
                    'Sample Duration':np.object, 'Pollutant Standard': np.object, 'Date Local': np.object,
                    'Units of Measure': np.object, 'Event Type': np.object, 'Observation Count': int,
                    'Observation Percent':np.float16, 'Arithmetic Mean': np.float64, '1st Max Value': np.float32,
                    '1st Max Hour': np.float32, 'AQI': np.object, 'Method Code': np.object, 'Method Name': np.object,
                    'Local Site Name': np.object, 'Address': np.object, 'State Name': np.object,'County Name': np.object,
                    'City Name': np.object, 'CBSA Name': np.object, 'Date of Last Change': np.object}
airdat_df = pd.read_csv(airdatall, dtype= airdat_df_dtypes, parse_dates=['Date Local', 'Date of Last Change'])

print('Compute unique ID for each monitoring site...')
airdat_df['UID'] = airdat_df['State Code'].astype(str).str.zfill(2) + \
                   airdat_df['County Code'].astype(str) + \
                   airdat_df['Site Num'].astype(str) + \
                   airdat_df["Latitude"].astype(str).str[-3:] + \
                   airdat_df["Longitude"].astype(str).str[-3:]

#-----------------------------------------------------------------------------------------------------------------------
# FORMAT AQI DATA
#-----------------------------------------------------------------------------------------------------------------------
if not os.path.exists(airdat_uniquedf_pickle):
    print('Compute df of unique UID-time records to extract meteorological variables with...')
    subcols = ['UID', 'Latitude', 'Longitude', 'Date Local']
    airdat_unique = pd.unique(airdat_df[subcols[0]].str.cat(airdat_df[subcols[1:]].astype(str), sep='._.')) #Concatenate columns

    airdat_uniquedf = pd.DataFrame(
        [str(i).split('._.') for i in list(airdat_unique)],
        columns = subcols) #Resplit and format to df
    #Faster than airdat_unique.str.split("._.", n = 4, expand = True) ?
    airdat_uniquedf['date'] = airdat_uniquedf['Date Local']

    # Get projected x and y coordinates for each site
    airdat_uniquedfproj = airdat_uniquedf.merge(pd.DataFrame(sites_gpd_lambers.drop(columns='geometry')), on='UID', how='inner')

    # # Expand df to extract 3-hourly meteorological averages
    # airdat_uniquedfexp = airdat_uniquedfproj.reindex(airdat_uniquedfproj.index.repeat(8))
    # airdat_uniquedfexp['datetime_dupli'] = pd.to_datetime(airdat_uniquedfexp['date'], format='%Y-%m-%d') + \
    #                                        pd.to_timedelta(airdat_uniquedfexp.groupby(level=0).cumcount() * 3, unit='h')

    pickle.dump(airdat_uniquedfproj, open(airdat_uniquedf_pickle, "wb"))
else:
    airdat_uniquedfproj = pickle.load(open(airdat_uniquedf_pickle, "rb"))

#-----------------------------------------------------------------------------------------------------------------------
# EXTRACT ALL METEOROLOGICAL VARIABLES FOR EACH STATION-DATE COMBINATION
#-----------------------------------------------------------------------------------------------------------------------
#Make dictionnary of netcdf files to extract to air quality stations
vardict = defaultdict(list)
for netc in getfilelist(NARRoutdir, re.compile('.*(mean|max|min|diff).*[.]nc')):
    vardict[re.sub(r'[0-9]{4}_', '', os.path.splitext(os.path.split(netc)[1])[0])].append(netc)

#Iterate through variables to extract
for var in vardict:
    try:
        print('Processing {}...'.format(var))
        extractCDFtoDF(indf = airdat_uniquedfproj, incdf= vardict[var], indir=NARRdir, varname=var,
                       datecol = 'date', keepcols = ['UID', 'date', 'x', 'y', var], overwrite=False,
                       outfile= os.path.join(rootdir, 'results/daily_SPEC_{}.p'.format(var.replace('.', '_'))))
    except Exception:
        traceback.print_exc()
        print('Skipping...')
        pass

#Collate all air data
try:
    airdat_climmerge = airdat_df
except:
    airdat_climmerge = pd.read_csv(airdatall, dtype= airdat_df_dtypes, parse_dates=['Date Local', 'Date of Last Change'])
    print('Compute unique ID for each monitoring site...')
    airdat_climmerge['UID'] = airdat_climmerge['State Code'].astype(str).str.zfill(2) + \
                              airdat_climmerge['County Code'].astype(str) + \
                              airdat_climmerge['Site Num'].astype(str) + \
                              airdat_climmerge["Latitude"].astype(str).str[-3:] + \
                              airdat_climmerge["Longitude"].astype(str).str[-3:]

plist = [os.path.join(rootdir, 'results/daily_SPEC_{}.p'.format(var.replace('.', '_'))) for var in vardict]
with open(plist[0], "rb") as input_file:
    narrjoin = pickle.load(input_file)
for pname in plist:
    if pname != plist[0]:
        print(pname)
        with open(pname, "rb") as input_file:
            indf = pickle.load(input_file)
        narrjoin = narrjoin.merge(indf.drop(columns=['x', 'y']), on=['UID', 'date'], how='outer')
        print(narrjoin.shape)
narrjoin['Date Local'] = narrjoin['date']
airdat_climmerge = airdat_climmerge.merge(narrjoin, on=['UID', 'Date Local'])

#-----------------------------------------------------------------------------------------------------------------------
# EXTRACT SMOKE DATA FOR EACH STATION-DATE COMBINATION
#-----------------------------------------------------------------------------------------------------------------------
#Intersect smoke polygons with sites
if not arcpy.Exists(sites_smokejointab):
    arcpy.Intersect_analysis([sites_out, smoke1419], sites_smokejoin, join_attributes='ALL')
    arcpy.CopyRows_management(sites_smokejoin, sites_smokejointab)
else:
    print('{} already exists. Skipping...'.format(sites_smokejointab))

#Subset table columns to ease manipulation
if not os.path.exists(sites_smokejoin_subpd):
    sites_smokejoin_gpd = gpd.read_file(AQIgdb, driver='FileGDB', layer=os.path.split(sites_smokejoin)[1])
    sites_smokejoin_sub = sites_smokejoin_gpd[['UID', 'x', 'y', 'ID', 'Density', 'fcname', 'Start', 'End']]
    pickle.dump(sites_smokejoin_sub, open(sites_smokejoin_subpd, "wb"))
else:
    sites_smokejoin_sub = pickle.load(open(sites_smokejoin_subpd, "rb"))

#Format dates
#Establish regex
def smokenamedate_regx(x):
    return re.search(re.compile("(?<=hms_smoke)[0-9]+"), x).group()

def smokehour_regx(x):
    searchout = re.search(re.compile("(?<=\s)[0-9]{4}$"), x)
    if not searchout is None:
        return int(searchout.group())
    else:
        return None

def smokedate_regx(x):
    searchout = re.search(re.compile("^[0-9]{7}(?=\s)"), x)
    if not searchout is None:
        return searchout.group()
    else:
        return None

sites_smokejoin_sub['fcdate'] = pd.to_datetime(sites_smokejoin_sub['fcname'].apply(smokenamedate_regx),
                                               infer_datetime_format=False)
sites_smokejoin_sub['starthour'] = sites_smokejoin_sub['Start'].apply(smokehour_regx)
sites_smokejoin_sub['startdate'] = pd.to_datetime(sites_smokejoin_sub['Start'].apply(smokedate_regx),
                                                  format='%Y%j')
sites_smokejoin_sub['endhour'] = sites_smokejoin_sub['End'].apply(smokehour_regx)
sites_smokejoin_sub['enddate'] = pd.to_datetime(sites_smokejoin_sub['End'].apply(smokedate_regx),
                                                format='%Y%j')

#Check records whose smoke startdate is different than the date of the layer
check = sites_smokejoin_sub[~((sites_smokejoin_sub['fcdate'] == sites_smokejoin_sub['startdate']) |
                              (pd.isnull(sites_smokejoin_sub['startdate'])))]
#Check those whose startdate is after the date of the layer
check[(check.fcdate - check.startdate) < np.timedelta64(-0,'D')][['fcdate', 'startdate', 'enddate', 'Start', 'End']]

#Check records whose smoke enddate is different than the date of the layer
check = sites_smokejoin_sub[~((sites_smokejoin_sub['fcdate'] == sites_smokejoin_sub['enddate']) |
                              (pd.isnull(sites_smokejoin_sub['enddate'])))]
#Check those whose endate is before the date of the layer
check[(check.enddate - check.fcdate) < np.timedelta64(-0,'D')][['fcdate', 'startdate', 'enddate', 'Start', 'End']]

#It seems that these records were assigned to the smoke layer for the wrong day. Only affects <3% of records. Correct them.
sites_smokejoin_sub.loc[((sites_smokejoin_sub.fcdate - sites_smokejoin_sub.startdate) < np.timedelta64(-0,'D')) |
                        ((sites_smokejoin_sub.enddate - sites_smokejoin_sub.fcdate) < np.timedelta64(-0,'D')),
                        'fcdate'] = \
    sites_smokejoin_sub.loc[((sites_smokejoin_sub.fcdate - sites_smokejoin_sub.startdate) < np.timedelta64(-0,'D')) |
                            ((sites_smokejoin_sub.enddate - sites_smokejoin_sub.fcdate) < np.timedelta64(-0,'D')),
                            'startdate']

#For each day, assign proportion of day when smoke was present at that location based on start and end datetime
# smoke density*proportion of hours that day with smoke
# Density indices are 5, 16, and 27 which are the center point of the following smoke categories:
# 0-10 micrograms/m3 for light , 10-21 micrograms/m3 for medium, 21-32 micrograms/m3 for thick
prestart = sites_smokejoin_sub.startdate < sites_smokejoin_sub.fcdate
sites_smokejoin_sub.loc[prestart, 'smokeindex_1d'] = \
    sites_smokejoin_sub.loc[prestart,'Density'].apply(lambda x: re.sub('\s', '', x)).astype(float)* \
    (sites_smokejoin_sub.loc[prestart,'endhour']- 0)/2400.0

postend = sites_smokejoin_sub.enddate > sites_smokejoin_sub.fcdate
sites_smokejoin_sub.loc[postend, 'smokeindex_1d'] = \
    sites_smokejoin_sub.loc[postend, 'Density'].apply(lambda x: re.sub('\s', '', x)).astype(float)* \
    (2400-sites_smokejoin_sub.loc[postend,'starthour'])/2400.0

inday = (sites_smokejoin_sub.enddate == sites_smokejoin_sub.fcdate) & \
        (sites_smokejoin_sub.startdate == sites_smokejoin_sub.fcdate)
sites_smokejoin_sub.loc[inday, 'smokeindex_1d'] = \
    sites_smokejoin_sub.loc[inday, 'Density'].apply(lambda x: re.sub('\s', '', x)).astype(float)* \
    (sites_smokejoin_sub.loc[inday, 'endhour'] - sites_smokejoin_sub.loc[inday, 'starthour'])/2400.0

#Then sum when multiple smoke plumes intersect a point in a day
sites_smoke_stat = sites_smokejoin_sub.groupby(['UID', 'fcdate']).smokeindex_1d.sum().reset_index()

#Get sum of index in the past 3 and 6 days (including current day)
roll3d = sites_smoke_stat.groupby('UID').rolling('3d', on='fcdate').smokeindex_1d.sum().rename('smokeindex_3d')
roll6d = sites_smoke_stat.groupby('UID').rolling('6d', on='fcdate').smokeindex_1d.sum().rename('smokeindex_6d')
sites_smoke_stat = sites_smoke_stat.set_index(['UID', 'fcdate']). \
    join(roll3d, on=['UID', 'fcdate']).join(roll6d, on=['UID', 'fcdate']).reset_index().rename(columns={'fcdate':'date'})

#Merge air quality station + NARR data with smoke variables
airdat_climmerge = pd.merge(airdat_climmerge, sites_smoke_stat, on=['UID', 'date'], how='left')

#Write out data to table (should use feather but some module conflicts and don't want to deal with it)
airdat_climmerge.to_csv(os.path.join(rootdir, 'results/airdat_NARRjoin.csv'))