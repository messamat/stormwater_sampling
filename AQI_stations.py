import arcpy
import mendeleev
import pandas as pd
import itertools
from arcpy.sa import *
import os
import re
import geopandas as gpd
from shapely.geometry import Point

#Set options
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

#Set up paths
rootdir = 'C:/Mathis/ICSL/stormwater'
AQIdir = os.path.join(rootdir, 'data/EPA_AirData_201902')

#Import variables
monitors = pd.read_csv(os.path.join(AQIdir, 'aqs_monitors/aqs_monitors.csv'))
sites = pd.read_csv(os.path.join(AQIdir, 'aqs_sites/aqs_sites.csv'))

#Output variables
outdir = os.path.join(rootdir, 'results/airdata')
if not os.path.isdir(outdir):
    os.mkdir(outdir)
sites_out = os.path.join(outdir, 'airsites.shp')
sites_outbuf = os.path.join(outdir, 'airsites_600buf.shp')
sites_outbufdis = os.path.join(outdir, 'airsites_600bufdis.shp')
sites_outbufunion = os.path.join(outdir, 'airsites_600bufunion.shp')
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
paramsel = ['11','12','14','22','65','82','84','85','86','88']
param_regex = re.compile('|'.join(['(^{}.*)'.format(p) for p in paramsel]))
monitors_chem = monitors.loc[(monitors['Parameter Code'].astype(str).str.contains(param_regex)) &
                             (monitors['Parameter Name'].str.contains(elems_regex)) &
                             ~(monitors['Parameter Name'].str.contains(elemsout_regex)),:]
#pd.unique(monitors_chem['Parameter Name'])

#Pad monitors' codes with 0s to match sites' codes
pd.unique(monitors_chem['County Code'])
monitors_chem['State Code'] = monitors_chem['State Code'].astype(str).str.pad(2, 'left', '0')

#Compute Unique Site Identifier and select sites that record chemical concentrations
monitors_chem['UID'] = monitors_chem['State Code'] + \
                       monitors_chem['County Code'].astype(str) + \
                       monitors_chem['Site Number'].astype(str)
sites['UID'] = sites['State Code'].astype(str).astype(str) + \
               sites['County Code'].astype(str) + \
               sites['Site Number'].astype(str)
sites_chem = sites[sites['UID'].isin(pd.unique(monitors_chem['UID']))]
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

#Create 600 m radius buffers and dissolve them
sites_buf = sites_gpd['geometry'].buffer(distance=600)
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
