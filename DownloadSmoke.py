import arcpy
import os
import re
from datetime import date, timedelta
import glob
import requests

#Custom functions
from Download_gist import *

#Set up paths
rootdir = 'D:/Mathis/ICSL/stormwater'
AQIdir = os.path.join(rootdir, 'data/EPA_AirData_201902')

AQIgdb = os.path.join(rootdir, 'results/airdata/AQI.gdb')
#Create gdb for analysis
if arcpy.Exists(AQIgdb):
    print('Geodatabase already exists')
else:
    arcpy.CreateFileGDB_management(os.path.join(rootdir,'results/airdata'), 'AQI.gdb')

NLCD_imp = os.path.join(rootdir, 'data/NLCD_2016_Impervious_L48_20190405.img') #Based on 2016 dara
cs_ref = arcpy.Describe(NLCD_imp).SpatialReference

smokedir = os.path.join(rootdir, 'data/HMS_ARCHIVE_SMOKE')
if not arcpy.Exists(smokedir):
    print('Creating directory {}...'.format(smokedir))
    os.mkdir(smokedir)
smoke1419 = os.path.join(AQIgdb, 'smoke1419')

# -----------------------------------------------------------------------------------------------------------------------
# DOWNLOAD AND MERGE SMOKE DATA
# -----------------------------------------------------------------------------------------------------------------------
smokedates = pd.date_range(start='2014-01-01', end=date.today() - timedelta(days=1)).strftime('%Y%m%d')
failedlist = []
for d in smokedates:
    try:
        for ext in ['shp', 'dbf', 'shx']:
            # Try downloading as .gzip
            smokeurl = "https://satepsanone.nesdis.noaa.gov/pub/volcano/" \
                       "FIRE/HMS_ARCHIVE/{0}/GIS/SMOKE/hms_smoke{1}.{2}".format(d[0:4], d, ext)
            dlfile(url="{}.gz".format(smokeurl), outpath=smokedir)

            # If not gz file (stopped zipping data after summer 2018)
            if not os.path.exists(os.path.join(smokedir, os.path.split(smokeurl)[1])):
                response = requests.get(smokeurl, stream=True)
            if response.status_code == 200:
                print "downloading " + smokeurl
            with open(os.path.join(smokedir, os.path.split(smokeurl)[1]), 'wb') as f:
                f.write(response.raw.read())
    except Exception:
        traceback.print_exc()
        failedlist.append(smokeurl)
        pass

# Get rid of empty shapefiles
# arcpy.env.workspace = smokedir; smokedatlist = arcpy.ListFeatureClasses() crashes with exit code
for t in glob.glob(os.path.join(smokedir, '*.shp')):
    if arcpy.Describe(t).shapeType == u'Null':
        print('Delete {}, null feature class...'.format(t))
        arcpy.Delete_management(t)  # Obnoxious bugging with locks when deleting. Delete manually if needed.
    else:
        nameformat = os.path.splitext(os.path.split(t)[1])[0]
        print('Add feature class name to fcname field for {}...'.format(nameformat))
        arcpy.AddField_management(t, 'fcname', 'TEXT', field_length=20)
        with arcpy.da.UpdateCursor(t, ['fcname']) as cursor:
            for row in cursor:
                row[0] = nameformat
                cursor.updateRow(row)

# smokefdic = defaultdict(list)
# for t in glob.glob(os.path.join(smokedir, '*shp')):
#     ftypes_proj = [[f.name, f.type, f.length, f.precision, arcpy.Describe(t).SpatialReference.name] for f in arcpy.ListFields(t)]
#     smokefdic[str(ftypes_proj)].append(t)


# Merge all smoke polygons for each year separately, making sure that field lengths are homogeneous
# Do not merge all years at the same time as default Merge crashes and fieldMappings.addTable(t) slows to a crawl
# after a few hundreds feature classes. Append_management is equally slow
for yr in range(2014, 2020):
    outsmokeyr = os.path.join(AQIgdb, 'smoke{}'.format(yr))
    fclist = glob.glob(os.path.join(smokedir, '*{}*.shp'.format(yr)))
    if not arcpy.Exists(outsmokeyr):
        fieldMappings = arcpy.FieldMappings()
        print('Adding tables to field mapping...')
        for fc in fclist:
            fieldMappings.addTable(fc)
        for fm in fieldMappings.fieldMappings:
            if fm.outputField.name == 'Start':
                if fm.outputField.length < 12:
                    fmstart = fm
                    outstartf = fmstart.outputField
                    outstartf.length = 12
                    fmstart.outputField = outstartf
                    fieldMappings.removeFieldMap(fieldMappings.findFieldMapIndex('Start'))
                    fieldMappings.addFieldMap(fmstart)

            if fm.outputField.name == 'End':
                if fm.outputField.length < 12:
                    fmend = fm
                    outendf = fmend.outputField
                    outendf.length = 12
                    fmend.outputField = outendf
                    fieldMappings.removeFieldMap(fieldMappings.findFieldMapIndex('End'))
                    fieldMappings.addFieldMap(fmend)
        print('Merging all daily smoke datasets for year {}...'.format(yr))
        arcpy.Merge_management(fclist, outsmokeyr, fieldMappings)
    else:
        print('{} already exists, skipping...'.format(outsmokeyr))

# Merge all years together
arcpy.Merge_management([os.path.join(AQIgdb, 'smoke{}'.format(yr)) for yr in range(2014, 2020)], smoke1419)

# Project to Albers Equal Area (same as rest of the analysis)
arcpy.DefineProjection_management(smoke1419, coor_system=arcpy.SpatialReference(4326))
arcpy.Project_management(smoke1419, '{}_aea'.format(smoke1419), out_coor_system=cs_ref)