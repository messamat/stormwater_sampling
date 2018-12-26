#Selected the Puget Sound polygon from https://viewer.nationalmap.gov/basic/#productGroupSearch
#Added all 1/9 arc-second NED tiles to cart, and exported URLs to text file

import os
import urllib
import zipfile
import arcpy
import re

UTM10 = arcpy.SpatialReference(26910)
arcpy.CheckOutExtension('Spatial')

rootdir = 'C:/Mathis/ICSL/stormwater'
def download_unzip_ned(tempdir, URLlist, outdir, outras, download=None, extract=None):
    os.chdir(tempdir)

    if not os.path.exists('zipped'):
        os.mkdir('zipped')

    if download:
        with open(URLlist, 'r') as tab:
            for l in tab.readlines():
                outzip = os.path.join(os.getcwd(), 'zipped', os.path.splitext(os.path.split(l)[1])[0]) + '.zip'
                if not os.path.exists(outzip):
                    print('Downloading {}'.format(l))
                    urllib.urlretrieve(l, os.path.join(outzip))

    if extract:
        for root, dir, files, in os.walk(os.path.join(os.getcwd(), 'zipped')):
            for f in files:
                fraw = os.path.splitext(f)[0]
                ras_regex = re.compile('.*{}.*(?![.]zip)$'.format(f))
                if not any([ras_regex.search(i) for i in os.listdir(tempdir)]):
                    print('Extracting {}'.format(f))
                    try:
                        zip = zipfile.ZipFile(os.path.join(root, f))
                        zip.extractall()
                    except:
                        print('Zipfile still loading')

    arcpy.env.workspace = os.getcwd()
    print('Mosaicking rasters')
    try:
        arcpy.MosaicToNewRaster_management(arcpy.ListRasters('*.img'), outdir, outras, number_of_bands= 1)
    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages(2))
        raslist = filter(re.compile('grd.*(?<!jpg)(?<!xml)$').search, os.listdir(os.getcwd()))
        arcpy.MosaicToNewRaster_management(raslist, outdir, outras, number_of_bands=1, pixel_type= '16_BIT_SIGNED')

    if arcpy.Exists(outras):
        os.remove('zipped')


arcpy.env.workspace = os.path.join(rootdir, 'results')
#Download and mosaic NED 1/9 arc-second
download_unzip_ned(tempdir = os.path.join(rootdir, 'data/NED19'), URLlist = 'cartExport_20181221_075329.txt',
                   outdir = 'D:/Processing', outras = 'ned19_ps', download = None, extract = None)
if not arcpy.Exists('D:/Processing/ned19_psproj'):
    arcpy.ProjectRaster_management('D:/Processing/ned19_ps', 'D:/Processing/ned19_psproj', UTM10, resampling_type ='BILINEAR')
if not arcpy.Exists('D:/Processing/slope19_ps'):
    slope = arcpy.sa.Slope('D:/Processing/ned19_psproj', 'DEGREE')
    slope.save('D:/Processing/slope19_ps')


#Download and mosaic NED 1/3 arc-second
if not arcpy.Exists('ned_ps'):
    download_unzip_ned(tempdir = os.path.join(rootdir, 'data/NED13'), URLlist = 'cartExport_20181220_114439.txt',
                       outdir = os.path.join(rootdir, 'results'), outras = 'ned_ps', download=None, extract=None)
if not arcpy.Exists('ned_psproj'):
    arcpy.ProjectRaster_management('ned_ps', 'ned_psproj', UTM10, resampling_type ='BILINEAR')
if not arcpy.Exists('slope_ps'):
    slope = arcpy.sa.Slope('ned_psproj', 'DEGREE')
    slope.save('slove_ps')