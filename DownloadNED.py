#Selected the Puget Sound polygon from https://viewer.nationalmap.gov/basic/#productGroupSearch
#Added all 1/9 arc-second NED tiles to cart, and exported URLs to text file

import os
import urllib
import zipfile
import arcpy
import re

arcpy.CheckOutExtension('Spatial')

rootdir = 'D:/Mathis/ICSL/stormwater'
NLCD_imp = os.path.join(rootdir, 'data/NLCD_2016_Impervious_L48_20190405.img') #Based on 2016 dara
ref_cs = arcpy.Describe(NLCD_imp).SpatialReference
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
if not arcpy.Exists('ned19_ps'):
    download_unzip_ned(tempdir = os.path.join(rootdir, 'data/NED19'), URLlist = 'cartExport_20181221_075329.txt',
                       outdir = os.path.join(rootdir, 'results'), outras = 'ned19_ps', download = True, extract = True)
if not arcpy.Exists('ned19_psproj'):
    arcpy.ProjectRaster_management('ned19_ps', 'ned19_psproj', ref_cs, resampling_type ='BILINEAR')
# if not arcpy.Exists('slope19_ps'):
#     slope = arcpy.sa.Slope('ned19_psproj', 'DEGREE')
#     slope.save('slope19_ps')


#Download and mosaic NED 1/3 arc-second
if not arcpy.Exists('ned13_ps'):
    download_unzip_ned(tempdir = os.path.join(rootdir, 'data/NED13'), URLlist = 'cartExport_20181220_114439.txt',
                       outdir = os.path.join(rootdir, 'results'), outras = 'ned13_ps', download=True, extract=True)
if not arcpy.Exists('ned13_psproj'):
    arcpy.ProjectRaster_management('ned13_ps', 'ned13_psproj', ref_cs, resampling_type ='BILINEAR')
# if not arcpy.Exists('slope13_ps'):
#     slope = arcpy.sa.Slope('ned13_psproj', 'DEGREE')
#     slope.save('slove13_ps')