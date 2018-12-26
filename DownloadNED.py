#Selected the Puget Sound polygon from https://viewer.nationalmap.gov/basic/#productGroupSearch
#Added all 1/9 arc-second NED tiles to cart, and exported URLs to text file

import os
import urllib
import zipfile
import arcpy
import re

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
    except:
        arcpy.MosaicToNewRaster_management(arcpy.ListRasters('grd*'), outdir, outras, number_of_bands=1,
                                           pixel_type= '16_BIT_SIGNED')

    if arcpy.Exists(os.path.join(rootdir, 'results/NED_PS')):
        os.remove('zipped')

#Download and mosaic NED 1/9 arc-second
download_unzip_ned(tempdir = os.path.join(rootdir, 'data/NED19'), URLlist = 'cartExport_20181219_165432.txt',
                   outdir = os.path.join(rootdir, 'results'), outras = 'NED19_PS')
#Download and mosaic NED 1/3 arc-second
download_unzip_ned(tempdir = os.path.join(rootdir, 'data/NED13'), URLlist = 'cartExport_20181220_114439.txt',
                   outdir = os.path.join('D:/Processing'), outras = 'NED13_PS', download=None, extract=None)