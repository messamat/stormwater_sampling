import arcpy
import os
import zipfile
import fnmatch

#Set options
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

#Set up paths
rootdir = 'C:/Mathis/ICSL/stormwater'
onlinedir = os.path.join(rootdir, 'results/map_forluwam')
gdb = os.path.join(rootdir,'results/Seattle_sampling.gdb')
PSgdb=os.path.join(rootdir,'results/PSOSM.gdb')
XRFsites_proj = os.path.join(gdb, 'XRFsites_proj')
trees = os.path.join(rootdir, 'data/CitySeattle_20180601/Trees/Trees.shp')

sites_out = os.path.join(onlinedir, 'XRFsites_proj.shp')
trees_out = os.path.join(onlinedir, 'trees.shp')
arcpy.CopyFeatures_management(XRFsites_proj, sites_out)
arcpy.CopyFeatures_management(trees, trees_out)

path = trees_out

def zipshp(path):
    outdir = os.path.split(path)[0]
    keyw = os.path.splitext(os.path.split(path)[1])[0]
    zip = zipfile.ZipFile(os.path.join(outdir, keyw+".zip"), "w", zipfile.ZIP_DEFLATED)
    filelist = [file for file in [os.path.join(outdir, keyw + ext)
                  for ext in ['.cpg', '.dbf', '.prj', '.sbn', '.sbx', '.shp', '.shx']] if os.path.isfile(file)]
    for f in filelist:
        zip.write(f, os.path.basename(f))
    zip.close()

zipshp(trees_out)
zipshp(sites_out)

arcpy.CreateMapTilePackage_management()