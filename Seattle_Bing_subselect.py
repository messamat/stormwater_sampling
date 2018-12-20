import re
import arcpy
from os import *
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

#Set up paths
rootdir = 'C:/Mathis/ICSL/stormwater/'
res= path.join(rootdir, 'results/bing/')
arcpy.env.workspace=res
mlclist = [filenames for (dirpath, dirnames, filenames) in walk(res)][0]
regex= re.compile(".+mlc\.tif$")
clasfiles = [m.group(0) for l in mlclist for m in [regex.search(l)] if m] #Look into this
kernel = NbrWeight(path.join(rootdir, 'results/logkernel100.txt')) #UPDATE