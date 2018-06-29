import arcpy
from arcpy.sa import *
from os import *
import numpy as np
from collections import defaultdict

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput=True

#Set up paths
rootdir = 'C:/Mathis/ICSL/stormwater/'
arcpy.env.workspace = path.join(rootdir, 'results')
#New variables
sgdb= path.join(rootdir,'results/Seattle_sampling.gdb')
roadstraffic_avg = path.join(sgdb,'Seattle_roadstraffic_AADT')
PSOSM_all='PSwtshd_OSMroads_all.shp'
UTM10 = arcpy.SpatialReference(26910)

