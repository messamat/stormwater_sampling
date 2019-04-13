'''
Author: Mathis Messager
Contact info: messamat@uw.edu
Creation date: December 2018

Required Arguments:
    - fc (feature class in gdb): line feature class to be divided into non-overlapping sets
    - tolerance (distance in fc coordinate system's linear unit): size of buffer around lines to determine intersection
Optional Argument:
    - keep (True/False): whether to keep (True, default) intermediate outputs or not (False)

Description: divide a line feature class into non-overlapping sets (e.g. for rasterizing and zonal statistics) by
            creating a buffer around the lines, intersecting the buffers, and then creating a new field called 'expl'
            containing which group each line belongs to. The output can be used to rasterize the feature class by
            iteratively selecting by 'expl' and mosaicking with a SUM or MEAN function.

Inspired from https://gis.stackexchange.com/questions/32217/exploding-overlapping-to-new-non-overlapping-polygons
'''

import arcpy
import os
from collections import defaultdict

arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

tolerance = (2.0**0.5)*float(res.getOutput(0))/2
fc = NTMsplitdiss
arcpy.env.workspace = os.path.join(rootdir, 'results/NTM.gdb')
fcbub = fc + 'buf'

def ExplodeOverlappingLines(fc, tolerance, keep=True):
    print('Buffering lines...')
    idName = "ORIG_FID"
    fcbuf = arcpy.Buffer_analysis(fc, fc+'buf', tolerance, line_side='FULL', line_end_type='FLAT')

    print('Intersecting buffers...')
    intersect = arcpy.Intersect_analysis(fcbuf,'intersect')

    print('Creating dictionary of overlaps...')
    #Find identical shapes and put them together in a dictionary
    segIDs = defaultdict(list)
    with arcpy.da.SearchCursor(intersect, ['Shape@WKT', idName]) as cursor:
        x=0
        for row in cursor:
            if x%100000 == 0:
                print(x)
            segIDs[row[0]].append(row[1])
            x+=1

    #Build dictionary of all buffers overlapping each buffer
    segIDs2 = defaultdict(list)
    for v in segIDs.values():
        for segID in v:
            segIDs2[segID].extend([k for k in v if k != segID and k not in segIDs2[segID]])

    print('Assigning lines to non-overlapping sets...')
    grpdict = {}
    # Mark all non-overlapping one to group 1
    for row in arcpy.da.SearchCursor(fcbuf, [idName]):
        if row[0] in segIDs2:
            grpdict[row[0]] = None
        else:
            grpdict[row[0]] = 1

    i = 1
    while None in grpdict.values():
        print(i)
        ovList = []
        for kv in sorted(segIDs2.items(), key=lambda kv: (len(kv[1]), kv[0])):
            if grpdict[kv[0]] is None:
                if kv[0] not in ovList:
                    grpdict[kv[0]] = i
                    ovList.extend(kv[1])
        i += 1

    print('Writing out results to "expl" field in...'.format(fc))
    arcpy.AddField_management(fc, 'expl', "SHORT")
    with arcpy.da.UpdateCursor(fc,
                               [arcpy.Describe(fc).OIDfieldName, 'expl']) as cursor:
        for row in cursor:
            row[1] = grpdict[row[0]]
            cursor.updateRow(row)

    if keep == False:
        print('Deleting intermediate outputs...')
        for fc in ['intersect', "explFindID"]:
            arcpy.Delete_management(fc)