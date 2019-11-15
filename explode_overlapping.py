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

import cPickle as pickle
# with open(os.path.join(rootdir, 'results/file.pkl'), 'rb') as f:
#     segIDs2 = pickle.load(f)

# f = open(os.path.join(rootdir, 'results/file.pkl'),"wb")
# pickle.dump(segIDs2,f)
# f.close()

def ExplodeOverlappingLines(fc, tolerance, keep=True, overwrite = False):
    idName = "ORIG_FID"

    if arcpy.Exists(fc+'buf') and overwrite==False:
        print('Line buffers already exist... skipping')
    else:
        print('Buffering lines...')
        arcpy.Buffer_analysis(fc, fc+'buf', tolerance, line_side='FULL', line_end_type='FLAT')

    if arcpy.Exists(fc+'intersect') and overwrite==False:
        print('Line buffer intersections already exist... skipping')
    else:
        print('Intersecting buffers...')
        arcpy.Intersect_analysis(fc+'buf', fc+'intersect')
        arcpy.RepairGeometry_management(fc+'intersect', delete_null='DELETE_NULL')

    print('Creating dictionary of overlaps...')
    #Find identical shapes and put them together in a dictionary, unique shapes will only have one value
    segIDs = defaultdict(list)
    #First delete corrupted shape
    keeplooping = True
    with arcpy.da.UpdateCursor(fc + 'intersect', ['SHAPE@WKT', 'OID@']) as cursor:
        while keeplooping == True:
            try:
                for row in cursor:
                    pass
                keeplooping = False
            except:
                print('Delete feature with OBJECTID {}'.format(row[1]))
                cursor.deleteRow()

    #Actually run through dataset
    with arcpy.da.UpdateCursor(fc+'intersect', ['SHAPE@WKT', idName]) as cursor:
        x=0
        for row in cursor:
            if x%100000 == 0:
                print('Processed {} records for duplicate shapes...'.format(x))
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
    for row in arcpy.da.SearchCursor(fc+'buf', [idName]):
        if row[0] in segIDs2:
            grpdict[row[0]] = None
        else:
            grpdict[row[0]] = 1

    segIDs2sort = sorted(segIDs2.items(), key=lambda x: (len(x[1]), x[0])) #Sort dictionary by number of overlapping features then by keys
    i = 2
    while None in grpdict.values(): #As long as there remain features not assigned to a group
        print(i)
        ovset = set()  # list of all features overlapping features within current group
        x = 0
        s_update = ovset.update
        for rec in segIDs2sort:
            if grpdict[rec[0]] is None: #If feature has not been assigned a group
                if rec[0] not in ovset: #If does not overlap with a feature in that group
                    print(rec[0])
                    grpdict[rec[0]] = i  # Assign current group to feature
                    s_update(rec[1])  # Add all overlapping feature to ovList
        i += 1 #Iterate to the next group

    print('Writing out results to "expl" field in...{}'.format(fc))
    arcpy.AddField_management(fc, 'expl', "SHORT")
    with arcpy.da.UpdateCursor(fc,
                               [arcpy.Describe(fc).OIDfieldName, 'expl']) as cursor:
        for row in cursor:
            if row[0] in grpdict:
                row[1] = grpdict[row[0]]
                cursor.updateRow(row)

    if keep == False:
        print('Deleting intermediate outputs...')
        for fc in [fc+'intersect', "explFindID"]:
            arcpy.Delete_management(fc)