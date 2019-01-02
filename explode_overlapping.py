'''
Inspired from https://gis.stackexchange.com/questions/32217/exploding-overlapping-to-new-non-overlapping-polygons
'''

import arcpy
import os
from collections import defaultdict

arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

def explodeOverlaps(fc, idName, tolerance):
    idName = "ORIG_FID"
    fcbuf = arcpy.Buffer_analysis(fc, fc+'buf', tolerance, line_side='FULL', line_end_type='FLAT')
    intersect = arcpy.Intersect_analysis(fcbuf,'intersect')
    findID = arcpy.FindIdentical_management(intersect,"explFindID","Shape")
    arcpy.MakeFeatureLayer_management('intersect',"intlyr")
    arcpy.AddJoin_management("intlyr",arcpy.Describe("intlyr").OIDfieldName,"explFindID","IN_FID","KEEP_ALL")
    featseqName = "explFindID.FEAT_SEQ"
    idNewName = "intersect." + idName
    segIDs = defaultdict(list)
    segIDs2 = defaultdict(list)

    for row in arcpy.da.SearchCursor("intlyr", [featseqName, idNewName]):
        segIDs[row[0]].append(row[1])
    for v in segIDs.values():
        for segID in v:
            segIDs2[segID].extend([k for k in v if k != segID and k not in segIDs2[segID]])

    arcpy.RemoveJoin_management("intlyr", arcpy.Describe("explFindID").name)

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

    arcpy.AddField_management(fc, 'expl', "SHORT")
    with arcpy.da.UpdateCursor(fc,
                               [arcpy.Describe(fc).OIDfieldName, 'expl']) as cursor:
        for row in cursor:
            row[1] = grpdict[row[0]]
            cursor.updateRow(row)