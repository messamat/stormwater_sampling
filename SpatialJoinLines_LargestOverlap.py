'''
Author: Mathis Messager
Contact info: messamat@uw.edu
Creation date: December 2018

Required Arguments:
        Target Features (Feature Layer)
        Join Features (Feature Layer)
        outgdb (File Geodatabase)
        Output Feature Class (Feature Class)
        bufsize (string input for arcpy buffer tool, must include unit)
        fields_select (list of fields from join features to add to target features table)
Optional Arguments:
        Keep All (True|False)
        Spatial Relationship (String)

Description: Joins attributes from one line feature class to another based on the spatial
      relationship between the two. The target features and the
      attributes from the join features are written to the output feature
      class. Both line feature classes are first split at self intersections.
      Then, the join feature line whose buffer of user-determined size shares the
      most overlapping area with a target feature line will be joined.

Inspired from Spatial Join Largest Overlap function from ESRI:
http://www.arcgis.com/home/item.html?id=e9cccd343bf84916bda1910c31e5eab2
'''

# Import system modules
import arcpy
import os

arcpy.env.overwriteOutput = True

def SpatialJoinLines_LargestOverlap(target_features, join_features, outgdb, out_fc, bufsize, keep_all, fields_select):
    arcpy.env.extent = target_features
    arcpy.env.workspace = outgdb

    #Split target and join lines at intersections
    arcpy.FeatureToLine_management(target_features, 'target_split')
    arcpy.FeatureToLine_management(join_features, 'joinfeat_split')

    #Bufferize both datasets
    arcpy.Buffer_analysis('target_split', 'target_buf', bufsize, method='GEODESIC')
    arcpy.Buffer_analysis('joinfeat_split', 'joinfeat_buf', bufsize, method='GEODESIC')
    #Get buffer area for target feature
    arcpy.AddGeometryAttributes_management('target_buf', 'AREA_GEODESIC', Area_Unit='SQUARE_METERS')

    #Spatial join with largest overlap
    # Calculate intersection between Target Feature and Join Features
    arcpy.Intersect_analysis(['target_buf', 'joinfeat_buf'], 'lines_intersect', join_attributes='ALL')
    arcpy.AlterField_management('lines_intersect','AREA_GEO','AREA_targetbuf', 'AREA_targetbuf')
    arcpy.AddGeometryAttributes_management('lines_intersect', 'AREA_GEODESIC', Area_Unit='SQUARE_METERS')
    arcpy.AlterField_management('lines_intersect','AREA_GEO','AREA_inters', 'AREA_inters')

    #Dissolve to sum intersecting area over
    arcpy.Statistics_analysis('lines_intersect', 'lines_intersect_stats',
                              statistics_fields = [['AREA_inters', 'SUM'], ['AREA_targetbuf', 'FIRST']],
                              case_field = ['FID_joinfeat_buf', 'FID_target_buf'])


    [f.name for f in arcpy.ListFields('lines_intersect_stats')]
    # Find which Join Feature has the largest overlap with each Target Feature
    fields = ['FID_target_buf', 'FID_joinfeat_buf', "SUM_AREA_inters", "FIRST_AREA_targetbuf"]
    overlap_dict = {}
    with arcpy.da.SearchCursor('lines_intersect_stats', fields) as scur:
        for row in scur:
            try:
                if row[2] > overlap_dict[row[0]][1]:
                    overlap_dict[row[0]] = [row[1], row[2], row[3]]
            except:
                overlap_dict[row[0]] = [row[1], row[2], row[3]]

    # Copy the target features and write the largest overlap join feature ID to each record
    # Set up all fields from the target features + ORIG_FID
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(target_features)
    fieldmap = arcpy.FieldMap()
    fieldmap.addInputField(target_features, arcpy.Describe(target_features).OIDFieldName)
    fld = fieldmap.outputField
    fld.type, fld.name, fld.aliasName = "LONG", "ORIG_FID", "ORIG_FID"
    fieldmap.outputField = fld
    fieldmappings.addFieldMap(fieldmap)
    # Perform the copy
    arcpy.conversion.FeatureClassToFeatureClass('target_split',
                                                os.path.dirname(out_fc),
                                                os.path.basename(out_fc), "", fieldmappings)
    # Add a new field JOIN_FID to contain the fid of the join feature with the largest overlap
    arcpy.AddField_management(out_fc, "JOIN_FID", "LONG")
    arcpy.AddField_management(out_fc, "AREA_inters", "DOUBLE")
    arcpy.AddField_management(out_fc, "AREA_targetbuf", "DOUBLE")
    arcpy.AddField_management(out_fc, "intersper", "DOUBLE")
    # Calculate the JOIN_FID field
    with arcpy.da.UpdateCursor(out_fc, ["ORIG_FID", "JOIN_FID", "AREA_inters", "AREA_targetbuf", "intersper"]) as ucur:
        for row in ucur:
            try:
                row[1] = overlap_dict[row[0]][0]
                row[2] = overlap_dict[row[0]][1]
                row[3] = overlap_dict[row[0]][2]
                row[4] = overlap_dict[row[0]][1]/ overlap_dict[row[0]][2]#intersper = AREA_inters/AREA_targetbuf
                ucur.updateRow(row)
            except:
                if not keep_all:
                    ucur.deleteRow()

    # Join all attributes from the join features to the output
    joinfields = [x.name for x in arcpy.ListFields('joinfeat_split') if not x.required and x.name in fields_select]
    arcpy.JoinField_management(out_fc, "JOIN_FID", 'joinfeat_split', arcpy.Describe('joinfeat_split').OIDFieldName, joinfields)

    for outlyr in ['joinfeat_split', 'target_slip','target_buf',
                   'joinfeat_buf', 'lines_intersect', 'lines_intersect_stats']:
        print('Deleting {}'.format(outlyr))
        arcpy.Delete_management(outlyr)