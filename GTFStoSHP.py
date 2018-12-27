'''
Author: Mathis Messager
Contact info: messamat@uw.edu
Creation date: December 2018

Required Arguments:
    - gtfs_dir (directory or .zip file): directory that contains GTFS* text file set or zip file. If zip file, it will
            be extracted to the directory where the .zip file is located
    - out_gdb (path of geodatabase where outputs will be written, must end in .gdb): if does not exist, it will be created
    - out_fc (name): rootname of feature class that will contain the transit shapes with the weekly number of buses that run it

Description: Takes a GTFS file set (see https://developers.google.com/transit/gtfs/reference/#agencytxt for description),
    convert shapes (unique bus line-daily schedule combination, e.g. 6:47 am MON-FRI Line 71) from text file to feature
    class and compute for each shape the number of buses that run on it every week.

*GTFS: General Transit Feed Specification
Inspired from http://www.stevencanplan.com/2010/10/how-to-convert-gtfs-to-shapefiles-and-kml/
'''

import arcpy
import os
import zipfile

arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

def GTFStoSHPweeklynumber(gtfs_dir, out_gdb, out_fc):
    if not arcpy.Exists(out_gdb):
        print('Creating {}...'.format(out_gdb))
        arcpy.CreateFileGDB_management(os.path.dirname(out_gdb), os.path.basename(out_gdb))
    else:
        print('{} already exists, using existing gdb...'.format(out_gdb))
    arcpy.env.workspace = out_gdb

    #

    if os.path.splitext(gtfs_dir)[1] == '.zip':
        zip = zipfile.ZipFile(gtfs_dir)
        zip.extractall(os.path.split(gtfs_dir)[0])
        gtfs_dir = os.path.split(gtfs_dir)[0]

    #Create shapefile from shapes.txt
    print('Creating point feature class...')
    newfc = arcpy.CreateFeatureclass_management(out_gdb, out_fc, "Point")
    arcpy.DefineProjection_management(newfc , 4326)
    
    txt = open(os.path.join(gtfs_dir, 'shapes.txt'), 'r')
    inputF = txt.readlines()
    
    arcpy.AddField_management(newfc, 'shape_id', 'LONG')
    arcpy.AddField_management(newfc, 'lat', 'FLOAT')
    arcpy.AddField_management(newfc, 'lon', 'FLOAT')
    arcpy.AddField_management(newfc, 'seq', 'LONG')
    arcpy.AddField_management(newfc, 'dist', 'DOUBLE')
    
    with arcpy.da.InsertCursor(newfc, ["shape_id","SHAPE@XY", "lat", "lon", "seq", "dist"]) as cursor:
        for line in inputF[1:]:
            id, ycoor, xcoor, seq, dist = line.split(",")
            xy = (float(xcoor), float(ycoor))
            newRow = (id, xy, ycoor, xcoor, seq, dist)
            cursor.insertRow(newRow)
    
    #Convert to polyline
    print('Converting points to lines...')
    newfc_lines = out_fc + '_lines'
    arcpy.PointsToLine_management(newfc, newfc_lines, Line_Field= 'shape_id')
    
    #Import GTFS text files to out_gdb
    print('Importing GTFS tables into gdb...')
    arcpy.CopyRows_management(os.path.join(gtfs_dir,'routes.txt'), 'routes')
    arcpy.CopyRows_management(os.path.join(gtfs_dir,'trips.txt'), 'trips')
    arcpy.CopyRows_management(os.path.join(gtfs_dir,'calendar.txt'), 'calendar')
    #arcpy.CopyRows_management(os.path.join(gtfs_dir,'stop_times.txt'), 'stop_times')
    
    #Each route has multiple trips
    #Eacj trip is a unique bus line-time-week day combination (e.g. 6:47 am MON-FRI Line 71)
    #Compute number of days/week for each service_id in calendar
    print('Summaryzing data...')
    arcpy.AddField_management('calendar', 'daynum', 'SHORT')
    cal_fields = ['daynum', 'monday', 'tuesday', 'wednesday', 'thursday', 'friday', 'saturday', 'sunday']
    with arcpy.da.UpdateCursor('calendar', cal_fields) as cursor:
        for row in cursor:
            row[0] = sum(row[1:])
            cursor.updateRow(row)
    
    #Join calendar to trips
    arcpy.MakeTableView_management('trips', 'trips_view')
    arcpy.AddJoin_management('trips_view', 'service_id', 'calendar', 'service_id', join_type='KEEP_ALL')
    arcpy.AddJoin_management('trips_view', 'route_id', 'routes', 'route_id', join_type='KEEP_ALL')
    arcpy.CopyRows_management('trips_view', 'trips_routes_calendar')
    
    #Summary trips by route_id & shape_id (here add
    arcpy.Statistics_analysis('trips_routes_calendar', 'trips_routes_count', statistics_fields=[['daynum', 'SUM']],
                              case_field = ['shape_id', 'route_id', 'route_short_name', 'route_desc', 'route_type'])
    #Join trip statistics to lines
    print('Joining and exporting data...')
    arcpy.MakeFeatureLayer_management(newfc_lines, 'shapes_lyr')
    arcpy.AddJoin_management('shapes_lyr', 'shape_id', 'trips_routes_count', 'shape_id')
    arcpy.CopyFeatures_management('shapes_lyr', 'kingtransit_routes')

    print('Deleting intermediate outputs...')
    for inter_lyr in ['routes', 'trips', 'calendar', newfc, newfc_lines]:
        arcpy.Delete_management(inter_lyr)

# rootdir =
# GTFStoSHPweeklynumber(gtfs_dir = os.path.join(rootdir, 'data/King_201812'),
#                       out_gdb = os.path.join(rootdir, 'results/transit.gdb'),
#                       out_fc = 'kingtransit')