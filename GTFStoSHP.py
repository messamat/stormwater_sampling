'''
Author: Mathis Messager
Contact info: messamat@uw.edu
Creation date: December 2018

Required Arguments:
    - gtfs_dir (directory or .zip file): directory that contains GTFS* text file set or zip file. If zip file, it will
            be extracted to the directory where the .zip file is located
    - out_gdb (path of geodatabase where outputs will be written, must end in .gdb): if does not exist, it will be created
    - out_fc (name): rootname of feature class that will contain the transit shapes with the weekly number of trips on it

Optional arguments:
    - current (True/False): whether to only take in account current trips (True) or also future trips (False)
                            based on calendar.txt start_date field and current date (datetime.datetime.now())
    - keep (True/False): whether to keep (True) intermediate outputs or not (False)

Description: Takes a GTFS file set (see https://developers.google.com/transit/gtfs/reference/#agencytxt for description),
    convert shapes (unique e.g. bus line-daily schedule combination, e.g. 6:47 am MON-FRI Line 71) from text file to feature
    class and compute for each shape the number of trips that run on it every week.

*GTFS: General Transit Feed Specification
Acknowledgement to http://www.stevencanplan.com/2010/10/how-to-convert-gtfs-to-shapefiles-and-kml/ for inspiration
'''

import arcpy
import os
import zipfile
import logging
from datetime import datetime

arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

# rootdir = 'C:/Mathis/ICSL/stormwater/'
# gtfs_dir = os.path.join(rootdir, 'data\SoundTransit_201812\gtfs_puget_sound_consolidated.zip')
# out_gdb = os.path.join(rootdir, 'results/transit.gdb')
# out_fc = 'PStransit'
# current = True

def GTFStoSHPweeklynumber(gtfs_dir, out_gdb, out_fc, current=True, keep=False):
    #Create gdb
    if not arcpy.Exists(out_gdb):
        print('Creating {}...'.format(out_gdb))
        arcpy.CreateFileGDB_management(os.path.dirname(out_gdb), os.path.basename(out_gdb))
    else:
        print('{} already exists, using existing gdb...'.format(out_gdb))
    arcpy.env.workspace = out_gdb

    #Unzip GTFS
    if os.path.splitext(gtfs_dir)[1] == '.zip':
        print('Unzipping {}...'.format(os.path.split(gtfs_dir)[1]))
        zip = zipfile.ZipFile(gtfs_dir)
        zip.extractall(os.path.split(gtfs_dir)[0])
        gtfs_dir = os.path.split(gtfs_dir)[0]

    #Create log to write out errors (https://docs.python.org/3/howto/logging.html#logging-basic-tutorial)
    logging.basicConfig(filename=os.path.join(gtfs_dir, datetime.now().strftime('errorlog_%Y%m%d%H%M.log')),
                        level=logging.WARNING)

    #Import GTFS text files to out_gdb
    print('Importing GTFS tables into gdb...')
    arcpy.CopyRows_management(os.path.join(gtfs_dir,'routes.txt'), 'routes')
    print('routes.txt imported')
    arcpy.CopyRows_management(os.path.join(gtfs_dir,'trips.txt'), 'trips')
    print('trips.txt imported')
    arcpy.CopyRows_management(os.path.join(gtfs_dir,'calendar.txt'), 'calendar')
    print('calendar.txt imported')
    #arcpy.CopyRows_management(os.path.join(gtfs_dir,'stop_times.txt'), 'stop_times')

    #Create shapefile from shapes.txt
    print('Creating point feature class...')
    newfc = arcpy.CreateFeatureclass_management(out_gdb, out_fc, "Point")
    arcpy.DefineProjection_management(newfc , 4326)
    
    txt = open(os.path.join(gtfs_dir, 'shapes.txt'), 'r')
    inputF = txt.readlines()

    #Create columns based on those in shape.txt
    headers = inputF[0].replace('\n', '').split(',')
    arcpy.AddField_management(newfc, 'shape_id', 'TEXT', field_length = 50)
    arcpy.AddField_management(newfc, 'shape_pt_lat', 'FLOAT')
    arcpy.AddField_management(newfc, 'shape_pt_lon', 'FLOAT')
    arcpy.AddField_management(newfc, 'shape_pt_sequence', 'LONG')
    if 'shape_dist_traveled' in headers:
        arcpy.AddField_management(newfc, 'shape_dist_traveled', 'FLOAT')

    #Get column number of coordinates
    iy = headers.index('shape_pt_lat')
    ix = headers.index('shape_pt_lon')

    #Fill in point feature class
    with arcpy.da.InsertCursor(newfc, headers + ["SHAPE@XY"]) as cursor:
        for line in inputF[1:]:
            try:
                lineformat = line.replace('\n', '').split(',')
                lineformat.append((float(lineformat[ix]), float(lineformat[iy]))) #Add coordinate tuple
                cursor.insertRow(tuple(lineformat))
            except Exception as e:
                print(e)
                logging.warning('{0} \n Error row: {1}'.format(e, line))

    #Convert to polyline
    print('Converting points to lines...')
    newfc_lines = out_fc + '_lines'
    arcpy.PointsToLine_management(newfc, newfc_lines, Line_Field= 'shape_id')

    #Each route has multiple trips
    #Eacj trip is a unique bus line-time-week day combination (e.g. 6:47 am MON-FRI Line 71)
    #Compute number of days/week for each service_id in calendar
    print('Summaryzing data...')
    arcpy.AddField_management('calendar', 'daynum', 'SHORT')
    cal_fields = ['daynum', 'start_date', 'monday', 'tuesday', 'wednesday', 'thursday', 'friday', 'saturday', 'sunday']
    with arcpy.da.UpdateCursor('calendar', cal_fields) as cursor:
        for row in cursor:
            #Make sure that start_date of trip is < present
            if datetime.strptime(str(row[1]), '%Y%m%d') < datetime.now():
                row[0] = sum(row[2:])
            else:
                if current == True:
                    row[0] = 0
            cursor.updateRow(row)

    #Join calendar to trips
    arcpy.MakeTableView_management('trips', 'trips_view')
    arcpy.AddJoin_management('trips_view', 'service_id', 'calendar', 'service_id', join_type='KEEP_ALL')
    arcpy.AddJoin_management('trips_view', 'route_id', 'routes', 'route_id', join_type='KEEP_ALL')
    arcpy.CopyRows_management('trips_view', 'trips_routes_calendar')
    
    #Summary trips by route_id & shape_id (here add
    arcpy.Statistics_analysis('trips_routes_calendar', 'trips_routes_count', statistics_fields=[['daynum', 'SUM']],
                              case_field = ['shape_id', 'route_id', 'route_short_name', 'route_desc', 'route_type'])

    print('Joining and exporting data...')
    #Make sure that shape_id in join table is in text form
    for f in arcpy.ListFields('trips_routes_count'):
        if f.name == 'shape_id':
            if f.type != 'String':
                print('Converting shape_id field to string...')
                arcpy.AddField_management('trips_routes_count', 'shape_idstr', 'TEXT')
                arcpy.CalculateField_management('trips_routes_count', field = 'shape_idstr',
                                                expression = 'str(!shape_id!)', expression_type='PYTHON')
                arcpy.DeleteField_management('trips_routes_count', 'shape_id')
                arcpy.AlterField_management('trips_routes_count', 'shape_idstr', 'shape_id', 'shape_id')

    # Join trip statistics to lines
    arcpy.MakeFeatureLayer_management(newfc_lines, 'shapes_lyr')
    arcpy.AddJoin_management('shapes_lyr', 'shape_id', 'trips_routes_count', 'shape_id')
    arcpy.CopyFeatures_management('shapes_lyr', out_fc + '_routes')

    if keep == False:
        print('Deleting intermediate outputs...')
        for inter_lyr in ['routes', 'trips', 'calendar', newfc, newfc_lines]:
            arcpy.Delete_management(inter_lyr)