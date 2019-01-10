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
    - keep (True/False): whether to keep (True) intermediate outputs or not (False). Note that the error log is deleted
                            only if it is empty.

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
from collections import defaultdict
from datetime import datetime

arcpy.env.overwriteOutput=True
arcpy.env.qualifiedFieldNames = False

# rootdir = 'C:/Mathis/ICSL/stormwater/'
# gtfs_dir = os.path.join(rootdir, 'data\\TransitWiki201812\\mason-wa-us.zip') #clallam_google_transit.zip
# out_gdb = os.path.join(rootdir, 'results/transit.gdb')
# out_fc = 'jefferson_google_transit'
# current = True
# keep = False

def delete_intoutput(dir, ziplist, log, intlist) :
    print('Deleting intermediate outputs...')
    if os.stat(log).st_size == 0L: #Delete log if empty
        os.remove(log)
    for file in ziplist:  # Delete zipped out files
        os.remove(os.path.join(dir, file))
    for inter_lyr in intlist:  # Delete intermediate GIS layers
        try:
            arcpy.Delete_management(inter_lyr)
        except Exception as e:
            print(e)
            pass

def GTFStoSHPweeklynumber(gtfs_dir, out_gdb, out_fc, current=True, keep=False):
    #Create gdb
    if not arcpy.Exists(out_gdb):
        print('Creating {}...'.format(out_gdb))
        arcpy.CreateFileGDB_management(os.path.dirname(out_gdb), os.path.basename(out_gdb))
    else:
        print('{} already exists, using existing gdb...'.format(out_gdb))
    arcpy.env.workspace = out_gdb
    gtfs_rootdir = os.path.split(gtfs_dir)[0]

    if os.path.exists(gtfs_dir) and \
            (os.path.isdir(gtfs_dir) or zipfile.is_zipfile(gtfs_dir)): #Make sure that path exists and is a dir
        try: #Try process, if fails, delete intermediate output if keep=FALSE
            #Unzip GTFS
            if zipfile.is_zipfile(gtfs_dir):
                print('Unzipping {}...'.format(os.path.split(gtfs_dir)[1]))
                with zipfile.ZipFile(gtfs_dir) as zipf:
                    zipfilelist = [info.filename for info in zipf.infolist()]
                    print('Overwriting {}...'.format(', '.join(
                        [f for f in zipfilelist if os.path.exists(os.path.join(gtfs_rootdir, f))])))
                    zipf.extractall(os.path.split(gtfs_dir)[0])
                del zipf

            #Create log to write out errors (https://docs.python.org/3/howto/logging.html#logging-basic-tutorial)
            errorlog = os.path.join(gtfs_rootdir, datetime.now().strftime('errorlog_%Y%m%d%H%M.log'))
            logging.basicConfig(filename= errorlog, level=logging.WARNING)

            #Import GTFS text files to out_gdb
            print('Importing GTFS tables into gdb...')
            arcpy.CopyRows_management(os.path.join(gtfs_rootdir,'routes.txt'), 'routes')
            print('routes.txt imported')
            arcpy.CopyRows_management(os.path.join(gtfs_rootdir,'trips.txt'), 'trips')
            print('trips.txt imported')
            arcpy.CopyRows_management(os.path.join(gtfs_rootdir,'calendar.txt'), 'calendar')
            print('calendar.txt imported')
            # arcpy.CopyRows_management(os.path.join(gtfs_rootdir, 'calendar_dates.txt'), 'calendar') #Imports in wrong format
            # print('calendar_dates.txt imported')
            #arcpy.CopyRows_management(os.path.join(gtfs_rootdir,'stop_times.txt'), 'stop_times')

            #-------------------------------------------------------------------------------------------------------------------
            # CREATE SHAPEFILE FROM SHAPES.TXT
            #-------------------------------------------------------------------------------------------------------------------
            print('Creating point feature class...')

            newfc = arcpy.CreateFeatureclass_management(out_gdb, out_fc, "Point")
            arcpy.DefineProjection_management(newfc , 4326)

            with open(os.path.join(gtfs_rootdir, 'shapes.txt'), 'r') as shp_txt:
                shp_F = shp_txt.readlines()

            #Create columns based on those in shape.txt
            headers = shp_F[0].replace('\n', '').split(',')
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
                for line in shp_F[1:]:
                    try:
                        lineformat = line.replace('\n', '').split(',')
                        lineformat.append((float(lineformat[ix]), float(lineformat[iy]))) #Add coordinate tuple
                        cursor.insertRow(tuple(lineformat))
                    except Exception as e:
                        print(e)
                        logging.warning('{0} \n Error row: {1}'.format(e, line))
            del shp_txt
            del shp_F

            #Convert to polyline
            print('Converting points to lines...')
            newfc_lines = out_fc + '_lines'
            arcpy.PointsToLine_management(newfc, newfc_lines, Line_Field= 'shape_id')

            #-------------------------------------------------------------------------------------------------------------------
            # COMPUTE NUMBER OF WEEKLY TRIPS PER SHAPE
            #-------------------------------------------------------------------------------------------------------------------
            print('Summaryzing data...')

            #Analyze calendar dates for trips with no calendar schedule
            '''
            In some cases, calendar.txt was entirely omitted for a trip, and ALL dates of service were included in 
            calendar_dates.txt (https://developers.google.com/transit/gtfs/reference/#calendar_datestxt).
            This happens when schedule varies most days of the month, or the agency wants to programmatically output service 
            dates without specifying a normal weekly schedule. In most cases however, calendar_dates.txt is used to define 
            exceptions to the default service categories defined in calendar.txt file. 
            '''
            with open(os.path.join(gtfs_rootdir, 'calendar_dates.txt'), 'r') as cal_txt:
                cal_F = cal_txt.readlines()
                cal_headers = cal_F[0].replace('\n', '').split(',')

            # Create a dictionary to store number of dates for each exception_type: 1 service added, 2 service removed
            cal_dic = defaultdict(lambda: [0]*2)
            for line in cal_F[1:]:
                lineformat = dict(zip(cal_headers, line.replace('\n', '').split(',')))
                cal_dic[lineformat['service_id']][int(lineformat['exception_type'])-1] += 1

            '''Each route has multiple trips
               Each trip is a unique bus line-time-week day combination (e.g. 6:47 am MON-FRI Line 71)
            '''

            #Compute number of days/week for each service_id in calendar
            arcpy.AddField_management('calendar', 'normalnum', 'SHORT')
            arcpy.AddField_management('calendar', 'date_added', 'SHORT')
            arcpy.AddField_management('calendar', 'date_removed', 'SHORT')
            arcpy.AddField_management('calendar', 'date_weekavg', 'FLOAT')
            arcpy.AddField_management('calendar', 'adjustnum', 'FLOAT')
            arcpy.AddField_management('calendar', 'service_len', 'SHORT')

            cal_fields = ['service_id',  'start_date','end_date', 'normalnum', #0, 1, 2, 3
                          'monday', 'tuesday', 'wednesday', 'thursday', 'friday', 'saturday', 'sunday', #4-10
                          'date_added', 'date_removed', 'date_weekavg', 'adjustnum', 'service_len'] #11, 12, 13, 14, 15

            with arcpy.da.UpdateCursor('calendar', cal_fields) as cursor:
                for row in cursor:
                    service_len = max((datetime.strptime(str(row[2]), '%Y%m%d') -
                                       datetime.strptime(str(row[1]), '%Y%m%d')).days,
                                      1)  # Compute length of service in days. If service runs for only one day, assign 1
                    row[15] = service_len

                    #print(row[1])
                    if row[0] in cal_dic: #if service_id in calendar_dates dictionary
                        row[11] = cal_dic[row[0]][0] #date_added = number of added service dates in calendar_dates.txt (type 1)
                        row[12] = cal_dic[row[0]][1] #date_removed = number of removed service dates in calendar_dates.txt (type 2)
                        row[13] = (cal_dic[row[0]][0]-cal_dic[row[0]][1])/(float(service_len)/7.0) #Compute average net number of exceptions/week
                    else:
                        row[11:14] = (0,0,0)
                    cursor.updateRow(row) #update row halfway to avoid having to retype row[13] equation

                    row[3] = sum(row[4:11]) #normalnum = sum(monday-sunday trips)
                    row[14] = max(sum(row[4:11]) + row[13], 0) #adjustnum = normalnum + date_weekavg

                    if datetime.strptime(str(row[1]), '%Y%m%d') > datetime.now() and current == True:  # If start_date > current date
                        row[3] = 0
                        row[14] = 0
                    cursor.updateRow(row)

            #Join calendar to trips
            arcpy.MakeTableView_management('trips', 'trips_view')
            arcpy.AddJoin_management('trips_view', 'service_id', 'calendar', 'service_id', join_type='KEEP_ALL')
            arcpy.AddJoin_management('trips_view', 'route_id', 'routes', 'route_id', join_type='KEEP_ALL')
            arcpy.CopyRows_management('trips_view', 'trips_routes_calendar')

            #Summarize trips by route_id & shape_id (here add
            arcpy.Statistics_analysis('trips_routes_calendar', 'trips_routes_count',
                                      statistics_fields=[['normalnum', 'SUM'], ['adjustnum', 'SUM'],
                                                         #['direction_id', 'FIRST'], ['direction_id', 'LAST'],
                                                         ['service_len', 'MIN'], ['service_len', 'MAX']],
                                      case_field = ['shape_id', 'route_id', 'route_short_name',
                                                    'route_desc', 'route_type'])

            #-------------------------------------------------------------------------------------------------------------------
            # JOIN AND EXPORT DATE
            #-------------------------------------------------------------------------------------------------------------------
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
                delete_intoutput(dir=gtfs_rootdir, ziplist=zipfilelist, log = errorlog,
                                 intlist=['routes', 'trips', 'trips_view', 'calendar', 'shapes_lyr',
                                            'trips_routes_calendar', 'trips_routes_count', out_fc, newfc, newfc_lines])

        #If an error is raised, delete intermediate outputs if keep == False
        except Exception as e:
            print(e)
            if keep == False:
                delete_intoutput(dir = gtfs_rootdir, ziplist = zipfilelist, log = errorlog,
                                 intlist = ['routes', 'trips', 'trips_view', 'calendar', 'shapes_lyr',
                                            'trips_routes_calendar', 'trips_routes_count', out_fc, newfc, newfc_lines])
    else:
        raise Exception("{} does not exist or is not a directory or zip file".format(gtfs_dir))