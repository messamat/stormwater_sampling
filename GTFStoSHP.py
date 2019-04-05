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
import pandas as pd
import fiona
import os
import re
import csv
import zipfile
import logging
import dateutil.parser
import traceback
from collections import defaultdict
from datetime import datetime

#--------------------------------------------------------------
# FOR TROUBLESHOOTING
# rootdir = 'D:/Mathis/ICSL/stormwater/'
# gtfs_dir= os.path.join(rootdir, 'data/NTM_0319')
# out_gdb=os.path.join(rootdir, 'results/NTM.gdb')
# out_fc = 'NTM'
# keep = True
# outdir = gtfs_dir
# infile = os.path.join(gtfs_dir,'NTMAPI_trips.csv')
# arcpy.env.workspace = out_gdb
# shapesf = 'D:\Mathis\ICSL\stormwater\data\NTM0319/NTM_shapes.shp'
#--------------------------------------------------------------

def format_schema_import(infile, outdir):
    filebase = os.path.split(infile)[1]

    #Read infile headers
    with open(infile, 'r') as txt:
        dialect = csv.Sniffer().sniff(txt.read())
        txt.seek(0)
        txtF = csv.DictReader(txt, delimiter=dialect.delimiter)
        txtH = txtF.fieldnames

        # Write schema.ini file
        idfields = ['route_id', 'service_id', 'trip_id', 'shape_id', 'block_id']
        for i in range(0, len(txtH)):
            fname = txtH[i]
            #If ID field in text file
            if fname in idfields:
                print(fname)
                maxflen = max([len(row[fname]) for row in txtF]) #maximum length of id field
                schema_header = '[{0}]\nFormat=CSVDelimited\nColNameHeader=True'.format(filebase)
                schema_addline = '\nCol{0}={1} Text Width {2}'.format(i + 1, fname, maxflen)

                #If a schema file already exists
                if 'schema.ini' in os.listdir(outdir):
                    with open(os.path.join(outdir, 'schema.ini'), 'r') as schematxt:
                        schemaF = schematxt.readlines()
                        allines = ('').join(schemaF)

                    postfile_regex = '[[]{}[]](?s).*(?=[[]|$)'.format(filebase)
                    col_regex = 'Col{0}.*(?=\n|$)'.format(i + 1)

                    # If txt file already in schema ini
                    if re.compile(postfile_regex).search(allines):
                        ftxt = re.compile(postfile_regex).search(allines).group()

                        with open(os.path.join(outdir, 'schema.ini'), 'w') as schematxt:
                            # If column in schema ini
                            if re.compile(col_regex).search(ftxt):
                                part_rep = re.sub(col_regex, 'Col{0}={1} Text Width {2}'.format(i + 1, fname, maxflen), ftxt)
                                full_rep = re.sub(postfile_regex, part_rep, allines)

                            else: #If column not in schema ini
                                full_rep = re.sub(postfile_regex, '{0}{1}'.format(ftxt, schema_addline), allines)
                            schematxt.write(full_rep)

                    else: #If txt file not already in schema ini
                        with open(os.path.join(outdir, 'schema.ini'), 'a') as schematxt:
                            schematxt.write('\n{0}{1}'.format(schema_header, schema_addline))

                #If a schema file doesn't already exist
                else:
                    with open(os.path.join(outdir, "schema.ini"), "w") as schematxt:
                        schematxt.write('{0}{1}'.format(schema_header, schema_addline))

                #Rewind table reader
                txt.seek(0)

def delete_intoutput(dir, ziplist, intlist) :
    print('Deleting intermediate outputs...')
    os.remove(os.path.join(dir, 'schema.ini'))
    for file in ziplist:  # Delete zipped out files
        os.remove(os.path.join(dir, file))
    for inter_lyr in intlist:  # Delete intermediate GIS layers
        try:
            arcpy.Delete_management(inter_lyr)
        except Exception:
            traceback.print_exc()
            pass

def GTFStoSHPweeklynumber(gtfs_dir, out_gdb, out_fc, current=True, keep=False):
    arcpy.env.overwriteOutput = True
    arcpy.env.qualifiedFieldNames = False

    #Create gdb
    if not arcpy.Exists(out_gdb):
        print('Creating {}...'.format(out_gdb))
        arcpy.CreateFileGDB_management(os.path.dirname(out_gdb), os.path.basename(out_gdb))
    else:
        print('{} already exists, using existing gdb...'.format(out_gdb))
    arcpy.env.workspace = out_gdb

    if os.path.exists(gtfs_dir) and \
            (os.path.isdir(gtfs_dir) or zipfile.is_zipfile(gtfs_dir)): #Make sure that path exists and is a dir
        try: #Try process, if fails, delete intermediate output if keep=FALSE
            #Unzip GTFS
            if zipfile.is_zipfile(gtfs_dir):
                print('Unzipping {}...'.format(os.path.split(gtfs_dir)[1]))
                with zipfile.ZipFile(gtfs_dir) as zipf:
                    zipfilelist = [info.filename for info in zipf.infolist()]
                    print('Overwriting {}...'.format(', '.join(
                        [f for f in zipfilelist if os.path.exists(os.path.join(gtfs_dir, f))])))
                    zipf.extractall(os.path.split(gtfs_dir)[0])
                del zipf
                gtfs_rootdir = os.path.split(gtfs_dir)[0]
            else:
                gtfs_rootdir = gtfs_dir

            #Import GTFS text files to out_gdb
            print('Importing GTFS tables into gdb...')
            gtfs_files = [os.path.join(dirpath, file)
                          for (dirpath, dirnames, filenames) in os.walk(gtfs_rootdir)
                          for file in filenames]

            routesf = filter(re.compile('.*routes[.](csv|txt)$').search, gtfs_files)[0]
            format_schema_import(infile=routesf, outdir=gtfs_rootdir)
            arcpy.TableToTable_conversion(routesf, out_gdb, 'routes')
            print('routes imported')
            tripsf = filter(re.compile('.*trips[.](csv|txt)$').search, gtfs_files)[0]
            format_schema_import(infile=tripsf, outdir=gtfs_rootdir)
            arcpy.TableToTable_conversion(tripsf, out_gdb, 'trips')
            print('trips imported')
            calendarf = filter(re.compile('.*calendar[.](csv|txt)$').search, gtfs_files)[0]
            format_schema_import(infile=calendarf, outdir=gtfs_rootdir)
            arcpy.TableToTable_conversion(calendarf, out_gdb, 'calendar')
            print('calendar imported')
            # arcpy.CopyRows_management(os.path.join(gtfs_rootdir, 'calendar_dates.txt'), 'calendar') #Imports in wrong format
            # print('calendar_dates.txt imported')
            #arcpy.CopyRows_management(os.path.join(gtfs_rootdir,'stop_times.txt'), 'stop_times')

            #-------------------------------------------------------------------------------------------------------------------
            # CREATE SHAPEFILE FROM SHAPES.TXT IF NEEDED
            #-------------------------------------------------------------------------------------------------------------------
            print('Creating point feature class...')

            shapesf = filter(re.compile('.*shapes[.](csv|txt|shp)$').search, gtfs_files)[0]

            if os.path.splitext(shapesf)[1] in ['.csv', '.txt']:
                newfc = arcpy.CreateFeatureclass_management(out_gdb, out_fc, "Point")
                arcpy.DefineProjection_management(newfc, 4326)

                with open(shapesf, 'r') as shp_txt:
                    dialect = csv.Sniffer().sniff(shp_txt.read())
                    shp_txt.seek(0)
                    txtF = csv.DictReader(shp_txt, delimiter=dialect.delimiter)
                    headers = txtF.fieldnames

                    #Create columns based on those in shape.txt
                    arcpy.AddField_management(newfc, 'shape_id', 'TEXT', field_length = 50)
                    arcpy.AddField_management(newfc, 'shape_pt_lat', 'FLOAT')
                    arcpy.AddField_management(newfc, 'shape_pt_lon', 'FLOAT')
                    arcpy.AddField_management(newfc, 'shape_pt_sequence', 'LONG')
                    if 'shape_dist_traveled' in headers:
                        arcpy.AddField_management(newfc, 'shape_dist_traveled', 'FLOAT')

                    #Fill in point feature class
                    with arcpy.da.InsertCursor(newfc, headers + ["SHAPE@XY"]) as cursor:
                        for line in txtF:
                            try:
                                lineformat = [line['shape_id'],
                                              line['shape_pt_lat'],
                                              line['shape_pt_lon'],
                                              line['shape_pt_sequence']]
                                if 'shape_dist_traveled' in headers:
                                    lineformat.append(line['shape_dist_traveled'])
                                lineformat.append((float(line['shape_pt_lon']), float(line['shape_pt_lat']))) #Add coordinate tuple
                                cursor.insertRow(tuple(lineformat))
                            except Exception as e:
                                traceback.print_exc()
                                logging.warning('{0} \n Error row: {1}'.format(e, line))

                #Convert to polyline
                print('Converting points to lines...')
                newfc_lines = out_fc + '_lines'
                arcpy.PointsToLine_management(newfc, newfc_lines, Line_Field= 'shape_id')
            else: #If already a shapefile simply create a copy
                newfc_lines = out_fc + '_lines'
                arcpy.CopyFeatures_management(shapesf, newfc_lines)

            #-------------------------------------------------------------------------------------------------------------------
            # COMPUTE NUMBER OF WEEKLY TRIPS PER SHAPE
            #-------------------------------------------------------------------------------------------------------------------
            print('Summaryzing data...')

            #Analyze calendar dates for trips with no calendar schedule
            '''
            In some cases, calendar was entirely omitted for a trip, and ALL dates of service were included in 
            calendar_dates (https://developers.google.com/transit/gtfs/reference/#calendar_datestxt).
            This happens when schedule varies most days of the month, or the agency wants to programmatically output service 
            dates without specifying a normal weekly schedule. In most cases however, calendar_dates is used to define 
            exceptions to the default service categories defined in calendar file. 
            '''
            calendar_datesf = filter(re.compile('.*calendar_dates[.](csv|txt)').search, gtfs_files)[0]
            with open(calendar_datesf, 'r') as cal_txt:
                dialect = csv.Sniffer().sniff(cal_txt.read())
                cal_txt.seek(0)
                cal_F = csv.DictReader(cal_txt, delimiter=dialect.delimiter)
                #cal_headers = cal_F.fieldnames

                # Create a dictionary to store number of dates for each exception_type (1 service added, 2 service removed) and dates
                cal_dic = defaultdict(lambda: [0, 0, datetime.max, datetime.min])
                #cal_dic record structure: {AgencyName._.service_id: # of added dates, # of removed dates, earliest date, latest date}
                for line in cal_F:
                    dkey = '{0}._.{1}'.format(line['AgencyName'], line['service_id'])
                    cal_dic[dkey][int(line['exception_type']) - 1] += 1
                    #print(dkey)

                    try:
                        service_date = datetime.strptime(line['date'], '%Y%m%d')
                        if service_date < cal_dic[dkey][2]:
                            cal_dic[dkey][2] = service_date
                        if service_date > cal_dic[dkey][3]:
                            cal_dic[dkey][3] = service_date
                    except:
                        pass
                        traceback.print_exc()
                        print('Continuing...')

            '''Each route has multiple trips
               Each trip is a unique bus line-time-week day combination (e.g. 6:47 am MON-FRI Line 71)
            '''

            #Compute number of days/week for each service_id in calendar
            arcpy.AddField_management('calendar', 'normalnum', 'SHORT')
            arcpy.AddField_management('calendar', 'date_added', 'SHORT')
            arcpy.AddField_management('calendar', 'date_removed', 'SHORT')
            arcpy.AddField_management('calendar', 'date_weekavg', 'FLOAT')
            arcpy.AddField_management('calendar', 'adjustnum', 'FLOAT')
            arcpy.AddField_management('calendar', 'service_len', 'LONG')

            cal_fields = ['AgencyName', 'service_id',  'start_date','end_date', 'normalnum', #0, 1, 2, 3
                          'monday', 'tuesday', 'wednesday', 'thursday', 'friday', 'saturday', 'sunday', #4-10
                          'date_added', 'date_removed', 'date_weekavg', 'adjustnum', 'service_len'] #11, 12, 13, 14, 15
            #fl = [f.name for f in arcpy.ListFields('calendar')]
            #[cf for cf in cal_fields if cf not in fl]

            #Correct bug in table whereby the end_date for a given agency includes a line break \n and the service_id
            with arcpy.da.UpdateCursor('calendar', ['AgencyName', 'end_date']) as cursor:
                for row in cursor:
                    if row[0] == 'Nashua_10087_132_914':
                        row[1] = row[1].split(',')[-1]
                        cursor.updateRow(row)

            cal_ids = [] #list to store list of service_ids in calendar table
            with arcpy.da.UpdateCursor('calendar', cal_fields) as cursor:
                for row in cursor:
                    #Create unique key AgencyName + service_id
                    dkey = '{0}._.{1}'.format(row[0], row[1])
                    #print(dkey)

                    if row[2] != None and row[3] != None: #Make sure that start and end date exist (i.e. not added for calendar_dates)
                        cal_ids.append(dkey)
                        # Compute length of service in days. If service runs for only one day, assign 1
                        try:
                            service_len = max((datetime.strptime(str(row[3]), '%Y%m%d') -
                                               datetime.strptime(str(row[2]), '%Y%m%d')).days,
                                              1)
                        except: #If fails, try automatically parsing the date
                            service_len = max((dateutil.parser.parse(row[3]) -
                                               dateutil.parser.parse(row[2])).days,
                                              1)
                        row[16] = service_len

                        #print(row[1])
                        if dkey in cal_dic: #if service_id in calendar_dates dictionary
                            row[12] = cal_dic[dkey][0] #date_added = number of added service dates in calendar_dates.txt (type 1)
                            row[13] = cal_dic[dkey][1] #date_removed = number of removed service dates in calendar_dates.txt (type 2)
                            row[14] = (cal_dic[dkey][0]-cal_dic[dkey][1])/(float(service_len)/7.0) #Compute average net number of exceptions/week
                        else:
                            row[12:15] = (0,0,0)
                        cursor.updateRow(row) #update row halfway to avoid having to retype row[14] equation

                        row[4] = sum(row[5:12]) #normalnum = sum(monday-sunday trips)
                        row[15] = max(sum(row[5:12]) + row[14], 0) #adjustnum = normalnum + date_weekavg

                        try:
                            sdate = datetime.strptime(str(row[2]), '%Y%m%d')
                        except:
                            sdate = dateutil.parser.parse(row[2])

                        if sdate > datetime.now() and current == True:  # If start_date > current date
                            row[4] = 0
                            row[15] = 0
                        cursor.updateRow(row)


            #For many service_id, calendar_dates contains all records while calendar does not contain the service_id
            #Making it impossible to join it to the calendar table. Therefore, we append these records to the calendar
            #table.
            #Compute average weekly number from calendar_dates by service_id but only keep those that span at least a month

            with arcpy.da.InsertCursor('calendar', ['AgencyName', 'service_id', 'adjustnum', 'service_len']) as cursor:
                for dkey in cal_dic:
                    if dkey not in cal_ids:
                        service_len = max((cal_dic[dkey][3] - cal_dic[dkey][2]).days, 1)
                        if service_len > 30:
                            newrow = (dkey.split('._.')[0],
                                      dkey.split('._.')[1],
                                      max((cal_dic[dkey][0]-cal_dic[dkey][1])/(float(service_len)/7.0),0),
                                      service_len)
                            cursor.insertRow(newrow)

            #Join calendar and routes to trips
            arcpy.AddField_management('trips', 'Agency_service_id', 'TEXT')
            arcpy.AddField_management('trips', 'Agency_route_id', 'TEXT')
            with arcpy.da.UpdateCursor('trips', ['Agency_service_id', 'AgencyName', 'service_id',
                                                 'Agency_route_id', 'route_id']) as cursor:
                for row in cursor:
                    row[0] = '{0}._.{1}'.format(row[1], row[2])
                    row[3] = '{0}._.{1}'.format(row[1], row[4])
                    cursor.updateRow(row)

            arcpy.AddField_management('calendar', 'Agency_service_id', 'TEXT')
            with arcpy.da.UpdateCursor('calendar', ['Agency_service_id', 'AgencyName', 'service_id']) as cursor:
                for row in cursor:
                    row[0] = '{0}._.{1}'.format(row[1], row[2])
                    cursor.updateRow(row)

            arcpy.AddField_management('routes', 'Agency_route_id', 'TEXT')
            with arcpy.da.UpdateCursor('routes', ['Agency_route_id', 'AgencyName', 'route_id']) as cursor:
                for row in cursor:
                    row[0] = '{0}._.{1}'.format(row[1], row[2])
                    cursor.updateRow(row)

            arcpy.MakeTableView_management('trips', 'trips_view')
            arcpy.AddJoin_management('trips_view', 'Agency_service_id', 'calendar', 'Agency_service_id', join_type='KEEP_ALL')
            arcpy.AddJoin_management('trips_view', 'Agency_route_id', 'routes', 'Agency_route_id', join_type='KEEP_ALL')
            arcpy.CopyRows_management('trips_view', 'trips_routes_calendar')

            ######### Summarize trips by route_id & shape_id #############
            #Compute statistics in dictionary
            trc_dic = defaultdict(lambda: [0, 0, 1000000, 0])
            cfields = ['shape_id', 'AgencyName', 'route_id', 'route_type', 'normalnum', 'adjustnum', 'service_len']
            fcount = arcpy.GetCount_management('trips_routes_calendar').getOutput(0)
            with arcpy.da.SearchCursor('trips_routes_calendar', cfields) as cursor:
                x=0
                for row in cursor:
                    ukey = '{0}._.{1}._.{2}._.{3}'.format(str(row[0]), str(row[1]), str(row[2]), str(row[3]))
                    # 'normalnum', 'SUM'
                    try:
                        trc_dic[ukey][0] += row[4]
                    except Exception:
                        traceback.print_exc()
                    # 'adjustnum', 'SUM'
                    try:
                        trc_dic[ukey][1] += row[5]
                    except Exception:
                        traceback.print_exc()
                    # service_len MIN
                    trc_dic[ukey][2] = row[6] if row[6] < trc_dic[ukey][2] else trc_dic[ukey][2]
                    # service_len MAX
                    trc_dic[ukey][3] = row[6] if row[6] > trc_dic[ukey][3] else trc_dic[ukey][3]

                    if x % 100 == 0:
                        print('Processing... {}%'.format(round(100*x/float(fcount), 2)))
                    x += 1

            #Convert dictionary to geodatabase table
            arcpy.CreateTable_management(out_gdb, 'trips_routes_count')
            tempfields = {}
            for f in arcpy.ListFields('trips_routes_calendar'):
                tempfields[f.name] = f.type
            outfields = ['shape_id', 'AgencyName', 'route_id', 'route_type'] + \
                     ['normalnum_SUM', 'adjustnum_SUM', 'service_len_MIN', 'service_len_MAX']
            for f in outfields:
                if f in tempfields:
                    arcpy.AddField_management('trips_routes_count', f, field_type = tempfields[f])
                else:
                    arcpy.AddField_management('trips_routes_count', f, field_type = 'FLOAT')

            with arcpy.da.InsertCursor('trips_routes_count', outfields) as cursor:
                for k in trc_dic:
                    try:
                        ksplit = k.split('._.')
                        lineformat = ksplit + trc_dic[k]
                        cursor.insertRow(tuple(lineformat))
                    except Exception as e:
                        traceback.print_exc()
                        logging.warning('{0} \n Error row: {1}'.format(e, line))
            #-------------------------------------------------------------------------------------------------------------------
            # JOIN AND EXPORT DATA
            #-------------------------------------------------------------------------------------------------------------------
            print('Joining and exporting data...')

            # Join trip statistics to lines by Agency Name and shape_id
            arcpy.AddField_management(newfc_lines, 'Agency_shape_id', 'TEXT')
            with arcpy.da.UpdateCursor(newfc_lines, ['Agency_shape_id', 'AgencyName', 'shape_id']) as cursor:
                for row in cursor:
                    row[0] = '{0}._.{1}'.format(row[1], row[2])
                    cursor.updateRow(row)

            arcpy.AddField_management('trips_routes_count', 'Agency_shape_id', 'TEXT')
            with arcpy.da.UpdateCursor('trips_routes_count', ['Agency_shape_id', 'AgencyName', 'shape_id']) as cursor:
                for row in cursor:
                    row[0] = '{0}._.{1}'.format(row[1], row[2])
                    cursor.updateRow(row)

            arcpy.MakeFeatureLayer_management(newfc_lines, 'shapes_lyr')
            arcpy.AddJoin_management('shapes_lyr', 'Agency_shape_id', 'trips_routes_count', 'Agency_shape_id')
            arcpy.CopyFeatures_management('shapes_lyr', out_fc + '_routes')

            # -------------------------------------------------------------------------------------------------------------------
            # DELETE INTERMEDIATE STUFF
            # -------------------------------------------------------------------------------------------------------------------
            for inter_lyr in ['trips_view', 'shapes_lyr']:  # Delete intermediate GIS layers
                try:
                    arcpy.Delete_management(inter_lyr)
                except Exception as e:
                    traceback.print_exc()
                    logging.error(e)
                    pass

            if keep == False:
                delete_intoutput(dir=gtfs_rootdir, ziplist=zipfilelist,
                                 intlist=['routes', 'trips', 'calendar', 'trips_routes_calendar',
                                          'trips_routes_count', out_fc, newfc, newfc_lines])

        #If an error is raised, delete intermediate outputs if keep == False
        except Exception as e:
            traceback.print_exc()
            logging.error(e)
            if keep == False:
                delete_intoutput(dir = gtfs_rootdir, ziplist = zipfilelist,
                                 intlist = ['routes', 'trips', 'trips_view', 'calendar', 'shapes_lyr',
                                            'trips_routes_calendar', 'trips_routes_count', out_fc, newfc, newfc_lines])
    else:
        raise Exception("{} does not exist or is not a directory or zip file".format(gtfs_dir))

# GTFStoSHPweeklynumber(gtfs_dir, out_gdb, out_fc, current=True, keep=False)