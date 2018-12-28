'''
Author: Mathis Messager
Contact info: messamat@uw.edu
Creation date: December 2018

Required Arguments:
    - kernel_dir (path): directory containing weighted kernel text files
    - in_raster (raster path): raster to run focal statistics on
    - out_gdb (path): file gdb where heatmap raster will be written
    - out_var (text): name of variable to name output file after (e.g. congestion, AADT)

Optional Arguments:
    - divnum (number, default = 1): number by which to divide heatmap values so as to keep numbers lower but at the
                        cost of precision (e.g. 10, 100, 1000)
    - keyw (text, default = ''): text pattern to subset kernels to run focal statistics on
                        (e.g. log, pow to run only on log kernels or power kernels)

Description: create integer heatmaps using arcpy focal statistics based on a set of weighted kernels contained in a
            directory
'''

import arcpy
import os
import re

def customheatmap(kernel_dir, in_raster, out_gdb, out_var, divnum=1, keyw=''):
    for kertxt in os.listdir(kernel_dir):
        if re.compile('kernel' + keyw).match(kertxt):
            outext = 'heat{0}{1}'.format(out_var, os.path.splitext(kertxt)[0][7:])
            if not arcpy.Exists(outext):
                print(outext)
                kernel = arcpy.sa.NbrWeight(os.path.join(kernel_dir, kertxt))
                heat_aadt = arcpy.sa.Int(
                    arcpy.sa.FocalStatistics(in_raster, neighborhood=kernel, statistics_type='SUM',
                                    ignore_nodata='DATA') / divnum + 0.5)
                heat_aadt.save(os.path.join(out_gdb, outext))