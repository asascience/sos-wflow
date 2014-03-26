#!/usr/bin/env python

import shutil
import requests
import time as clock
from datetime import datetime
from datetime import date
from datetime import timedelta
import json
import os
import bisect
import sys
from xml.etree import ElementTree as ET
import argparse
import fileinput


import netCDF4
import numpy as np

import logging
logger = logging.getLogger("pytools")
logger.addHandler(logging.NullHandler())




def Create(start_date, minlon, minlat, maxlon, maxlat, root_directory, ago):
    #-------------------------------------------------------------------------------------------------------

    #Sets start and end date.

    if start_date == 'DEFAULT':
    	start_date = date.today()
    date_of_concern = start_date-timedelta(int(ago))
    end_date = str(date_of_concern)+" "+"23"+":"+"59"+":"+"59"+"Z"
    start_date = str(date_of_concern)+" "+"00"+":"+"00"+":"+"01"+"Z"


    #-------------------------------------------------------------------------------------------------------
    #Scrapes the data
    print "Gathering SLDMB data from %s to %s" %(start_date,end_date)
    print "Bounded by the following coordinates:"
    print "Minimum Longitude: %s" %minlon
    print "Maximum Longitude: %s" %maxlon
    print "Minimum Loatitude: %s" %minlat
    print "Maximum Latitude: %s" %maxlat
    payload = {
               "MinLon" : minlon,
               "MinLat"  : minlat, 
               "MaxLon"  : maxlon, 
               "MaxLat"  : maxlat, 
               "Start"   : start_date, 
               "End"     : end_date
               }
    url = "http://staging.coastmap.com/sldmb" 
    r=requests.get(url, params=payload)
    r.raise_for_status()
    json_data=json.loads(r.text)
    json_data=json_data[u'usp_buoydataboxandhoursResults'][u'usp_buoydataboxandhoursResult']
    #-------------------------------------------------------------------------------------------------------
    #Finds every unique sensor and if the height was recorded improperly, set it to 0
    json_buoy = []
    for each_json in json_data:
        json_buoy.append(each_json['BuoyNum'])
        if np.float32(each_json["AltSensor"]) > 10:
            each_json["AltSensor"] = 0
        
    unique_buoys = {}
    for buoy in json_buoy:
        if buoy not in unique_buoys:
            unique_buoys[buoy] = []
    for each_json_2 in json_data:
        for each_buoy in unique_buoys:
            if each_buoy == each_json_2['BuoyNum']:
                unique_buoys[each_buoy].append(each_json_2)
                
                
    for buoy in unique_buoys:
        output_directory = root_directory+"/%s" %buoy
        output_filename = "%s_%s.nc" %(buoy, date_of_concern)
        full_station_urn = "urn:ioos:station:maracoos:%s" %buoy
        full_sensor_urn = "urn:ioos:sensor:maracoos:%s" %buoy
        global_attributes = []
        attributes = []
        _FillValue = -999.9
        
        buoy_time = []
        buoy_sst = []
        buoy_verts = []
        buoy_lon = []
        buoy_lat = []
        buoy_cycle = []
        buoy_id = []
        buoy_fom = []
        buoy_quality = []
        
        for each_unique_buoy in (unique_buoys[buoy]):
            buoy_time.append(each_unique_buoy['DataDate'])
            buoy_verts.append(each_unique_buoy['AltSensor'])
            buoy_sst.append(float(each_unique_buoy['SST'])+273.15)
            buoy_lat.append(each_unique_buoy['Lat'])
            buoy_lon.append(each_unique_buoy['Lon'])
            buoy_cycle.append(each_unique_buoy['Cycle'])
            buoy_id.append(each_unique_buoy['DataID'])
            buoy_fom.append(each_unique_buoy['FOM'])
            buoy_quality.append(each_unique_buoy['Quality'])
        for p in range(len(buoy_time)):
            buoy_time[p] = clock.mktime(clock.strptime(buoy_time[p], "%Y-%m-%d %H:%M:%SZ"))
        longitude = np.mean(np.float32(buoy_lon))
        latitude = np.mean(np.float32(buoy_lat))
        precise_longitude = np.asarray(buoy_lon).astype(np.float)
        precise_latitude = np.asarray(buoy_lat).astype(np.float)
        data = (buoy_time, buoy_verts, buoy_sst)
        
        
        if data is not None:
            try:
                assert len(data) > 0
                # Data was passed in as a tuple of lists (time, vertical, value). Conver to numpy arrays.
                times     = np.asarray(data[0]).astype(np.float)
                verticals = np.ma.asarray(data[1]).astype(np.float)
                values    = np.ma.asarray(data[2]).astype(np.float)
        
            except (AssertionError, IndexError):
                logger.warn("No data passed in for '%s'.  Skipping file creation" % full_sensor_urn)
        
        
        assert times.size == verticals.size == values.size
        
        # Get unique time and verticals (the data used for the each variable)
    
        
        
        starting = (times[0])
        ending   = (times[-1])
        
        # Make directory
        if not os.path.exists(output_directory):
            os.makedirs(output_directory) 
            new_buoy_agg = '    <dataset name="Buoy %s Sea Surface Temperatures" ID="buoy%s-agg" urlPath="buoy%s-agg.ncml" dataType="Point">\n      <metadata inherited="true">\n        <serviceName>all</serviceName>\n        <dataType>Point</dataType>\n      </metadata>\n      <netcdf xmlns="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2">\n        <attribute name="title" value="SLDMB Buoy %s Aggregated Sea Surface Temperatures" />\n        <attribute name="summary" value="Aggregated Trajectory Dataset of Sea Surface Temperatures of SLDMB Buoy %s" />\n        <aggregation dimName="time" type="joinExisting" >\n          <scan location="file:/data/thredds-data/sldmb/%s/" suffix=".nc" />\n        </aggregation>\n      </netcdf>\n    </dataset>' %(buoy,buoy,buoy,buoy,buoy,buoy)
            catalog = "/var/www/tomcats/stable/content/thredds/agg_catalogs/sldmb_agg_catalog.xml"
            line_pattern =  '<!-- INSERT HERE -->\n'
            for line in fileinput.input(catalog, inplace=True):
                if line == line_pattern:
                    print new_buoy_agg
                print line,
            #shutil.copyfile("/opt/sldmb/aggregate_ncml_template.ncml","/data/thredds-data/sldmb/%s/%s_agg.ncml"%(buoy,buoy))
        # Make file
        #if output_filename is None:
        #output_filename = "buoy%s_%s_TO_%s.nc" % (i, starting, ending)
        
        filepath = os.path.join(output_directory, output_filename)
        if os.path.exists(filepath):
            os.unlink(filepath)
        
        #variable_name = full_sensor_urn.split(":")[2]
        variable_name = "temperature"
        
        nc = netCDF4.Dataset(filepath, "w")
        logger.debug("Opened file for writing: %s" % filepath)
    
        nc.setncattr("Conventions", "CF-1.6")
        nc.setncattr("date_created", datetime.utcnow().strftime("%Y-%m-%dT%H:%M:00Z"))
        
        # Station_ID
        nc.createDimension("name_strlen", len(buoy))
        trajectory = nc.createVariable("trajectory_id", "c", ("name_strlen"))
        trajectory.cf_role = "trajectory_id"
        trajectory.standard_name = buoy
        trajectory.long_name = buoy
        trajectory[:] = buoy
        
        
        logger.debug("Setting up time...")
        # Time extents
        nc.setncattr("time_coverage_start",    (starting))
        nc.setncattr("time_coverage_end",     (ending))
        # duration (ISO8601 format)
        nc.setncattr("time_coverage_duration", "P%sS" % unicode(int(round((ending - starting)))))
        # resolution (ISO8601 format)
        # subtract adjacent times to produce an array of differences, then get the most common occurance
    
        if len(data[0]) >1:
            time_diffs = times[1]- times[0]
        else:
            time_diffs = 0
        nc.setncattr("time_coverage_resolution", "P%sS" % unicode(int(round(time_diffs))))
        
        # Time - 32-bit unsigned integer
        nc.createDimension("time")
        time = nc.createVariable("time",    "f8", ( "time",))
        time.units          = "seconds since 1970-01-01T00:00:00Z"
        time.standard_name  = "time"
        time.long_name      = "time of measurement"
        time.calendar       = "gregorian"
        logger.debug("Setting data array...")
        time._CoordianteAxisType = "Time"
        time[:] = times
        
        # Station Location
        #nc.createDimension("lat")
        lat = nc.createVariable("lat", "f4", ("time"))
        lat.units           = "degrees_north"
        lat.standard_name   = "latitude"
        lat.long_name       = "station latitude"
        lat._CoordianteAxisType = "Lat"
        lat[:] = precise_latitude
        
        #nc.createDimension("lon")
        lon = nc.createVariable("lon", "f4", ("time"))
        lon.units           = "degrees_east"
        lon.standard_name   = "longitude"
        lon.long_name       = "station longitude"
        lon._CoordianteAxisType = "Lon"
        lon[:] = precise_longitude
    
        
        # Metadata variables
        crs = nc.createVariable("crs", "i4")
        crs.long_name           = "http://www.opengis.net/def/crs/EPSG/0/4326"
        crs.grid_mapping_name   = "latitude_longitude"
        crs.epsg_code           = "EPSG:4326"
        crs.semi_major_axis     = float(6378137.0)
        crs.inverse_flattening  = float(298.257223563)
        
        platform = nc.createVariable("platform", "i4")
        platform.ioos_code      = full_station_urn
        platform.short_name     = "title", full_station_urn
        platform.long_name      = "description", full_station_urn
        
        instrument = nc.createVariable("instrument", "i4")
        instrument.definition   = "http://mmisw.org/ont/ioos/definition/sensorID"
        instrument.long_name    = full_sensor_urn
        
        # Sync file structure
        nc.sync()
        
        # The coordinates attribute.  This may get appended to below before being written to the sensor variable.
        coordinates = ["time", "lat", "lon", "precise_lat", "precise_lon"]
        
        # This should always be a TimeSeries, but the functionaliety is kept incase that fact changes.
        # TIMESERIES
        nc.setncattr("featureType", "trajectory")
    
        logger.debug("Setting values...")
        var = nc.createVariable(variable_name,    "f4", ( "time" ), fill_value=_FillValue)
            # Set 'coordinates' attribute
        setattr(var, "coordinates", " ".join(coordinates))
        setattr(var, "standard_name", variable_name)
        setattr(var, "units", "K")
        # Set data
        logger.debug("Setting data array...")
        for i in range(len(values)):
            if values[i] == -1000+273.15:
                values[i] = _FillValue
    
        var[:] = values
        nc.project = 'MARACOOS SLDMB'
        nc.title = 'Mid-Atlantic Regional Association Coastal Ocean Observing System Self-Locating Datum Marker Buoy'
        nc.institution = 'MARACOOS'
        nc.keywords = 'MARACOOS, SLDMB, Temperature, SST'
        nc.references = 'http://www.MARACOOS.org'
        nc.geospatial_lat_min = min(precise_latitude)
        nc.geospatial_lat_max = max(precise_latitude)
        nc.geospatial_lon_min = min(precise_longitude)
        nc.geospatial_lon_max = max(precise_longitude)
        nc.publisher_name = 'RPS ASA on behalf of MARACOOS.'
        nc.publisher_phone = '(401) 789-6224'
        nc.publisher_email = 'devops@asascience.com'
        nc.publisher_url = 'http://www.asascience.com/'

    
        nc.close()	
	






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('start_date', nargs='?', help= "Defines the start time desired. Defaults to current date.", default='DEFAULT')
    parser.add_argument('minlon', help= "Defines the Minimum Longitude of Concern  Defaults to -180.", nargs='?', default='-180')
    parser.add_argument('minlat', help= "Defines the Minimum Latitude of Concern  Defaults to -90.", nargs='?', default='-90')
    parser.add_argument('maxlon', help= "Defines the Maximum Longitude of Concern  Defaults to 180.", nargs='?', default='180')
    parser.add_argument('maxlat', help= "Defines the Maximum Latitude of Concern  Defaults to 90.", nargs='?', default='90')
    parser.add_argument('--root_directory', '-rd', help= "Define the root directory for saving files.  Defaults to /data/thredds-data/sldmb", nargs='?', default='/data/thredds-data/sldmb')
    parser.add_argument('--ago', '-a', nargs='?', help= "Defines how many days ago from the 'start_date' field specified the data should be scraped. 1 refers to the last ful day.", default = 1)

    args = parser.parse_args()

    Create(args.start_date, args.minlon, args.minlat, args.maxlon, args.maxlat, args.root_directory, args.ago)
