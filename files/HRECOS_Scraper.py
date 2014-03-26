#!/usr/bin/env python

import fileinput
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

import netCDF4
import numpy as np
import zipfile
import StringIO
import logging
logger = logging.getLogger("pytools")
logger.addHandler(logging.NullHandler())



def Create(start_date, tz, units, root_directory, ago):
	#-------------------------------------------------------------------------------------------------------
    sid_list= [
                'HRLCK8H',
                'HRLCK8M',
                'HRALBPH',
                'HRALBPM',
                'HRTVBM',
                'HRMARPH',
                'HRWSTPTH',
                'HRPMNTH',
                'HRPMNTM',       
                'HRPIER84'
                ]

    standard_name_dict = {
                 'Dissolved Oxygen'           : 'mass_concentration_of_oxygen_in_sea_water',
                 'Water Temp'                 : 'sea_water_temperature',
                 'Acidity'                    : 'sea_water_ph_reported_on_total_scale',
                 'Dissolved Oxygen Percent'   : 'fractional_saturation_of_oxygen_in_sea_water',
                 'Turbidity'                  : 'sea_water_turbidity',
                 'Specific Conductivity'      : 'sea_water_electrical_conductivity',
                 'Rainfall'                   : 'thickness_of_rainfall_amount',
                 'Soil Temp'                  : 'soil_temp',
                 'Wind Speed'                 : 'wind_speed',
                 'Air Temp'                   : 'air_temperature',
                 'Wind Gusts'                 : 'wind_speed_of_gust',
                 'Wind Direction'             : 'wind_from_direction',
                 'Relative Humidity'          : 'relative_humidity',
                 'Dew Point'                  : 'dew_point_temperature',
                 'Rainfall Daily Accumulation': 'Rainfall Daily Accumulation',
                 'Radiation (PAR)'            : 'surface_upwelling_photosynthetic_photon_flux_in_air',
                 'Radiation (Total)'          : 'surface_upwelling_photosynthetic_photon_flux_in_air',
                 'Wind Dir Std Dev'           : 'Wind Dir Std Dev',
                 'Barometric Pressure'        : 'air_pressure',
                 'Water Elevation (NAVD88)'   : 'water_surface_height_above_reference_datum',
                 'Depth'                      : 'depth',
                 'Salinity'                   : 'sea_water_salinity'
                 }
 
    attributes_dict = {
                 'Dissolved Oxygen'           : [['comment','These values use the approximation that ppm = mg/l.']],
                 'Water Temp'                 : [],
                 'Acidity'                    : [],
                 'Dissolved Oxygen Percent'   : [],
                 'Turbidity'                  : [],
                 'Specific Conductivity'      : [],
                 'Rainfall'                   : [],
                 'Soil Temp'                  : [],
                 'Wind Speed'                 : [],
                 'Air Temp'                   : [],
                 'Wind Gusts'                 : [],
                 'Wind Direction'             : [],
                 'Relative Humidity'          : [],
                 'Dew Point'                  : [],
                 'Rainfall Daily Accumulation': [['comment','No CF Standard_Name']],
                 'Radiation (PAR)'            : [['comment','This value is the total Radiation since the last sample (15 minute time window']],
                 'Radiation (Total)'          : [['comment','This value is the total Radiation since the last sample (15 minute time window']],
                 'Wind Dir Std Dev'           : [['comment','No CF Standard_Name']],
                 'Barometric Pressure'        : [],
                 'Water Elevation (NAVD88)'   : [['comment','Datum is NAVD88']],
                 'Depth'                      : [],
                 'Salinity'                   : []
                 }
 
    longitude_dict = {
                'HRLCK8H' : -73.9909,
                'HRLCK8M' : -73.9909,
                'HRALBPH' : -73.7581,
                'HRALBPM' : -73.7581,
                'HRTVBM'  : -73.9170,
                'HRMARPH' : -73.9390,
                'HRWSTPTH': -73.9553,
                'HRPMNTH' : -73.8960,
                'HRPMNTM' : -73.8960,       
                'HRPIER84': -74.0031
                }
    latitude_dict = {
                'HRLCK8H' : 42.8290,
                'HRLCK8M' : 42.8290,
                'HRALBPH' : 42.6196,
                'HRALBPM' : 42.6196,
                'HRTVBM'  : 43.0180,
                'HRMARPH' : 41.7207,
                'HRWSTPTH': 41.3860,
                'HRPMNTH' : 41.0431,
                'HRPMNTM' : 41.0431,       
                'HRPIER84': 40.7646
                }
	#Sets start and end date.

    if start_date == 'DEFAULT':
    	start_date = date.today()
    start = start_date-timedelta(int(ago))
    end = str(start)
    start = str(start)

    for sid in sid_list:

        latitude = latitude_dict[sid]
        longitude = longitude_dict[sid]
    
        payload = {
               "sid" : sid,
               "start"  : start, 
               "end"  : end, 
               "tz"  : tz, 
               "units"   : units, 
               }
    
        url = "http://hudson.dl.stevens-tech.edu/cgi-bin/csv_all_dual_query_zip.pl" 
    
        r=requests.get(url, params=payload)
        r.raise_for_status()
    
        First_split = r.content.split("\x00\x00\x00")
        Split = First_split[2].split(".csv")
        full_name = Split[0]
        try:
            print "Data Requested is from: ", full_name
            DATA_FLAG = True
            z = zipfile.ZipFile(StringIO.StringIO(r.content))
            zip_data = z.read(Split[0]+".csv")
            Data_split = zip_data.split("\n")
            i = 0
            dataset = []
            for line in Data_split:
                dataset.append(line)
                dataset[i] = dataset[i].split(",")
                i = i+1
        except:
            DATA_FLAG = False
            print "There is no data for this day for station: %s" %sid
        
        if DATA_FLAG == True:
            data_dict = {}
            data_title_dict = {}
        
            i = 0
            for lines in dataset:
                if len(lines) == 7:
                    fix = lines
                    lines = [fix[0], fix[1] + fix [2], fix[3], fix[4], fix[5], fix[6]]
                if len(lines)>3 and lines[1] != " Param":
                    if not lines[1] in data_dict.keys():
                        value_array = []
                        time_array = []
                        unit_array = []
                        data_dict[lines[1]] = value_array
                        data_title_dict[lines[1]] = lines[1]
                    value_array.append(lines[3])
                    time_array.append(lines[0])
                    unit_array.append(lines[2])
                    data_dict[lines[1]] = value_array
                    data_dict["%s Time" %lines[1]] = time_array
                    data_dict["%s Units" %lines[1]] = unit_array
                i = i+1
        
            station = sid
            output_directory = root_directory+"/%s" %sid
            output_filename = "%s_%s.nc" %(sid, start)
            full_station_urn = "urn:ioos:station:hrecos:%s" %sid
            full_sensor_urn = "urn:ioos:sensor:hrecos:%s" %sid
            global_attributes = []
            attributes = []
            fillvalue = -999.9
    
            time_sort_list = []
            i = 0
            for this in data_title_dict:
                time_sort_list.append(data_dict["%s Time" %this])
            time_sort_list.sort(key=len)
            times = time_sort_list[-1]


            for i in range(len(times)):
                times[i] = clock.mktime(clock.strptime(times[i], "%Y-%m-%d %H:%M:%S"))
            #Make path and copy aggregation template
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
                new_station_agg = '    <dataset name="Station %s Aggregated Data" ID="station%s-agg" urlPath="station%s-agg.ncml" dataType="Point">\n      <metadata inherited="true">\n        <serviceName>all</serviceName>\n        <dataType>Point</dataType>\n      </metadata>\n      <netcdf xmlns="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2">\n        <attribute name="title" value="HRECOS Aggregated Station %s Data" />\n        <attribute name="summary" value="Aggregated Dataset of HRECOS Station %s" />\n        <aggregation dimName="time" type="joinExisting" >\n          <scan location="file:/data/thredds-data/hrecos/%s/" suffix=".nc" />\n        </aggregation>\n      </netcdf>\n    </dataset>' %(station, station, station,station,station,station)
                catalog = "/var/www/tomcats/stable/content/thredds/agg_catalogs/hrecos_agg_catalog.xml"
                line_pattern =  '<!-- INSERT HERE -->\n'
                for line in fileinput.input(catalog, inplace=True):
                    if line == line_pattern:
                        print new_station_agg
                    print line,
                 
    	   # Make file
            if output_filename is None:
                output_filename = "station_%s_%s_TO_%s.nc" % (key, start, end)
                    
                    
            filepath = os.path.join(output_directory, output_filename)
        
            if os.path.exists(filepath):
                os.unlink(filepath)
    	   
            nc = netCDF4.Dataset(filepath, "w")
        
            lat = nc.createVariable("lat", "f4")    	
            lat.units           = "degrees_north"        	
            lat.standard_name   = "latitude"        	
            lat.long_name       = "station latitude"        	
            lat._CoordianteAxisType = "Lat"        	
            lat[:] = latitude      
    
            lon = nc.createVariable("lon", "f4")   	
            lon.units           = "degrees_east"        	
            lon.standard_name   = "longitude"        	
            lon.long_name       = "station longitude"        	
            lon._CoordianteAxisType = "Lon"        	
            lon[:] = longitude
    
    
            nc.createDimension("name_strlen", len(station))        	
            station_name = nc.createVariable("station_id", "c", ("name_strlen"))        	
            station_name.cf_role = "timeseries_id"        	
            station_name.standard_name = sid        	
            station_name.long_name = full_name        	
            station_name[:] = sid
    
            
            nc.createDimension("time")        	
            time = nc.createVariable("time",    "f8", ( "time",))        	
            time.units          = "seconds since 1970-01-01T00:00:00Z"        	
            time.standard_name  = "time"        	
            time.long_name      = "time of measurement"        	
            time.calendar       = "gregorian"        	
            logger.debug("Setting data array...")        	
            time._CoordianteAxisType = "Time"        	
            time[:] = times      
    
            nc.setncattr("time_coverage_start",    (times[0]))        	
            nc.setncattr("time_coverage_end",     (times[-1])) 
            nc.setncattr("time_coverage_duration", "P%sS" % unicode(int(round((times[-1] - times[0])))))        	
            diffs = times[1]-times[0]
            time_diffs = diffs
            nc.setncattr("time_coverage_resolution", "P%sS" % unicode(int(round(time_diffs))))  
            nc.setncattr("Conventions", "CF-1.6")
            nc.setncattr("date_created", datetime.utcnow().strftime("%Y-%m-%dT%H:%M:00Z"))
    	   
            #Metadata variables	
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
    
            nc.project = 'HRECOS'
            nc.title = 'Hudson River Environmental Conditions Observing System'
            nc.institution = 'HRECOS.'
            nc.keywords = 'Hudson, Weather, Wind, Hydro'
            nc.references = 'http://www.HRECOS.org'
            nc.geospatial_lat_min = latitude
            nc.geospatial_lat_max = latitude
            nc.geospatial_lon_min = longitude
            nc.geospatial_lon_max = longitude
            nc.publisher_name = 'RPS ASA on behalf of HRECOS.'
            nc.publisher_phone = '(401) 789-6224'
            nc.publisher_email = 'devops@asascience.com'
            nc.publisher_url = 'http://www.asascience.com/'
    
    
            # Sync file structure           
            nc.sync()
    
    
    
            for key in data_title_dict:
                values = data_dict["%s" %key]
                variable_name = standard_name_dict[key]
                units = data_dict["%s Units" %key][0] 
                           
                #Unit Conversion
                if key == 'Specific Conductivity':
                    for i in range(len(values)):
                        values[i] = float(values[i])/10
                        units = 'S/m'
                        
                if key == 'Water Temp' or key == 'Air Temp' or key == 'Soil Temp' or key == 'Dew Point':
                    for i in range(len(values)):
                        values[i] = float(values[i]) + 273.15
                        units = 'K'
                        
                if key == 'Rainfall' or key == 'Rainfall Daily Accumulation':
                    for i in range(len(values)):
                        values[i] = float(values[i])/1000
                        units = 'm'
                        
                if key == 'Dissolved Oxygen':
                    for i in range(len(values)):
                        values[i] = float(values[i])/1000
                        units = 'kg/m^3'
                        
                if key == 'Barometric Pressure':
                    for i in range(len(values)):
                        values[i] = float(values[i])/100
                        units = 'Pa'       
                	
                coordinates = ["time", "lat", "lon"]
                logger.debug("Opened file for writing: %s" % filepath)
        
                nc.setncattr("featureType", "timeSeries")
    
    	       #logger.debug("Setting values...")	    
                var = nc.createVariable(variable_name,    "f4", ("time" ), fill_value=fillvalue)	   
    	       # Set 'coordinates' attribute	        
                setattr(var, "coordinates", " ".join(coordinates))	        
                setattr(var, "standard_name", variable_name)	        
                setattr(var, "units", units)
                setattr(var, "description", key)
                #setting up extra attributes
                for attribute in attributes_dict[key]:
                    setattr(var,attribute[0],attribute[1])
    	      # Set data	        
                logger.debug("Setting data array...")
                var[:] = values
                nc.sync()
                
            nc.close()	

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('start_date', nargs='?', help= "Defines the start time desired. Defaults to current date.", default='DEFAULT')
	parser.add_argument('tz', help= "Defines the timezone to use.  Defaults to Eastern Standard.", nargs='?', default='ET')
	parser.add_argument('units', help= "Defines the unit system to use.  Defaults to English.", nargs='?', default='english')
	parser.add_argument('--root_directory', '-rd', help= "Define the root directory for saving files.  Defaults to /data/thredds-data/hrecos", nargs='?', default='/data/thredds-data/hrecos')
	parser.add_argument('--ago', '-a', nargs='?', help= "Defines how many days ago from the 'start_date' field specified the data should be scraped. 1 refers to the last ful day.", default = 1)

	args = parser.parse_args()

	Create(args.start_date, args.tz, args.units, args.root_directory, args.ago)#



