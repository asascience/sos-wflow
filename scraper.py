#!/usr/bin/env python
'''
scraper.py

Scrapes weatherflow data
'''
from xml.etree import ElementTree as ET
from bs4 import BeautifulSoup
from netCDF4 import Dataset
import dateutil.parser
import requests
import calendar
from datetime import datetime
import numpy as np



sos_ns = '{http://www.opengis.net/sos/1.0}'
om_ns = '{http://www.opengis.net/om/1.0}'
gml_ns = '{http://www.opengis.net/gml/3.2}'
ioos_ns = '{http://www.noaa.gov/ioos/0.6.1}'

class WeatherFlowSOS:
    base_url = 'http://d.weatherflow.com/sos/sos.pl'
    response_format = 'text/xml;schema="ioos/0.6.1"'
    stations = {
     'station-1127': 'urn:x-weatherflow:def:station:weatherflow.com::1127',
     'station-397': 'urn:x-weatherflow:def:station:weatherflow.com::397',
     'station-44554': 'urn:x-weatherflow:def:station:weatherflow.com::44554',
     'station-46182': 'urn:x-weatherflow:def:station:weatherflow.com::46182',
     'station-47091': 'urn:x-weatherflow:def:station:weatherflow.com::47091',
     'station-48242': 'urn:x-weatherflow:def:station:weatherflow.com::48242',
     'station-48328': 'urn:x-weatherflow:def:station:weatherflow.com::48328',
     'station-697': 'urn:x-weatherflow:def:station:weatherflow.com::697',
     'station-700': 'urn:x-weatherflow:def:station:weatherflow.com::700'}
    observations = 'Winds'
    variables = {'WindSpeed':'m/s', 'WindDirection':'deg', 'WindVerticalVelocity':'m/s', 'WindGust':'m/s'}

    def __init__(self):
        self.response = None

    def get_capabilities(self):
        params = {
            'request': 'GetCapabilities',
            'service': 'SOS',
            'version': '1.0.0'
        }
        response = requests.get(self.base_url, params=params)
        return response.text

    def verify_capabilities(self):
        response = self.get_capabilities()
        tree = ET.fromstring(response)

        contents = tree.find(sos_ns + 'Contents')
        observation_offering_list = contents.find(sos_ns + 'ObservationOfferingList')
        if len(observation_offering_list):
            return True
        return False

    def get_observations(self, station, t0, t1):
        station_name = self.stations[station]
    #    http://d.weatherflow.com/sos/sos.pl?eventTime=2013-12-06T12%3A18%3A47Z%2F2013-12-16T12%3A18%3A47Z&
    #        service=SOS&
    #        offering=urn%3Ax-weatherflow%3Adef%3Astation%3Aweatherflow.com%3A%3A700&
    #        request=GetObservation&
    #        version=1.0.0&
    #        responseFormat=text%2Fxml%3Bschema%3D%22ioos%2F0.6.1%22&
    #        observedProperty=Winds
        params = {
            'request' : 'GetObservation',
            'service' : 'SOS',
            'version' : '1.0.0',
            'offering' : station_name,
            'resonseFormat' : self.response_format,
            'observedProperty' : self.observations,
            'eventtime' : t0.isoformat() + '/' + t1.isoformat()
        }
        self.station = station

        response = requests.get(self.base_url, params=params)
        if response.status_code != 200:
            response.raise_for_status()

        self.response = BeautifulSoup(response.text, 'xml')
        obs_list = self.response.find_all(attrs={'name':'Station1NumberOfObservationsTimes'})
        if obs_list is None or len(obs_list) < 1:
            raise Exception('No observations')

        self.num_obs = int(obs_list[0].text)

    def get_latlon(self):
        pos = self.response.find_all('Point', attrs={'gml:id' : 'Station1LatLon'})
        if pos is None or len(pos) < 1:
            raise Exception("Station Position omitted")
        lat, lon = (float(i) for i in pos[0].find('pos').text.split(' '))
        return lat, lon


    def get_variable(self, variable):
        if variable not in self.variables:
            raise KeyError("%s is not a variable of %s" % (variable, self.__class__.__name__))
        uom = self.variables[variable]
        obs = self.response.find_all('Quantity', attrs={'name':variable})
        if obs is None or len(obs) < 1:
            raise Exception('No observations for %s' % variable)

        if len(obs) != self.num_obs:
            raise Exception("%s observations doesn't match the station's number of observations" % variable)

        points = []
        for o in obs:
            if o.attrs['uom'] != uom:
                raise Exception("%s is expecting %s units, received %s" % (variable, uom, o.attrs['uom']))
            if 'xsi:nil' in o.attrs:
                points.append(None)
                continue
            try:
                points.append(float(o.text))
            except ValueError:
                points.append(None)
                continue
        return points

    def get_times(self):
        times = self.response.find_all('timePosition')

        if times is None or len(times) < 1:
            raise Exception('No time positions')

        if len(times) != self.num_obs:
            raise Exception('Number of timesteps does not match number of observations')

        tpoints = []
        for t in times:
            tstr = t.text
            dtg = dateutil.parser.parse(tstr)
            timestamp = calendar.timegm(dtg.utctimetuple())
            tpoints.append(timestamp)
        return tpoints

    def get_time_period(self):
        period = self.response.find('TimePeriod')
        if period is None:
            raise Exception("No time period defined")
        t0 = period.find('beginPosition').text
        t1 = period.find('endPosition').text
        return t0, t1

            
class Scraper:
    def __init__(self, filename, sos):
        self.filename = filename
        self.sos = sos

    def create_nc(self):

        today = datetime.utcnow()
        nc = Dataset(self.filename, 'w')
        nc.project = 'Weatherflow'
        nc.Conventions = 'CF-1.6'
        nc.featureType = 'timeSeries'
        nc.cdm_data_type = 'Station'
        nc.title = 'Weatherflow Observed Buoy Water Velocity'
        nc.institution = 'Weatherflow, Inc.'
        nc.date_created = today.isoformat()
        nc.keywords = 'Weather, Winds'
        nc.references = 'http://www.weatherflow.com'
        lat, lon = self.sos.get_latlon()
        nc.geospatial_lat_min = lat
        nc.geospatial_lat_max = lat
        nc.geospatial_lon_min = lon
        nc.geospatial_lon_max = lon
        nc.publisher_name = 'RPS ASA on behalf of Weatherflow, Inc.'
        nc.publisher_phone = '(401) 789-6224'
        nc.publisher_email = 'devops@asascience.com'
        nc.publisher_url = 'http://www.asascience.com/'

        nc.time_coverage_start, nc.time_coverage_end = self.sos.get_time_period()

        nc.createDimension('time')
        nc.createDimension('name_strlen', len(self.sos.station))
        station_name = nc.createVariable('station_name', 'c', ('name_strlen'))
        station_name.long_name = self.sos.station
        station_name.cf_role = 'timeseries_id'
        station_name[:] = self.sos.station

        lat_var = nc.createVariable('lat', 'f4', tuple())
        lat_var.long_name = 'Latitude'
        lat_var.standard_name = 'latitude'
        lat_var.axis = 'Y'
        lat_var.units = 'degrees_north'

        lon_var = nc.createVariable('lon', 'f4', tuple())
        lon_var.long_name = 'Longitude'
        lon_var.standard_name = 'longitude'
        lon_var.axis = 'X'
        lon_var.units = 'degrees_north'

        time_var = nc.createVariable('time', 'f8', ('time',))
        time_var.long_name = 'Time'
        time_var.standard_name = 'time'
        time_var.units = 'seconds since 1970-01-01'
        time_var.axis = 'T'

        wind_speed = nc.createVariable('wind_speed', 'f4', ('time',), fill_value=-9999.)
        wind_speed.long_name = 'Wind Speed'
        wind_speed.standard_name = 'wind_speed'
        wind_speed.units = 'm/s'
        wind_speed.observation_type = 'measured'
        wind_speed.coordinates = 'time lat lon'

        wind_direction = nc.createVariable('wind_direction', 'f4', ('time',), fill_value=-9999.)
        wind_direction.long_name = 'Wind Direction'
        wind_direction.standard_name = 'wind_direction'
        wind_direction.units = 'deg'
        wind_direction.observation_type = 'measured'
        wind_direction.coordinates = 'time lat lon'

        wind_gust = nc.createVariable('wind_gust', 'f4', ('time',), fill_value=-9999.)
        wind_gust.long_name = 'Wind Gust'
        wind_gust.standard_name = 'wind_speed_of_gust'
        wind_gust.units = 'm/s'
        wind_gust.observation_type = 'measured'
        wind_gust.coordinates = 'time lat lon'


        times = np.array(self.sos.get_times())
        time_var[:] = times

        for i,val in enumerate(self.sos.get_variable('WindSpeed')):
            if val:
                wind_speed[i] = val

        for i,val in enumerate(self.sos.get_variable("WindDirection")):
            if val:
                wind_direction[i] = val

        for i,val in enumerate(self.sos.get_variable('WindGust')):
            if val:
                wind_gust[i] = val

        lat_var[:] = lat
        lon_var[:] = lon

        nc.close()



