# Application Layout on DAP01

### Python

Python is installed at the system level, but the version is too old for some of the required libraries, so I installed it to `/opt/python` with version 2.7.6.  I installed all the required libraries at the system-level, mostly because I don't have time to setup a virtual environment.

- To run python `/opt/python/bin/python`

### Tomcat Instances

Apache `httpd` is running on port 80 and acting as a proxy to three tomcat instances, detailed below.

#### STABLE THREDDS

Stable version of TDS (4.3.18) and latest stable NcSOS (RC8), displaying Weatherflow data.

- root: `/var/www/tomcats/stable`
- catalog: `/var/www/tomcats/stable/content/thredds/catalog.xml`
- service-script: `service thredds-stable`
- web: http://sos.maracoos.org/stable/catalog/catalog.html
- sflow, hrecos, sldmb catalog: see specific subsections below

#### EDS THREDDS

EDS THREDDS is running on DAP01. No configuration yet.

- root: `/var/www/tomcats/eds`
- catalog: `/var/www/tomcats/eds/content/thredds/catalog.xml`
- service-script: `service thredds-eds`
- web: http://sos.maracoos.org/eds/catalog/catalog.html

#### Dev THREDDS (4.4.x)

Development THREDDS has NcSOS 4.4 series branch installed (asascience-open/ncsos@44dd72239b3b365d2081410c9064ec22a48e0dce), displaying Weatherflow data.

- root: `/var/www/tomcats/thredds`
- catalog: `/var/www/tomcats/thredds/content/thredds/catalog.xml`
- service-script: `service thredds-dev`
- web: http://sos.maracoos.org/dev/catalog/catalog.html

### SOS Weatehrflow Scraper Instance

- tomcat user is running scraper hourly via cron
- python: `/opt/python/bin/python`
- script: `/opt/sos-wflow/scraper.py`
- data directory: `/data/thredds-data/weatherflow`
- daily catalog: `/var/www/tomcats/stable/content/thredds/daily_catalogs/weatherflow_catalog.xml`
- agg catalog: `/var/www/tomcats/stable/content/thredds/agg_catalogs/weatherflow_agg_catalog.xml`

### SLDM Weatehrflow Scraper Instance

- tomcat user is running scraper daily via cron
- python: `/opt/python/bin/python`
- script: `/opt/sldmb/SLDMB.py`
- data directory: `/data/thredds-data/sldmb`
- daily catalog: `/var/www/tomcats/stable/content/thredds/daily_catalogs/sldm_catalog.xml`
- agg catalog: `/var/www/tomcats/stable/content/thredds/agg_catalogs/sldm_agg_catalog.xml`
- agg template: `/var/www/tomcats/stable/content/thredds/agg_catalogs/sldm_agg_catalog.nodata`

### HRECOS Weatehrflow Scraper Instance

- tomcat user is running scraper daily via cron
- python: `/opt/python/bin/python`
- script: `/opt/hrecos/HRECOS_Scraper.py`
- data directory: `/data/thredds-data/hrecos`
- daily catalog: `/var/www/tomcats/stable/content/thredds/daily_catalogs/hrecos_catalog.xml`
- agg catalog: `/var/www/tomcats/stable/content/thredds/agg_catalogs/hrecos_agg_catalog.xml`
- agg template: `/var/www/tomcats/stable/content/thredds/agg_catalogs/hrecos_agg_catalog.nodata`



