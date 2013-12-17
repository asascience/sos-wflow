# Application Layout on DAP01

### Python

Python is installed at the system level, but the version is too old for some of the required libraries, so I installed it to `/opt/python` with version 2.7.6.  I installed all the required libraries at the system-level, mostly because I don't have time to setup a virtual environment.

- To run python `/opt/python/bin/python`

### Tomcat Instances

Apache `httpd` is running on port 80 and acting as a proxy to two tomcat instances, one a production instance and the other a development instance.

#### Production THREDDS

Production THREDDS is running on DAP01. To my knowledge, it doesn't currently have NcSOS installed. It has the default configurations and no new datasets, yet.

- root: `/var/www/tomcats/eds`
- catalog: `/var/www/tomcats/eds/content/thredds/catalog.xml`
- service-script: `service eds-thredds`

#### Dev THREDDS

Development THREDDS has RC7 NcSOS installed. It also has Weatherflow aggregated, static and updating datasets. One thing to note is that the Thredds servlet is running as "dev" instead of "thredds" so the URL paths all start with `/dev/`

- root: `/var/www/tomcats/thredds`
- catalog: `/var/www/tomcats/thredds/content/thredds/catalog.xml`
- service-script: `service dev-thredds`

### Scraper Instance

I haven't been able to write a proper daemon or service script so I just ran the scraper as-is using the tomcat user, with `nohup`. It should run fine as long as the server doesn't get restarted.

- script: `/opt/scraper/scraper.py`
- data directory: `/data/thredds-data/weatherflow`
