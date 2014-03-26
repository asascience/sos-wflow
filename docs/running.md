# Required Libraries

        netCDF4
        requests
        python-dateutil
        BeautifulSoup4
        lxml

# Getting the Source

        git clone https://github.com/asascience/sos-wflow/

# Running the Scrapers

        python scraper.py <path-to-weatherflow root> <# hours to scrape> <optional starting time>
        python SLDMB.py <starting-date> <minimum-longitude> <minimum-latitude> <maximum-longitude> <maximum-latitude> <path-to-sldmb root> <amounf-of-days-to-go-back-for-data>
        python HRECOS_Scraper.py <starting-date> <time-zone-to-report-in> <units-for-data> <path-to-hrecos root> <amounf-of-days-to-go-back-for-data>



