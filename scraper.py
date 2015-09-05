#!/usr/bin/env python

# This is a template for a Python scraper on morph.io (https://morph.io)

from __future__ import print_function

# import scraperwiki
# import lxml.html

import re

import requests
import scraperwiki


def main():
    routes = """
      505_LIN_1_2
      98_MNL_1_1
    """.split()

    for route in routes:
        scrape1(route)


def scrape1(route):
    url = "http://maps.travelsouthyorkshire.com/RouteWkt/{}.wkt".format(route)
    r = requests.get(url)

    print(r.text)

    # Some routes have more than one line, for example
    # http://maps.travelsouthyorkshire.com/RouteWkt/70_YTC_2_1.wkt
    for line_index, m in enumerate(re.findall(r'(\([^()]*\))', r.text)):
        print(m)
        line_text = m

        for i,m in enumerate(re.findall(r'(\d+) +(\d+)', line_text)):
            easting, northing = map(int, m)
            print(route, line_index, i, easting, northing)
            scraperwiki.sqlite.save(
              unique_keys=['route', 'line_index', 'point_index'],
              data=dict(route=route,
                line_index=line_index, point_index=i,
                easting=easting, northing=northing))


if __name__ == '__main__':
    main()
