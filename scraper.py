#!/usr/bin/env python

# This is a template for a Python scraper on morph.io (https://morph.io)

from __future__ import print_function

# import scraperwiki
# import lxml.html

import re

import requests
import scraperwiki


def main():
    # We can't use just a number to identify the route:
    # - routes come in 2 variants ("there" and "back");
    # - same number, different routes
    # - same number, similar route, different provider
    # (for example, search for "70")
    # So the identifier we use should probably be the stuff
    # after the last '/'.
    bus = "505"
    r = requests.get('http://maps.travelsouthyorkshire.com/RouteWkt/505_LIN_1_2.wkt')

    print(r.text)

    for i,m in enumerate(re.findall(r'(\d+) +(\d+)', r.text)):
        easting, northing = map(int, m)
        print(bus, i, easting, northing)
        # Some routes have more than one line, and we might need
        # to model that. For example,
        # http://maps.travelsouthyorkshire.com/RouteWkt/70_YTC_2_1.wkt
        scraperwiki.sqlite.save(unique_keys=['bus', 'point_index'],
          data=dict(bus=bus,
            point_index=i, easting=easting, northing=northing))


if __name__ == '__main__':
    main()
