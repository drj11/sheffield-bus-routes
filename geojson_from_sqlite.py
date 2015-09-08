#!/usr/bin/env python3

import itertools
import json
import sqlite3

import latlon

def main():
    connection = sqlite3.connect('scraperwiki.sqlite')
    cursor = connection.cursor()

    cursor.execute(
      "SELECT route, line_index, point_index, easting, northing FROM swdata")

    # :todo: one set of lines per route, not one set for all routes
    geojson_lines = []
    for route, lines in itertools.groupby(cursor, lambda r: r[0]):
        print("============ ROUTE", route)
        for line, points in itertools.groupby(lines, lambda r: r[1]):
            print("========= LINE", line)
            geojson_line = dict(type='LineString', coordinates=[])
            for row in points:
                route, line_index, point_index, easting, northing = row
                p = latlon.osgrid_to_wgs84((easting, northing))
                lat = p.latitude
                lon = p.longitude
                geojson_line['coordinates'].append([lon, lat])
            print(geojson_line)
            geojson_lines.append(geojson_line)

    with open("routes.geojson", 'w') as w:
        w.write("var routes = ")
        json.dump(geojson_lines, w, indent=2)


if __name__ == '__main__':
    main()
