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

    # One feature per route, each feature being a MultiLineString.
    features = []
    for route, lines in itertools.groupby(cursor, lambda r: r[0]):
        # Will be the coordinate member of a MultiLineString.
        # 3D-list.
        geojson_route = []
        print("============ ROUTE", route)
        for line, points in itertools.groupby(lines, lambda r: r[1]):
            print("========= LINE", line)
            line_string = []
            for row in points:
                route, line_index, point_index, easting, northing = row
                p = latlon.osgrid_to_wgs84((easting, northing))
                lat = p.latitude
                lon = p.longitude
                line_string.append([lon, lat])
            print(line_string)
            geojson_route.append(line_string)
        feature = dict(type="Feature",
          properties=None,
          geometry=dict(type="MultiLineString",
          coordinates=geojson_route))
        features.append(feature)
    geojson_collection = dict(type='FeatureCollection',
      features=features)

    with open("routes.geojson", 'w') as w:
        # w.write("var routes = ")
        json.dump(geojson_collection, w, indent=2)


if __name__ == '__main__':
    main()
