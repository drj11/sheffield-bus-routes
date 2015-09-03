# This is a template for a Python scraper on morph.io (https://morph.io)
# including some code snippets below that you should find helpful

from __future__ import print_function

# import scraperwiki
# import lxml.html

import requests


def main():
    r = requests.get('http://maps.travelsouthyorkshire.com/RouteWkt/505_LIN_1_2.wkt')

    print(r.text)

print(__name__)
main()
