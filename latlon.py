#!/usr/bin/env python
# coding=utf-8

from __future__ import division

import math
import sys

from collections import namedtuple

LatLon = namedtuple('LatLon', 'latitude longitude datum')
Point = namedtuple('Point', 'x y z')

# Converted to Python from
# https://github.com/chrisveness/geodesy/blob/master/osgridref.js
def osgrid_to_latlon(gridref):
    datum = 'OSGB36'

    (E, N) = gridref

    # Airy 1830 major & minor semi-axes
    a = 6377563.396
    b = 6356256.909
    # NatGrid scale factor on central meridian
    F0 = 0.9996012717
    # NatGrid true origin
    phi0 = 49*math.pi/180
    lamda0 = -2*math.pi/180
    # northing & easting of true origin, metres
    N0 = -100000
    E0 = 400000

    # eccentricity squared
    e2 = 1 - (b*b)/(a*a)
    # n, n², n³
    n = (a-b)/(a+b)
    n2 = n**2
    n3 = n**3

    phi=phi0
    M=0
    # ie until < 0.01mm
    while N-N0-M >= 0.00001:
        phi = (N-N0-M)/(a*F0) + phi

        Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (phi-phi0)
        Mb = (3*n + 3*n*n + (21/8)*n3) * math.sin(phi-phi0) * math.cos(phi+phi0)
        Mc = ((15/8)*n2 + (15/8)*n3) * math.sin(2*(phi-phi0)) * math.cos(2*(phi+phi0))
        Md = (35/24)*n3 * math.sin(3*(phi-phi0)) * math.cos(3*(phi+phi0))
        # meridional arc
        M = b * F0 * (Ma - Mb + Mc - Md)

    cosphi = math.cos(phi)
    sinphi = math.sin(phi)
    # nu = transverse radius of curvature
    nu = a*F0/math.sqrt(1-e2*sinphi*sinphi)
    # rho = meridional radius of curvature
    rho = a*F0*(1-e2)/math.pow(1-e2*sinphi*sinphi, 1.5)
    # eta = ?
    eta2 = nu/rho-1

    tanphi = math.tan(phi)
    tan2phi = tanphi*tanphi
    tan4phi = tan2phi*tan2phi
    tan6phi = tan4phi*tan2phi
    secphi = 1/cosphi
    nu3 = nu**3
    nu5 = nu3*nu*nu
    nu7 = nu5*nu*nu
    VII = tanphi/(2*rho*nu)
    VIII = tanphi/(24*rho*nu3)*(5+3*tan2phi+eta2-9*tan2phi*eta2)
    IX = tanphi/(720*rho*nu5)*(61+90*tan2phi+45*tan4phi)
    X = secphi/nu
    XI = secphi/(6*nu3)*(nu/rho+2*tan2phi)
    XII = secphi/(120*nu5)*(5+28*tan2phi+24*tan4phi)
    XIIA = secphi/(5040*nu7)*(61+662*tan2phi+1320*tan4phi+720*tan6phi)

    dE = (E-E0)
    dE2 = dE*dE
    dE3 = dE2*dE
    dE4 = dE2*dE2
    dE5 = dE3*dE2
    dE6 = dE4*dE2
    dE7 = dE5*dE2
    phi = phi - VII*dE2 + VIII*dE4 - IX*dE6
    lamda = lamda0 + X*dE - XI*dE3 + XII*dE5 - XIIA*dE7

    point = LatLon(math.degrees(phi), math.degrees(lamda), 'OSGB36')
    if datum != 'OSGB36':
        point = convertDatum(point, datum)

    return point

# Ellipsoids
Ellipsoid = namedtuple('Ellipsoid', 'a b f')
Airy1830 = Ellipsoid(a=6377563.396, b=6356256.909, f=1/299.3249646)
WGS84 = Ellipsoid(a=6378137, b=6356752.31425, f=1/298.257223563)

Datum = namedtuple('Datum', 'ellipsoid transform')
Transform = namedtuple('Transform', 'tx ty tz rx ry rz s')

DATUM = {
    "OSGB36": Datum(Airy1830,
        Transform(tx=-446.448, ty=125.157, tz=-542.060,  # m
                  rx=-0.1502, ry=-0.2470, rz=-0.8421,    # sec
                  s=20.4894                              # ppm
    )),
    "WGS84": Datum(WGS84,
        Transform(tx=0, ty=0, tz=0,
                  rx=0, ry=0, rz=0,
                  s=0
    )),
    }

def convert_datum(latlon, target_datum):
    """
    Convert latlon to use the datum *target_datum*.
    """

    assert 'WGS84' == target_datum
    transform = DATUM[latlon.datum].transform

    cartesian = latlon_to_cartesian(latlon)
    cartesian = apply_inverse_transform(cartesian, transform)
    target_latlon = cartesian_to_latlon(cartesian, target_datum)

    return target_latlon

def apply_inverse_transform(cartesian, transform):
    transform = Transform(*[-v for v in transform])
    return apply_transform(cartesian, transform)

def apply_transform(cartesian, transform):
    x, y, z = cartesian
    tx = transform.tx
    ty = transform.ty
    tz = transform.tz
    rx = math.radians(transform.rx/3600)
    ry = math.radians(transform.ry/3600)
    rz = math.radians(transform.rz/3600)
    s = transform.s * 1e-6 + 1

    # apply
    x2 = tx + x*s - y*rz + z*ry
    y2 = ty + x*rz + y*s - z*rx
    z2 = tz - x*ry + y*rx + z*s

    return Point(x2, y2, z2)


def latlon_to_cartesian(latlon):
    phi = math.radians(latlon.latitude)
    lamda = math.radians(latlon.longitude)
    h = 0 # ellipsoid height, unused

    # :todo: Fetch ellipsoid from datum
    assert latlon.datum == 'OSGB36'
    ellipsoid = DATUM[latlon.datum].ellipsoid
    a = ellipsoid.a
    b = ellipsoid.b

    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    sinlamda = math.sin(lamda)
    coslamda = math.cos(lamda)

    e2 = (a*a - b*b) / (a*a)
    nu = a / math.sqrt(1 - e2*sinphi*sinphi)

    x = (nu+h) * cosphi * coslamda
    y = (nu+h) * cosphi * sinlamda
    z = ((1-e2)*nu + h) * sinphi

    return Point(x, y, z)

def cartesian_to_latlon(cartesian, target_datum):
    assert 'WGS84' == target_datum

    x, y, z = cartesian

    ellipsoid = DATUM[target_datum].ellipsoid
    a = ellipsoid.a
    b = ellipsoid.b

    e2 = (a*a-b*b) / (a*a)       # 1st eccentricity squared
    eps2 = (a*a-b*b) / (b*b)     # 2nd eccentricity squared
    p = math.sqrt(x*x + y*y)     # distance from minor axis
    R = math.sqrt(p*p + z*z)     # polar radius

    # parametric latitude
    # (Bowring eqn 17, replacing tanbeta = z·a / p·b)
    tanbeta = (b*z)/(a*p) * (1+eps2*b/R)
    sinbeta = tanbeta / math.sqrt(1+tanbeta*tanbeta)
    cosbeta = sinbeta / tanbeta

    # geodetic latitude (Bowring eqn 18)
    phi = math.atan2(z + eps2*b*sinbeta*sinbeta*sinbeta,
                       p - e2*a*cosbeta*cosbeta*cosbeta)

    # longitude
    lamda = math.atan2(y, x)

    # height above ellipsoid (Bowring eqn 7) [not currently used]
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    # length of the normal terminated by the minor axis
    v = a*math.sqrt(1-e2*sinphi*sinphi)
    h = p*cosphi + z*sinphi - (a*a/v)

    point = LatLon(
      math.degrees(phi), math.degrees(lamda), target_datum)

    return point


def osgrid_to_wgs84(gridref):
    """
    *gridref* is given as a 2-tuple of (easting, northing) with
    easting and northing both being in metres.

    The equivalent point specified by latitude and longitude on
    the WGS84 system will be returned.
    """

    latlon = osgrid_to_latlon(gridref)
    latlon = convert_datum(latlon, 'WGS84')
    return latlon

def main(argv=None):
    if argv is None:
        argv = sys.argv

    (e, n) = map(float, argv[1:])
    latlon = osgrid_to_wgs84((e, n))
    print(latlon)

if __name__ == '__main__':
    main()
